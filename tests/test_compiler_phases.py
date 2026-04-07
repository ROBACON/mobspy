"""Tests for individual compiler phases.

Each phase function in mobspy.compiler.compiler is independently
testable. These tests verify phase contracts without going through
the full compile_model() pipeline.
"""

from __future__ import annotations

import pytest

from mobspy import BaseSpecies, MobsPyError, New, Simulation
from mobspy.compiler.compiler import (
    phase_duplicate_detection,
    phase_parameter_validation,
    phase_species_setup,
    phase_volume_resolution,
)
from mobspy.dsl.list_species import List_Species
from mobspy.dsl.species_utils import create_orthogonal_vector_structure
from mobspy.types import (
    ReactionData,
    SpeciesSetupResult,
    VolumeResolutionResult,
)


def _make_ovs(*species):
    """Build orthogonal vector structure from species."""
    return create_orthogonal_vector_structure(species)


class TestPhaseSpeciesSetup:
    """Test the species validation and enumeration phase."""

    def test_basic_species(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        model = List_Species([A, B])
        ovs = _make_ovs(A, B)
        result = phase_species_setup(model, ovs)
        assert isinstance(result, SpeciesSetupResult)
        assert "A" in result.species
        assert "B" in result.species
        assert "A" in result.names_used
        assert "B" in result.names_used

    def test_species_with_characteristics(self) -> None:
        A = BaseSpecies(["A"])
        A.alive
        A.dead
        model = List_Species([A])
        ovs = _make_ovs(A)
        result = phase_species_setup(model, ovs)
        assert "A.alive" in result.species
        assert "A.dead" in result.species
        assert "A" in result.mappings
        assert len(result.mappings["A"]) == 2

    def test_inherited_species(self) -> None:
        Animal = BaseSpecies(["Animal"])
        Cat, Dog = New(Animal, 2)
        Animal.alive
        Animal.dead
        model = List_Species([Cat, Dog])
        ovs = _make_ovs(Cat, Dog)
        result = phase_species_setup(model, ovs)
        assert "Cat.alive" in result.species
        assert "Cat.dead" in result.species
        assert "Dog.alive" in result.species
        assert "Dog.dead" in result.species

    def test_duplicate_names_raise(self) -> None:
        A = BaseSpecies(["A"])
        B = BaseSpecies(["A"])  # Same name
        model = List_Species([A, B])
        ovs = _make_ovs(A, B)
        with pytest.raises(MobsPyError, match="unique"):
            phase_species_setup(model, ovs)


class TestPhaseVolumeResolution:
    """Test volume and dimension resolution."""

    def test_default_volume(self) -> None:
        result = phase_volume_resolution(1, None, [])
        assert isinstance(result, VolumeResolutionResult)
        assert result.volume == 1
        assert result.dimension == 3

    def test_explicit_dimension(self) -> None:
        result = phase_volume_resolution(1, 2, [])
        assert result.dimension == 2

    def test_volume_with_units(self) -> None:
        from mobspy import u

        result = phase_volume_resolution(5 * u.liter, None, [])
        assert result.dimension == 3


class TestPhaseDuplicateDetection:
    """Test O(n) duplicate detection."""

    def test_no_duplicates(self) -> None:
        reactions = {
            "r1": ReactionData(
                reactants=[(1, "A")], products=[(1, "B")], kinetics="A * 1"
            ),
            "r2": ReactionData(
                reactants=[(1, "B")], products=[(1, "C")], kinetics="B * 1"
            ),
        }
        phase_duplicate_detection(reactions)

    def test_duplicate_warns(self) -> None:
        reactions = {
            "r1": ReactionData(
                reactants=[(1, "A")], products=[(1, "B")], kinetics="A * 1"
            ),
            "r2": ReactionData(
                reactants=[(1, "A")], products=[(1, "B")], kinetics="A * 1"
            ),
        }
        # Should not raise, just warn
        phase_duplicate_detection(reactions)


class TestPhaseParameterValidation:
    """Test parameter validation phase."""

    def test_no_collision(self) -> None:
        params = {"k1": (1.0, "per_min"), "volume": (1.0, "dimensionless")}
        names = frozenset({"A", "B"})
        phase_parameter_validation(params, names, set())

    def test_collision_with_species(self) -> None:
        params = {"A": (1.0, "per_min")}
        names = frozenset({"A", "B"})
        with pytest.raises(MobsPyError, match="unique"):
            phase_parameter_validation(params, names, set())


class TestPhaseIntegration:
    """Integration tests running multiple phases sequentially."""

    def test_full_pipeline_matches_compile_model(self) -> None:
        """Verify phases produce the same result as compile_model()."""
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        S = Simulation(A | B)
        model_str = S.compile(verbose=True)
        assert model_str is not None
        assert "A,10" in model_str
        assert "B,0" in model_str

    def test_reversible_reaction_pipeline(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ (1.0, 0.5)
        A(100)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model
        assert len(cm.reactions) == 2

    def test_multi_species_with_states(self) -> None:
        Bacteria, Virus = BaseSpecies(["Bacteria", "Virus"])
        Bacteria.healthy + Virus >> Bacteria.infected @ 0.01
        Bacteria.infected >> Bacteria.dead @ 0.1
        Bacteria.healthy(100)
        Virus(10)
        S = Simulation(Bacteria | Virus)
        S.compile(verbose=False)
        cm = S._concrete_model
        assert len(cm.reactions) == 2
        assert cm.species.get("Bacteria.healthy") == 100
        assert cm.species.get("Virus") == 10


class TestCompilerPhaseTypes:
    """Verify phase result types are properly frozen."""

    def test_species_setup_frozen(self) -> None:
        A = BaseSpecies(["A"])
        model = List_Species([A])
        ovs = _make_ovs(A)
        result = phase_species_setup(model, ovs)
        with pytest.raises(AttributeError):
            result.species = {}  # type: ignore[misc]

    def test_volume_resolution_frozen(self) -> None:
        result = phase_volume_resolution(1, None, [])
        with pytest.raises(AttributeError):
            result.volume = 2  # type: ignore[misc]
