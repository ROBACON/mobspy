"""Tests for the ConcreteModel backend-agnostic IR."""

from __future__ import annotations

import pytest

from mobspy import BaseSpecies, Simulation
from mobspy.types import (
    CompilerResult,
    ConcreteModel,
    ReactionData,
)


class TestConcreteModelConstruction:
    """Verify ConcreteModel can be built from compilation output."""

    def test_compile_produces_concrete_model(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        S = Simulation(A | B)
        S.compile(verbose=False)
        assert S._concrete_model is not None
        assert isinstance(S._concrete_model, ConcreteModel)

    def test_concrete_model_species(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model
        assert "A" in cm.species
        assert "B" in cm.species
        assert cm.species["A"] == 10
        assert cm.species["B"] == 0

    def test_concrete_model_reactions(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model
        assert len(cm.reactions) == 1
        rxn = next(iter(cm.reactions.values()))
        assert isinstance(rxn, ReactionData)
        assert len(rxn.reactants) == 1
        assert len(rxn.products) == 1

    def test_concrete_model_mappings(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model
        assert "A" in cm.mappings
        assert "B" in cm.mappings

    def test_concrete_model_parameters(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model
        assert "volume" in cm.parameters

    def test_concrete_model_has_mole_false(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        S = Simulation(A | B)
        S.compile(verbose=False)
        assert S._concrete_model.has_mole is False

    def test_concrete_model_model_string(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        S = Simulation(A | B)
        S.compile(verbose=True)
        cm = S._concrete_model
        assert "Species" in cm.model_string
        assert "Reactions" in cm.model_string


class TestConcreteModelBridges:
    """Verify round-trip conversion between ConcreteModel and CompilerResult."""

    def test_to_compiler_result(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model
        cr = cm.to_compiler_result()
        assert isinstance(cr, CompilerResult)
        assert cr.species_for_sbml == cm.species
        assert cr.reactions_for_sbml == cm.reactions
        assert cr.parameters_for_sbml == cm.parameters

    def test_from_compiler_result_roundtrip(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm1 = S._concrete_model
        cr = cm1.to_compiler_result()
        cm2 = ConcreteModel.from_compiler_result(cr)
        assert cm2.species == cm1.species
        assert cm2.reactions == cm1.reactions
        assert cm2.parameters == cm1.parameters
        assert cm2.assigned_species == cm1.assigned_species

    def test_to_compiled_model(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model
        compiled = cm.to_compiled_model(
            species_not_mapped={"A": 10, "B": 0},
            mappings=dict(cm.mappings),
        )
        assert compiled.species_for_sbml == cm.species
        assert compiled.reactions_for_sbml == cm.reactions


class TestConcreteModelWithCharacteristics:
    """Test ConcreteModel with species that have characteristics."""

    def test_species_with_states(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A.alive >> A.dead @ 0.1
        A.alive(100)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model
        assert "A.alive" in cm.species
        assert "A.dead" in cm.species
        assert cm.species["A.alive"] == 100
        assert cm.species["A.dead"] == 0

    def test_inherited_species(self) -> None:
        from mobspy import New

        Bacteria = BaseSpecies(["Bacteria"])
        B1, B2 = New(Bacteria, 2)
        Bacteria.state1 >> Bacteria.state2 @ 1
        B1.state1(50)
        S = Simulation(B1 | B2)
        S.compile(verbose=False)
        cm = S._concrete_model
        assert "B1.state1" in cm.species
        assert "B1.state2" in cm.species
        assert "B2.state1" in cm.species
        assert "B2.state2" in cm.species

    def test_reversible_reaction(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ (1.0, 0.5)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model
        assert len(cm.reactions) == 2


class TestConcreteModelEvents:
    """Test ConcreteModel with events."""

    def test_event_at(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        S = Simulation(A | B)
        S.at(10, {A: 50})
        S.compile(verbose=False)
        cm = S._concrete_model
        assert len(cm.events) > 0


class TestConcreteModelImmutability:
    """Verify ConcreteModel is frozen."""

    def test_frozen(self) -> None:
        cm = ConcreteModel()
        with pytest.raises(AttributeError):
            cm.species = {}  # type: ignore[misc]

    def test_default_values(self) -> None:
        cm = ConcreteModel()
        assert cm.species == {}
        assert cm.reactions == {}
        assert cm.parameters == {}
        assert cm.events == {}
        assert cm.assignments == {}
        assert cm.mappings == {}
        assert cm.assigned_species == ()
        assert cm.model_string == ""
        assert cm.has_mole is False
        assert cm.unit_context is None
