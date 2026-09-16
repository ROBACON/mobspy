"""Tests for error paths in update_model and related methods."""

from __future__ import annotations

import pytest

import mobspy
from mobspy import BaseSpecies, ModelParameters, Simulation
from mobspy.exceptions import SimulationError


class TestUpdateModelErrors:
    """Error paths for Simulation.update_model."""

    def test_update_before_compilation_raises(self) -> None:
        """Calling update_model on an uncompiled simulation must raise."""
        A = BaseSpecies()
        A >> mobspy.Zero @ 1
        A(10)
        S = Simulation(A)
        S.level = -1

        with pytest.raises(SimulationError, match="not compiled yet"):
            S.update_model(["k1", 1])

    def test_malformed_arg_wrong_tuple_length_raises(self) -> None:
        """Args that are not (name, value) pairs must raise."""
        A = BaseSpecies()
        k1 = ModelParameters(1)
        A >> mobspy.Zero @ k1
        A(10)
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)

        with pytest.raises(SimulationError, match="name, value"):
            S.update_model(["k1", 1, "extra"])

    def test_malformed_arg_single_element_raises(self) -> None:
        """A single-element list must also raise."""
        A = BaseSpecies()
        k1 = ModelParameters(1)
        A >> mobspy.Zero @ k1
        A(10)
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)

        with pytest.raises(SimulationError, match="name, value"):
            S.update_model(["k1"])

    def test_unsupported_type_raises(self) -> None:
        """Passing a non-species, non-parameter, non-string key must raise."""
        A = BaseSpecies()
        A >> mobspy.Zero @ 1
        A(10)
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)

        with pytest.raises(SimulationError, match="Unsupported argument type"):
            S.update_model([42, 100])

    def test_string_not_in_parameters_or_species_raises(self) -> None:
        """A string key that matches neither parameter nor species must raise."""
        A = BaseSpecies()
        k1 = ModelParameters(1)
        A >> mobspy.Zero @ k1
        A(10)
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)

        with pytest.raises(
            SimulationError, match="not found either in parameters or species"
        ):
            S.update_model(["nonexistent_name", 42])

    def test_multiple_args_first_malformed_raises(self) -> None:
        """If the first of multiple args is malformed, it raises before processing the rest."""
        A = BaseSpecies()
        k1 = ModelParameters(1)
        A >> mobspy.Zero @ k1
        A(10)
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)

        with pytest.raises(SimulationError, match="name, value"):
            S.update_model(["too", "many", "items"], ["k1", 2])


class TestUpdateModelHappyPaths:
    """Complement existing happy-path coverage for update_model."""

    def test_update_parameter_by_string(self) -> None:
        """Updating a parameter via its string name propagates to SBML."""
        A = BaseSpecies()
        k1 = ModelParameters(0.5)
        A >> mobspy.Zero @ k1
        A(10)
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)

        S.update_model(["k1", 2.0])
        sbml = S.generate_sbml()[0]
        assert '<parameter id="k1" value="2"' in sbml

    def test_update_parameter_by_object(self) -> None:
        """Updating a parameter via the ModelParameters object works."""
        A = BaseSpecies()
        k1 = ModelParameters(0.5)
        A >> mobspy.Zero @ k1
        A(10)
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)

        S.update_model([k1, 3.0])
        sbml = S.generate_sbml()[0]
        assert '<parameter id="k1" value="3"' in sbml

    def test_update_species_count(self) -> None:
        """Updating a species count via the species object works."""
        A = BaseSpecies()
        A >> mobspy.Zero @ 1
        A(10)
        S = Simulation(A)
        S.level = -1
        S.compile(verbose=False)

        S.update_model([A, 50])
        sbml = S.generate_sbml()[0]
        assert 'initialAmount="50"' in sbml


class TestSimRemoveReaction:
    """Tests for sim_remove_reaction utility."""

    def test_remove_reaction_creates_new_sim(self) -> None:
        """Removing a reaction produces a new Simulation without that reaction."""
        A, B = BaseSpecies()
        rxn = A >> B @ 1
        A(10)
        S = Simulation(A | B)

        new_sim = S - rxn

        assert new_sim is not S
        assert rxn not in new_sim._reactions_set

    def test_remove_reaction_via_subtraction(self) -> None:
        """The __sub__ operator delegates to sim_remove_reaction."""
        A, B, C = BaseSpecies()
        rxn1 = A >> B @ 1
        _ = B >> C @ 0.5
        A(10)
        S = Simulation(A | B | C)

        new_sim = S - rxn1

        assert rxn1 not in new_sim._reactions_set

    def test_remove_nonexistent_reaction_raises(self) -> None:
        """Removing a reaction that isn't in the set should raise SimulationError."""
        A, B, C = BaseSpecies()
        _ = A >> B @ 1
        orphan_rxn = B >> C @ 0.5
        A(10)

        S = Simulation(A)

        with pytest.raises(SimulationError):
            S - orphan_rxn
