"""Semantic compilation tests: assert on model structure, not string output.

These tests complement the golden-file tests in test_compilation.py.
They assert on species counts, reaction counts, kinetics patterns,
and structural properties, making them resilient to formatting changes.
"""

from __future__ import annotations

import pytest

import mobspy
from mobspy import All, BaseSpecies, New, Simulation
from tests.conftest import CompiledModelAssertions


@pytest.mark.compilation
class TestBasicModelStructure:
    """Semantic versions of TestBasicModels golden-file tests."""

    def test_simple_bimolecular(self) -> None:
        """model_1: A + B >> C[1]."""
        A, B, C = BaseSpecies(["A", "B", "C"])
        A + B >> C @ 1
        S = Simulation(A | B | C)
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_species("A", "B", "C")
        m.has_n_reactions(1)
        m.kinetics_contains("1")

    def test_inheritance_expansion(self) -> None:
        """model_4 pattern: parent reaction expands to children."""
        Bacteria, Virus = BaseSpecies(["Bacteria", "Virus"])
        B1, B2 = New(Bacteria, ["B1", "B2"])
        V1, V2 = New(Virus, ["V1", "V2"])
        Bacteria.not_infected + Virus >> Bacteria.infected @ 1
        S = Simulation(B1 | B2 | V1 | V2)
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_species(
            "B1_dot_not_infected",
            "B1_dot_infected",
            "B2_dot_not_infected",
            "B2_dot_infected",
            "V1",
            "V2",
        )
        # 2 bacteria * 2 viruses = 4 reactions
        m.has_n_reactions(4)

    def test_self_replication(self) -> None:
        """model_5 pattern: A >> 2*A and 2*A >> 3*A with inheritance."""
        A = BaseSpecies(["A"])
        B, C = New(A, ["B", "C"])
        A >> 2 * A @ 1
        2 * A >> 3 * A @ 1
        S = Simulation(B | C)
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_species("B", "C")
        # A>>2A expands to 2 (B,C), 2A>>3A expands to 4 (BB,BC,CB,CC) = 6
        m.has_n_reactions(6)

    def test_zero_production(self) -> None:
        """model_6 pattern: Zero >> 2*A with characteristics."""
        A = BaseSpecies(["A"])
        B = New(A, ["B"])
        C = New(A, ["C"])
        B.b1
        B.b2
        C.c1
        C.c2
        mobspy.Zero >> 2 * A @ 1
        S = Simulation(B | C)
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_species("B_dot_b1", "B_dot_b2", "C_dot_c1", "C_dot_c2")
        # Born species: 2 children * 2 chars each = 4 default products
        m.has_n_reactions(4)


@pytest.mark.compilation
class TestCharacteristicStructure:
    """Tests verifying characteristic expansion semantics."""

    def test_single_characteristic_axis(self) -> None:
        A = BaseSpecies(["A"])
        A.alive
        A.dead
        A.alive >> A.dead @ 0.5
        A.alive(100)
        S = Simulation(A)
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_species("A_dot_alive", "A_dot_dead")
        m.species_count("A_dot_alive", 100)
        m.species_count("A_dot_dead", 0)
        m.has_n_reactions(1)

    def test_multi_axis_species(self) -> None:
        Color, Size = BaseSpecies(["Color", "Size"])
        Color.red
        Color.blue
        Size.big
        Size.small
        Thing = Color * Size
        Thing(10)
        Color.red >> Color.blue @ 0.1
        S = Simulation(Thing)
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        # 2 colors * 2 sizes = 4 concrete species (order: Size then Color)
        m.has_species(
            "Thing_dot_big_dot_blue",
            "Thing_dot_big_dot_red",
            "Thing_dot_small_dot_blue",
            "Thing_dot_small_dot_red",
        )
        # red >> blue expands to 2 reactions (big and small variants)
        m.has_n_reactions(2)

    def test_all_operator_counts(self) -> None:
        A = BaseSpecies(["A"])
        A.x
        A.y
        A.z
        All[A](10)
        S = Simulation(A)
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.species_count("A_dot_x", 10)
        m.species_count("A_dot_y", 10)
        m.species_count("A_dot_z", 10)


@pytest.mark.compilation
class TestAtSyntaxEquivalence:
    """Verify @ syntax produces consistent structure across models."""

    def test_at_syntax_consistency(self) -> None:
        A1, B1 = BaseSpecies(["A1", "B1"])
        A1 >> B1 @ 2.5
        A1(10)
        S1 = Simulation(A1 | B1)
        S1.compile(verbose=False)

        A2, B2 = BaseSpecies(["A2", "B2"])
        A2 >> B2 @ 2.5
        A2(10)
        S2 = Simulation(A2 | B2)
        S2.compile(verbose=False)

        # Same structure
        assert len(S1._reactions_for_sbml) == len(S2._reactions_for_sbml)
        r1 = next(iter(S1._reactions_for_sbml.values()))
        r2 = next(iter(S2._reactions_for_sbml.values()))
        assert r1.kinetics == r2.kinetics.replace("A2", "A1").replace("B2", "B1")

    def test_at_reversible_vs_rev(self) -> None:
        A1, B1 = BaseSpecies(["A1", "B1"])
        A1 >> B1 @ (1.0, 0.5)
        A1(10)
        S1 = Simulation(A1 | B1)
        S1.compile(verbose=False)
        m1 = CompiledModelAssertions(S1)
        m1.has_n_reactions(2)

        A2, B2 = BaseSpecies(["A2", "B2"])
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            mobspy.Rev[A2 >> B2][1.0, 0.5]
        A2(10)
        S2 = Simulation(A2 | B2)
        S2.compile(verbose=False)
        m2 = CompiledModelAssertions(S2)
        m2.has_n_reactions(2)


@pytest.mark.compilation
class TestEventStructure:
    """Semantic tests for event compilation."""

    def test_time_event_structure(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 0.1
        A(100)
        S = Simulation(A | B)
        S.duration = 50
        with S.event_time(10):
            A(50)
        S.compile(verbose=False)
        assert S._events_for_sbml is not None
        assert len(S._events_for_sbml) == 1
        event = next(iter(S._events_for_sbml.values()))
        assert event.trigger == "true"

    def test_explicit_at_structure(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 0.1
        A(100)
        S = Simulation(A | B)
        S.duration = 50
        S.at(10, {A: 50})
        S.compile(verbose=False)
        assert S._events_for_sbml is not None
        assert len(S._events_for_sbml) == 1
        event = next(iter(S._events_for_sbml.values()))
        assert event.trigger == "true"

    def test_condition_event_structure(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 0.1
        A(100)
        S = Simulation(A | B)
        S.duration = 50
        S.when(A <= 20, {B: 200})
        S.compile(verbose=False)
        assert S._events_for_sbml is not None
        assert len(S._events_for_sbml) == 1
