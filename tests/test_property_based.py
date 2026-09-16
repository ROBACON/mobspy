"""Property-based tests for MobsPy compiler invariants.

Uses Hypothesis to generate random model configurations and verify
structural invariants hold for all compiled models.
"""

from __future__ import annotations

from hypothesis import given, settings
from hypothesis import strategies as st

from mobspy import BaseSpecies, New, Simulation

# --- Strategies ---

species_name_st = st.text(
    alphabet=st.sampled_from("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"),
    min_size=1,
    max_size=8,
).filter(lambda s: s not in {"Time", "Rev", "All"})

rate_st = st.floats(
    min_value=0.01, max_value=100.0, allow_nan=False, allow_infinity=False
)

count_st = st.integers(min_value=0, max_value=1000)


# --- Invariant tests ---


class TestCompilationInvariants:
    """Property-based tests for compiler invariants."""

    @given(rate=rate_st, count=count_st)
    @settings(max_examples=20, deadline=5000)
    def test_single_reaction_invariants(self, rate: float, count: int) -> None:
        """A >> B @ rate always produces exactly one reaction and two species."""
        A, B = BaseSpecies(["A", "B"])
        A >> B @ rate
        A(count)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == 2
        assert len(cm.reactions) == 1
        assert cm.species["A"] == count
        assert cm.species["B"] == 0

    @given(rate=rate_st)
    @settings(max_examples=20, deadline=5000)
    def test_reversible_produces_two_reactions(self, rate: float) -> None:
        """A >> B @ (k_fwd, k_rev) always produces exactly two reactions."""
        A, B = BaseSpecies(["A", "B"])
        A >> B @ (rate, rate / 2)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.reactions) == 2

    @given(count=count_st)
    @settings(max_examples=20, deadline=5000)
    def test_count_conservation(self, count: int) -> None:
        """Total initial count equals what was assigned."""
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(count)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model

        total = sum(cm.species.values())
        assert total == count

    @given(n_states=st.integers(min_value=2, max_value=5))
    @settings(max_examples=10, deadline=5000)
    def test_characteristic_expansion_count(self, n_states: int) -> None:
        """A species with n characteristics produces exactly n concrete species."""
        A = BaseSpecies(["A"])
        state_names = [f"s{i}" for i in range(n_states)]
        for name in state_names:
            A.c(name)

        S = Simulation(A)
        S.compile(verbose=False)
        cm = S._concrete_model

        assert len(cm.species) == n_states

    @given(rate=rate_st, count=count_st)
    @settings(max_examples=10, deadline=5000)
    def test_volume_parameter_always_present(self, rate: float, count: int) -> None:
        """Every compiled model has a 'volume' parameter."""
        A, B = BaseSpecies(["A", "B"])
        A >> B @ rate
        A(count)
        S = Simulation(A | B)
        S.compile(verbose=False)
        cm = S._concrete_model

        assert "volume" in cm.parameters

    @given(rate=rate_st)
    @settings(max_examples=10, deadline=5000)
    def test_concrete_model_roundtrip(self, rate: float) -> None:
        """ConcreteModel -> CompilerResult -> ConcreteModel is lossless."""
        from mobspy.types import ConcreteModel

        A, B = BaseSpecies(["A", "B"])
        A >> B @ rate
        S = Simulation(A | B)
        S.compile(verbose=False)

        cm1 = S._concrete_model
        cr = cm1.to_compiler_result()
        cm2 = ConcreteModel.from_compiler_result(cr)

        assert cm2.species == cm1.species
        assert cm2.reactions == cm1.reactions
        assert cm2.parameters == cm1.parameters
        assert cm2.events == cm1.events
        assert cm2.assignments == cm1.assignments
        assert cm2.mappings == cm1.mappings


class TestInheritanceInvariants:
    """Property-based tests for species inheritance."""

    def test_new_inherits_characteristics(self) -> None:
        """New(Parent, n) creates species with parent's characteristics."""
        Parent = BaseSpecies(["Parent"])
        C1, C2, C3 = New(Parent, 3)

        Parent.on
        Parent.off

        S = Simulation(C1 | C2 | C3)
        S.compile(verbose=False)
        cm = S._concrete_model

        for child in (C1, C2, C3):
            name = child.get_name()
            child_species = [k for k in cm.species if k.startswith(name)]
            assert len(child_species) == 2, (
                f"Expected 2 species for {name}, got {len(child_species)}: {child_species}"
            )

    @given(rate=rate_st)
    @settings(max_examples=10, deadline=5000)
    def test_inherited_reaction_expansion(self, rate: float) -> None:
        """Reactions defined on parent expand to all children."""
        Parent = BaseSpecies(["Parent"])
        C1, C2 = New(Parent, 2)
        Parent.a >> Parent.b @ rate
        S = Simulation(C1 | C2)
        S.compile(verbose=False)
        cm = S._concrete_model

        # Each child gets the reaction: 2 children * 1 reaction = 2 reactions
        assert len(cm.reactions) == 2
