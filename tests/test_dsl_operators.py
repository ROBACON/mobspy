"""Tests for the codebase refactor: new syntax, declarations, structured IR."""

from __future__ import annotations

import warnings

import pytest

from mobspy import BaseSpecies, New, Simulation, Zero
from mobspy.execution import generate_sbml_from_compiled
from mobspy.modules.compiler import compile_model as compile_model_direct
from mobspy.modules.declarations import (
    CountAssignment,
    ModelRegistry,
    RatedProduct,
    ReactionDecl,
    get_registry,
    snapshot_and_clear_registry,
)
from mobspy.modules.list_species import List_Species
from mobspy.modules.species_utils import (
    create_orthogonal_vector_structure,
)
from mobspy.types import CharacteristicQuery, ConcreteSpeciesId
from tests.conftest import CompiledModelAssertions

# ---------------------------------------------------------------
# @ rate syntax
# ---------------------------------------------------------------


class TestAtRateSyntax:
    """Tests for the ``@`` rate operator."""

    def test_simple_forward_rate(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.5
        A(10)
        S = Simulation(A | B)
        S.duration = 5
        S.compile(verbose=False)
        reactions = S._reactions_for_sbml
        assert reactions is not None
        assert len(reactions) == 1
        rxn = next(iter(reactions.values()))
        assert rxn.reactants == [(1, "A")]
        assert rxn.products == [(1, "B")]
        assert "1.5" in rxn.kinetics

    def test_reversible_tuple_rate(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ (2.0, 0.5)
        A(10)
        S = Simulation(A | B)
        S.duration = 5
        S.compile(verbose=False)
        reactions = S._reactions_for_sbml
        assert reactions is not None
        assert len(reactions) == 2
        kinetics = sorted(r.kinetics for r in reactions.values())
        assert any("2.0" in k or "2" in k for k in kinetics)
        assert any("0.5" in k for k in kinetics)

    def test_at_rate_with_characteristics(self) -> None:
        A = BaseSpecies(["A"])
        A.alive
        A.dead
        A.alive >> A.dead @ 0.1
        A.alive(100)
        S = Simulation(A)
        S.duration = 5
        S.compile(verbose=False)
        reactions = S._reactions_for_sbml
        assert reactions is not None
        rxn = next(iter(reactions.values()))
        assert rxn.reactants == [(1, "A_dot_alive")]
        assert rxn.products == [(1, "A_dot_dead")]

    def test_at_rate_with_lambda(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A + B >> Zero @ (lambda r1, r2: 0.5 * r1 * r2)
        A(10)
        B(5)
        S = Simulation(A | B)
        S.duration = 5
        S.compile(verbose=False)
        reactions = S._reactions_for_sbml
        assert reactions is not None
        assert len(reactions) == 1

    def test_at_rate_multi_product(self) -> None:
        A, B, C = BaseSpecies(["A", "B", "C"])
        A >> (B + C) @ 0.3
        A(10)
        S = Simulation(A | B | C)
        S.duration = 5
        S.compile(verbose=False)
        reactions = S._reactions_for_sbml
        assert reactions is not None
        assert len(reactions) == 1

    def test_at_rate_multiple_reactions(self) -> None:
        """Multiple reactions with @ syntax in the same model."""
        A, B, C = BaseSpecies(["A", "B", "C"])
        A >> B @ 1.0
        B >> C @ 0.5
        A(10)
        S = Simulation(A | B | C)
        S.duration = 5
        S.compile(verbose=False)
        reactions = S._reactions_for_sbml
        assert reactions is not None
        assert len(reactions) == 2

    def test_rated_product_is_frozen(self) -> None:
        """RatedProduct is a frozen dataclass."""
        A = BaseSpecies(["A"])
        rp = A @ 1.5
        assert isinstance(rp, RatedProduct)
        assert rp.rate == 1.5
        assert rp.is_reversible is False

    def test_rated_product_reversible(self) -> None:
        A = BaseSpecies(["A"])
        rp = A @ (1.0, 0.5)
        assert isinstance(rp, RatedProduct)
        assert rp.is_reversible is True
        assert rp.rate == 1.0
        assert rp.reverse_rate == 0.5

    def test_at_rate_with_inheritance(self) -> None:
        Base = BaseSpecies(["Base"])
        Child = New(Base, ["Child"])
        Base >> Zero @ 0.1
        Base(10)
        Child(5)
        S = Simulation(Base | Child)
        S.duration = 5
        S.compile(verbose=False)
        reactions = S._reactions_for_sbml
        assert reactions is not None
        assert len(reactions) == 2


# ---------------------------------------------------------------
# Model registry
# ---------------------------------------------------------------


class TestModelRegistry:
    """Tests for the ModelRegistry and declaration types."""

    def test_registry_captures_reactions(self) -> None:
        snapshot_and_clear_registry()
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.5
        reg = get_registry()
        assert len(reg.reactions) >= 1
        decl = reg.reactions[-1]
        assert isinstance(decl, ReactionDecl)
        assert decl.rate == 1.5
        snapshot_and_clear_registry()

    def test_registry_captures_counts(self) -> None:
        snapshot_and_clear_registry()
        A = BaseSpecies(["A"])
        A(100)
        reg = get_registry()
        assert len(reg.counts) >= 1
        count = reg.counts[-1]
        assert isinstance(count, CountAssignment)
        assert count.quantity == 100
        snapshot_and_clear_registry()

    def test_snapshot_and_clear(self) -> None:
        snapshot_and_clear_registry()
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        snap = snapshot_and_clear_registry()
        assert len(snap.reactions) >= 1
        assert len(snap.counts) >= 1
        fresh = get_registry()
        assert len(fresh.reactions) == 0
        assert len(fresh.counts) == 0

    def test_simulation_snapshots_registry(self) -> None:
        snapshot_and_clear_registry()
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        S = Simulation(A | B)
        assert hasattr(S, "_declarations")
        assert isinstance(S._declarations, ModelRegistry)
        assert len(S._declarations.reactions) >= 1


# ---------------------------------------------------------------
# Structured IR types
# ---------------------------------------------------------------


class TestConcreteSpeciesId:
    """Tests for ConcreteSpeciesId."""

    def test_simple_species(self) -> None:
        sid = ConcreteSpeciesId("A")
        assert sid.to_sbml_id() == "A"
        assert str(sid) == "A"

    def test_species_with_chars(self) -> None:
        sid = ConcreteSpeciesId("A", ("alive", "red"))
        assert sid.to_sbml_id() == "A_dot_alive_dot_red"
        assert sid.to_display() == "A.alive.red"
        assert str(sid) == "A.alive.red"

    def test_roundtrip(self) -> None:
        original = ConcreteSpeciesId("Ecoli", ("alive", "motile"))
        sbml_str = original.to_sbml_id()
        parsed = ConcreteSpeciesId.from_sbml_id(sbml_str)
        assert parsed == original

    def test_frozen(self) -> None:
        sid = ConcreteSpeciesId("A", ("x",))
        with pytest.raises(AttributeError):
            sid.base = "B"  # type: ignore[misc]

    def test_custom_separator(self) -> None:
        sid = ConcreteSpeciesId("A", ("x", "y"))
        assert sid.to_sbml_id(separator="__") == "A__x__y"


class TestCharacteristicQuery:
    """Tests for CharacteristicQuery."""

    def test_from_std_char(self) -> None:
        q = CharacteristicQuery.from_legacy_set("std$")
        assert q.include == frozenset()
        assert q.match_all is False

    def test_from_all_char(self) -> None:
        q = CharacteristicQuery.from_legacy_set("all$")
        assert q.match_all is True

    def test_from_set(self) -> None:
        q = CharacteristicQuery.from_legacy_set({"alive", "red"})
        assert q.include == frozenset({"alive", "red"})
        assert q.match_all is False

    def test_from_set_with_all(self) -> None:
        q = CharacteristicQuery.from_legacy_set({"alive", "all$"})
        assert q.include == frozenset({"alive"})
        assert q.match_all is True


# ---------------------------------------------------------------
# compile_model standalone function
# ---------------------------------------------------------------


class TestCompileModel:
    """Tests for the standalone compile_model function."""

    def test_basic_compilation(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        model = List_Species([A, B])
        ovs = create_orthogonal_vector_structure(model)
        species_counts = [
            {
                "object": spe,
                "characteristics": c["characteristics"],
                "quantity": c["quantity"],
            }
            for spe in model
            for c in spe.get_quantities()
        ]
        reactions_set: set[object] = set()
        for spe in model:
            for ref in spe.get_references():
                reactions_set = reactions_set.union(ref.get_reactions())
        result = compile_model_direct(
            model,
            reactions_set=reactions_set,
            species_counts=species_counts,
            orthogonal_vector_structure=ovs,
            verbose=False,
        )
        assert "A" in result.species_for_sbml
        assert "B" in result.species_for_sbml
        assert result.species_for_sbml["A"] == 10
        assert len(result.reactions_for_sbml) == 1

    def test_sbml_generation_from_compiled(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.0
        A(10)
        S = Simulation(A | B)
        S.duration = 5
        S.compile(verbose=False)
        model = List_Species([A, B])
        ovs = create_orthogonal_vector_structure(model)
        species_counts = [
            {
                "object": spe,
                "characteristics": c["characteristics"],
                "quantity": c["quantity"],
            }
            for spe in model
            for c in spe.get_quantities()
        ]
        reactions_set: set[object] = set()
        for spe in model:
            for ref in spe.get_references():
                reactions_set = reactions_set.union(ref.get_reactions())
        result = compile_model_direct(
            model,
            reactions_set=reactions_set,
            species_counts=species_counts,
            orthogonal_vector_structure=ovs,
            verbose=False,
        )
        sbml = generate_sbml_from_compiled(result)
        assert "<?xml" in sbml
        assert "species" in sbml.lower()


# ---------------------------------------------------------------
# Explicit event API
# ---------------------------------------------------------------


class TestExplicitEventAPI:
    """Tests for S.at() and S.when()."""

    def test_at_basic(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 0.1
        A(100)
        S = Simulation(A | B)
        S.duration = 50
        S.at(10, {A: 50})
        S.compile(verbose=False)
        events = S._events_for_sbml
        assert events is not None
        assert len(events) >= 1
        event = next(iter(events.values()))
        assert event.trigger == "true"

    def test_when_basic(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 0.1
        A(100)
        S = Simulation(A | B)
        S.duration = 50
        S.when(A <= 20, {B: 200})
        S.compile(verbose=False)
        events = S._events_for_sbml
        assert events is not None
        assert len(events) >= 1

    def test_at_with_characteristics(self) -> None:
        A = BaseSpecies(["A"])
        A.alive
        A.dead
        A.alive >> A.dead @ 0.1
        A.alive(100)
        S = Simulation(A)
        S.duration = 50
        S.at(10, {A.alive: 50, A.dead: 0})
        S.compile(verbose=False)
        events = S._events_for_sbml
        assert events is not None
        event = next(iter(events.values()))
        assignments = event.assignments
        species_names = [a[0] for a in assignments]
        assert "A_dot_alive" in species_names
        assert "A_dot_dead" in species_names

    def test_when_with_delay(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 0.1
        A(100)
        S = Simulation(A | B)
        S.duration = 50
        S.when(A <= 10, {B: 500}, delay=5)
        S.compile(verbose=False)
        events = S._events_for_sbml
        assert events is not None

    def test_when_rejects_bool(self) -> None:
        from mobspy.exceptions import ValidationError

        A, B = BaseSpecies(["A", "B"])
        A >> B @ 0.1
        A(100)
        S = Simulation(A | B)
        with pytest.raises(ValidationError):
            S.when(True, {A: 50})  # type: ignore[arg-type]

    def test_at_coexists_with_context_manager(self) -> None:
        """Both S.at() and context manager events work together."""
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 0.1
        A(100)
        S = Simulation(A | B)
        S.duration = 50
        S.at(10, {A: 50})
        with S.event_time(20):
            B(200)
        S.compile(verbose=False)
        events = S._events_for_sbml
        assert events is not None
        assert len(events) >= 2


# ---------------------------------------------------------------
# Deprecation warnings
# ---------------------------------------------------------------


class TestDeprecationWarnings:
    """Tests for deprecation warnings on old patterns."""

    def test_rev_emits_deprecation(self) -> None:
        from mobspy import Rev

        A, B = BaseSpecies(["A", "B"])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            Rev[A >> B][1.0, 0.5]
            deprecation_warnings = [
                x for x in w if issubclass(x.category, DeprecationWarning)
            ]
            assert len(deprecation_warnings) >= 1
            assert "deprecated" in str(deprecation_warnings[0].message).lower()


# ---------------------------------------------------------------
# Semantic test helpers
# ---------------------------------------------------------------


class TestSemanticAssertions:
    """Tests using CompiledModelAssertions instead of golden files."""

    def test_simple_model_structure(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 1.5
        A(10)
        S = Simulation(A | B)
        S.duration = 5
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_species("A", "B")
        m.species_count("A", 10)
        m.species_count("B", 0)
        m.has_n_reactions(1)
        m.has_reaction_involving("A")
        m.kinetics_contains("1.5")

    def test_reversible_structure(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ (2.0, 0.5)
        A(10)
        S = Simulation(A | B)
        S.duration = 5
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_n_reactions(2)
        m.kinetics_contains("2")
        m.kinetics_contains("0.5")

    def test_inheritance_expansion(self) -> None:
        Base = BaseSpecies(["Base"])
        Child1 = New(Base, ["Child1"])
        Child2 = New(Base, ["Child2"])
        Base >> Zero @ 0.1
        Base(10)
        Child1(5)
        Child2(3)
        S = Simulation(Base | Child1 | Child2)
        S.duration = 5
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_species("Base", "Child1", "Child2")
        m.has_n_reactions(3)

    def test_characteristics_expand_species(self) -> None:
        A = BaseSpecies(["A"])
        A.alive
        A.dead
        A.alive >> A.dead @ 0.1
        A.alive(100)
        S = Simulation(A)
        S.duration = 5
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_species("A.alive", "A.dead")
        m.species_count("A.alive", 100)
        m.has_n_reactions(1)

    def test_event_adds_event_data(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        A >> B @ 0.1
        A(100)
        S = Simulation(A | B)
        S.duration = 50
        S.at(10, {A: 50})
        S.compile(verbose=False)
        assert S._events_for_sbml is not None
        assert len(S._events_for_sbml) >= 1


# ---------------------------------------------------------------
# Rate expression builder
# ---------------------------------------------------------------


class TestRateBuilder:
    """Tests for the RateExpression builder API."""

    def test_basic_expression(self) -> None:
        from mobspy.modules.rate_builder import param_ref, species_ref

        rate = param_ref("k") * species_ref("A")
        assert str(rate) == "(k*A)"

    def test_compound_expression(self) -> None:
        from mobspy.modules.rate_builder import species_ref

        rate = 0.5 * species_ref("A") / (1 + species_ref("A"))
        assert "0.5" in str(rate)
        assert "A" in str(rate)

    def test_power_expression(self) -> None:
        from mobspy.modules.rate_builder import species_ref

        rate = species_ref("A") ** 2
        assert str(rate) == "(A^2)"

    def test_hill_function(self) -> None:
        from mobspy.modules.rate_builder import hill

        h = hill("S", "Vmax", "Km", n=2)
        rendered = str(h)
        assert "Vmax" in rendered
        assert "Km" in rendered
        assert "S" in rendered

    def test_builder_rate_in_reaction(self) -> None:
        from mobspy.modules.rate_builder import species_ref

        A, B = BaseSpecies(["A", "B"])
        rate = 0.5 * species_ref(A)
        A >> B @ rate
        A(10)
        S = Simulation(A | B)
        S.duration = 5
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_n_reactions(1)
        m.kinetics_contains("0.5")
        m.kinetics_contains("A")

    def test_builder_with_param_ref(self) -> None:
        from mobspy.modules.rate_builder import param_ref, species_ref

        A, B = BaseSpecies(["A", "B"])
        rate = param_ref("k") * species_ref(A)
        A >> B @ rate
        A(10)
        S = Simulation(A | B)
        S.duration = 5
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_n_reactions(1)
        m.kinetics_contains("k")

    def test_where_conditional(self) -> None:
        from mobspy.modules.rate_builder import where

        rate = where("A > 5", 0.5, 1.0)
        assert "piecewise" in str(rate)
        assert "0.5" in str(rate)
        assert "1.0" in str(rate)

    def test_where_method(self) -> None:
        from mobspy.modules.rate_builder import species_ref

        rate = species_ref("A").where("A > 5", 1.0)
        assert "piecewise" in str(rate)

    def test_conditional_in_reaction(self) -> None:
        from mobspy.modules.rate_builder import where

        A, B = BaseSpecies(["A", "B"])
        rate = where("A > 5", 0.5, 1.0)
        A >> B @ rate
        A(10)
        S = Simulation(A | B)
        S.duration = 5
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_n_reactions(1)
        m.kinetics_contains("piecewise")


# ---------------------------------------------------------------
# .named() and [] / @ equivalence
# ---------------------------------------------------------------


class TestBracketAtEquivalence:
    """Tests that [] and @ produce the same result, and [] is deprecated."""

    def test_bracket_emits_deprecation(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            A >> B[1.0]
            dep = [x for x in w if issubclass(x.category, DeprecationWarning)]
            assert len(dep) >= 1
            assert "@" in str(dep[0].message)

    def test_bracket_still_works_with_stoichiometry(self) -> None:
        """[] supports 2*A >> 3*B[rate] which @ cannot."""
        A, B = BaseSpecies(["A", "B"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            2 * A >> 3 * B[0.5]
        A(10)
        S = Simulation(A | B)
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_n_reactions(1)

    def test_bracket_with_lambda(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            A >> B[lambda r1: 0.5 * r1]
        A(10)
        S = Simulation(A | B)
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_n_reactions(1)

    def test_bracket_with_characteristics(self) -> None:
        A = BaseSpecies(["A"])
        A.alive
        A.dead
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            A.alive >> A.dead[0.1]
        A.alive(100)
        S = Simulation(A)
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_n_reactions(1)

    def test_bracket_and_at_same_result(self) -> None:
        A1, B1 = BaseSpecies(["A1", "B1"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            A1 >> B1[1.5]
        A1(10)
        S1 = Simulation(A1 | B1)
        S1.compile(verbose=False)

        A2, B2 = BaseSpecies(["A2", "B2"])
        A2 >> B2 @ 1.5
        A2(10)
        S2 = Simulation(A2 | B2)
        S2.compile(verbose=False)

        r1 = next(iter(S1._reactions_for_sbml.values()))
        r2 = next(iter(S2._reactions_for_sbml.values()))
        assert r1.kinetics == r2.kinetics.replace("A2", "A1").replace("B2", "B1")


class TestNamedMethod:
    """Tests for Species.named() explicit naming."""

    def test_named_sets_name(self) -> None:
        Color, Size = BaseSpecies(["Color", "Size"])
        Color.red
        Color.blue
        Size.big
        Size.small
        Thing = (Color * Size).named("MyThing")
        assert Thing.get_name() == "MyThing"

    def test_named_in_simulation(self) -> None:
        Color, Size = BaseSpecies(["Color", "Size"])
        Color.red
        Color.blue
        Size.big
        Size.small
        Thing = (Color * Size).named("Widget")
        Thing(10)
        Color.red >> Color.blue @ 0.1
        S = Simulation(Thing)
        S.compile(verbose=False)
        m = CompiledModelAssertions(S)
        m.has_species("Widget.big.blue")
        m.has_n_reactions(2)
