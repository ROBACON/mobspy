"""Unit tests for compiler_operators: NOT-expansion logic."""

from __future__ import annotations

from typing import Any

import pytest

from mobspy import BaseSpecies
from mobspy.constants import NOT_CHAR
from mobspy.modules.compiler_operators import (
    create_all_not_reactions,
    get_all_non_listed_characteristics,
    new_reaction_with_new_characteristics,
)
from mobspy.modules.reactions import Reacting_Species, Reactions

# ---------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------


def _build_reaction(
    reactant_species: Any,
    reactant_chars: set[str],
    product_species: Any,
    product_chars: set[str],
    rate: float = 1.0,
) -> Reactions:
    """Build a Reactions object from species and characteristics."""
    react = Reacting_Species(reactant_species, reactant_chars)
    prod = Reacting_Species(product_species, product_chars)
    return react >> prod[rate]


def _reactant_species_names(reaction: Reactions) -> list[tuple[Any, set[str]]]:
    """Extract (species, characteristics) from reactants."""
    return [(r["object"], r["characteristics"]) for r in reaction.reactants]


def _product_species_names(reaction: Reactions) -> list[tuple[Any, set[str]]]:
    """Extract (species, characteristics) from products."""
    return [(p["object"], p["characteristics"]) for p in reaction.products]


# ---------------------------------------------------------------
# get_all_non_listed_characteristics
# ---------------------------------------------------------------


class TestGetAllNonListedCharacteristics:
    """Tests for get_all_non_listed_characteristics."""

    def test_single_dimension_excludes_listed(self) -> None:
        """With one dimension, excluding one yields the rest."""
        (A,) = BaseSpecies(["A"])
        A.red
        A.blue
        A.green

        result = get_all_non_listed_characteristics(A, {"red"})
        flat = {frozenset(s) for s in result}
        assert frozenset({"blue"}) in flat
        assert frozenset({"green"}) in flat
        assert len(flat) == 2

    def test_excludes_multiple_characteristics(self) -> None:
        """Excluding two from three yields one."""
        (A,) = BaseSpecies(["A"])
        A.red
        A.blue
        A.green

        result = get_all_non_listed_characteristics(A, {"red", "blue"})
        assert len(result) == 1
        assert "green" in result[0]

    def test_empty_characteristics_yields_empty(self) -> None:
        """Excluding empty set yields empty dimension (no characteristics added)."""
        (A,) = BaseSpecies(["A"])
        A.x
        A.y

        result = get_all_non_listed_characteristics(A, set())
        assert len(result) == 1
        assert result[0] == set()

    def test_species_with_no_characteristics(self) -> None:
        """A species with no characteristics yields a single empty set."""
        (A,) = BaseSpecies(["A"])

        result = get_all_non_listed_characteristics(A, {"anything"})
        assert len(result) == 1
        assert result[0] == set()

    def test_preserves_operator_characteristics(self) -> None:
        """Operator markers (containing $) are preserved in every combination."""
        (A,) = BaseSpecies(["A"])
        A.x
        A.y
        A.z

        result = get_all_non_listed_characteristics(A, {"x", "all$"})
        for combo in result:
            assert "all$" in combo

    def test_multi_dimension_via_inheritance(self) -> None:
        """Inherited species with separate dimensions produce cartesian product."""
        Color, Size = BaseSpecies(["Color", "Size"])
        Color.red
        Color.blue
        Size.big
        Size.small

        Thing = (Color * Size).named("Thing")

        result = get_all_non_listed_characteristics(Thing, {"red"})
        combos = {frozenset(s) for s in result}

        assert frozenset({"blue", "big"}) in combos
        assert frozenset({"blue", "small"}) in combos
        assert len(combos) == 2

    def test_multi_dimension_exclude_from_each(self) -> None:
        """Excluding one from each dimension narrows both."""
        Color, Size = BaseSpecies(["Color", "Size"])
        Color.red
        Color.blue
        Color.green
        Size.big
        Size.small

        Thing = (Color * Size).named("Thing")

        result = get_all_non_listed_characteristics(Thing, {"red", "big"})
        combos = {frozenset(s) for s in result}

        assert frozenset({"blue", "small"}) in combos
        assert frozenset({"green", "small"}) in combos
        assert len(combos) == 2


# ---------------------------------------------------------------
# new_reaction_with_new_characteristics
# ---------------------------------------------------------------


class TestNewReactionWithNewCharacteristics:
    """Tests for new_reaction_with_new_characteristics."""

    def test_replaces_reactant_characteristics(self) -> None:
        (A,) = BaseSpecies(["A"])
        A.x
        A.y
        (B,) = BaseSpecies(["B"])

        reaction = _build_reaction(A, {"x"}, B, set(), rate=2.0)
        new_r = new_reaction_with_new_characteristics(reaction, A, {"y"})

        reactants = _reactant_species_names(new_r)
        assert len(reactants) == 1
        assert reactants[0][0] is A
        assert reactants[0][1] == {"y"}

    def test_replaces_product_characteristics(self) -> None:
        (A,) = BaseSpecies(["A"])
        A.x
        A.y
        (B,) = BaseSpecies(["B"])

        reaction = _build_reaction(B, set(), A, {"x"}, rate=1.0)
        new_r = new_reaction_with_new_characteristics(reaction, A, {"y"})

        products = _product_species_names(new_r)
        assert len(products) == 1
        assert products[0][0] is A
        assert products[0][1] == {"y"}

    def test_preserves_rate(self) -> None:
        A, B = BaseSpecies(["A", "B"])
        reaction = _build_reaction(A, set(), B, set(), rate=3.14)
        new_r = new_reaction_with_new_characteristics(reaction, A, set())

        assert new_r.rate == pytest.approx(3.14)

    def test_unrelated_species_unchanged(self) -> None:
        """Species not matching spe_to_modify keep original characteristics."""
        A, B = BaseSpecies(["A", "B"])
        A.x
        A.y
        B.p
        B.q

        reaction = _build_reaction(A, {"x"}, B, {"p"}, rate=1.0)
        new_r = new_reaction_with_new_characteristics(reaction, A, {"y"})

        reactants = _reactant_species_names(new_r)
        products = _product_species_names(new_r)
        assert reactants[0][1] == {"y"}
        assert products[0][1] == {"p"}

    def test_both_sides_modified_when_same_species(self) -> None:
        """When the same species appears on both sides, both are updated."""
        (A,) = BaseSpecies(["A"])
        A.alive
        A.dead

        reaction = _build_reaction(A, {"alive"}, A, {"alive"}, rate=0.5)
        new_r = new_reaction_with_new_characteristics(reaction, A, {"dead"})

        reactants = _reactant_species_names(new_r)
        products = _product_species_names(new_r)
        assert reactants[0][1] == {"dead"}
        assert products[0][1] == {"dead"}


# ---------------------------------------------------------------
# create_all_not_reactions
# ---------------------------------------------------------------


class TestCreateAllNotReactions:
    """Tests for create_all_not_reactions."""

    def test_no_not_reactions_unchanged(self) -> None:
        """Reactions without NOT markers pass through unmodified."""
        A, B = BaseSpecies(["A", "B"])
        reaction = _build_reaction(A, set(), B, set(), rate=1.0)
        reactions = {reaction}

        result = create_all_not_reactions(reactions)
        assert len(result) == 1
        assert reaction in result

    def test_simple_not_expansion(self) -> None:
        """NOT on a single dimension expands to all alternatives."""
        (A,) = BaseSpecies(["A"])
        A.healthy
        A.sick
        A.dead
        (B,) = BaseSpecies(["B"])

        reaction = _build_reaction(A, {"sick", NOT_CHAR}, B, set(), rate=1.0)
        reactions = {reaction}

        result = create_all_not_reactions(reactions)

        assert reaction not in result
        assert len(result) == 2

        reactant_chars = {frozenset(r.reactants[0]["characteristics"]) for r in result}
        assert frozenset({"healthy"}) in reactant_chars
        assert frozenset({"dead"}) in reactant_chars

    def test_not_with_single_alternative(self) -> None:
        """NOT excluding all but one yields a single expanded reaction."""
        (A,) = BaseSpecies(["A"])
        A.on
        A.off
        (B,) = BaseSpecies(["B"])

        reaction = _build_reaction(A, {"on", NOT_CHAR}, B, set(), rate=1.0)
        reactions = {reaction}

        result = create_all_not_reactions(reactions)
        assert len(result) == 1
        rxn = next(iter(result))
        assert rxn.reactants[0]["characteristics"] == {"off"}

    def test_not_on_product_side(self) -> None:
        """NOT on the product side also expands correctly."""
        (A,) = BaseSpecies(["A"])
        A.x
        A.y
        A.z
        (B,) = BaseSpecies(["B"])

        reaction = _build_reaction(B, set(), A, {"x", NOT_CHAR}, rate=1.0)
        reactions = {reaction}

        result = create_all_not_reactions(reactions)
        assert len(result) == 2

        product_chars = {frozenset(r.products[0]["characteristics"]) for r in result}
        assert frozenset({"y"}) in product_chars
        assert frozenset({"z"}) in product_chars

    def test_not_preserves_other_reactions(self) -> None:
        """Non-NOT reactions in the set are preserved alongside expansions."""
        A, B = BaseSpecies(["A", "B"])
        A.p
        A.q
        A.r

        normal = _build_reaction(B, set(), A, set(), rate=0.5)
        not_rxn = _build_reaction(A, {"p", NOT_CHAR}, B, set(), rate=1.0)
        reactions = {normal, not_rxn}

        result = create_all_not_reactions(reactions)

        assert normal in result
        assert len(result) == 3
