from __future__ import annotations

from itertools import product as ite_product
from typing import Any

from mobspy.modules.meta_class import Reacting_Species, Zero

"""
    Species and meta-species supporting scripts
"""


def create_all_not_reactions(reactions: set[Any]) -> set[Any]:
    new_reactions_set = reactions

    def include_new_combinations(r: Any, attribute_to_get: str) -> bool:
        ignore_flag = True
        for x in getattr(r, attribute_to_get):
            if "not$" in x["characteristics"]:
                new_reactions_set.remove(r)
                ignore_flag = False
                combinations = get_all_non_listed_characteristics(
                    x["object"],
                    x["characteristics"] - {"not$"},
                )

                for comb in combinations:
                    x["characteristics"] = comb
                    new_r = new_reaction_with_new_characteristics(r, x["object"], comb)
                    new_reactions_set.add(new_r)

        return ignore_flag

    has_not = True
    while has_not:
        has_not = False
        loop_reactions = set(new_reactions_set)  # Snapshot
        for r in loop_reactions:
            if not include_new_combinations(r, "reactants"):
                has_not = True
            if not include_new_combinations(r, "products"):
                has_not = True

    return new_reactions_set


def get_all_non_listed_characteristics(
    species: Any,
    characteristics: set[str],
) -> list[set[str]]:
    """
    Gets all characteristics related to a species except the
    ones listed below. It uses inheritance to find all
    characteristics linked to a species
    """
    operator_characteristics = {c for c in characteristics if "$" in c}

    multi_list_char_struct: list[set[str]] = []
    for spe in species.get_references():
        if characteristics:
            dimension_cha = {
                e for e in spe.get_characteristics() if e not in characteristics
            }
        else:
            dimension_cha = set()

        if dimension_cha:
            multi_list_char_struct.append(dimension_cha)

    combinations = [
        set(combo) | operator_characteristics
        for combo in ite_product(*multi_list_char_struct)
    ]

    return combinations


def new_reaction_with_new_characteristics(
    r: Any,
    spe_to_modify: Any,
    new_characteristics: set[str],
) -> Any:
    """
    Create a copy of this reaction, replacing characteristics of a specific dict

    :param dict_to_modify: The specific reactant/product dict to modify
    :param new_characteristics: New characteristics set for that dict
    :return: New Reactions object
    """
    # Build new reactants list
    React: Any = Zero
    for reactant in r.reactants:
        if reactant["object"] == spe_to_modify:
            rs = Reacting_Species(
                reactant["object"],
                new_characteristics,
                reactant["stoichiometry"],
                reactant["label"],
            )
        else:
            rs = Reacting_Species(
                reactant["object"],
                reactant["characteristics"],
                reactant["stoichiometry"],
                reactant["label"],
            )

        React = React + rs

    # Build new products list
    Product: Any = Zero
    for product in r.products:
        if product["object"] == spe_to_modify:
            ps = Reacting_Species(
                product["object"],
                new_characteristics,
                product["stoichiometry"],
                product["label"],
            )
        else:
            ps = Reacting_Species(
                product["object"],
                product["characteristics"],
                product["stoichiometry"],
                product["label"],
            )

        Product = Product + ps

    # Create new reaction
    return React >> Product[r.rate]
