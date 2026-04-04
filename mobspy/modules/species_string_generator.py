"""Generate concrete species name strings from characteristic combinations."""

from __future__ import annotations

from itertools import product as itertools_product
from typing import Any

from mobspy.constants import STD_CHAR
from mobspy.types import ConcreteSpeciesId


def characteristics_dictionary(
    characteristics: set[str],
    characteristics_to_object: dict[str, Any],
) -> dict[Any, str]:
    """
    This function constructs a dictionary that leads from the
    species string characteristics to their respective object.
    This structure allows to easily find where the
    characteristics in a species string are locates

    Args:
        characteristics: Of characteristics.
        characteristics_to_object: With characteristics as keys and their respective
            base meta-species object as value.
    """
    object_to_characteristic: dict[Any, str] = {}
    for characteristic in characteristics:
        if "$" not in characteristic:
            object_to_characteristic[characteristics_to_object[characteristic]] = (
                characteristic
            )
    return object_to_characteristic


def construct_species_char_list(
    spe_or_reactiong_object: Any,
    characteristics: set[str] | str,
    characteristics_to_object: dict[str, Any],
    symbol: str | None = None,
) -> list[Any] | str:
    """
    This function constructs a list in the format
    ['species_name', 'char1', 'char2', ...]. It generates
    this list for a given meta-species and the specified
    characteristics. Values of characteristics not specified
    in a particular position are replaced by their default
    value. If a symbol is given it generates a string from
    the list using the symbol to join it.

    Args:
        spe_object: Meta-species object to be used.
        characteristics: Of characteristics given.
        characteristics_to_object: With characteristics as keys and their respective
            base meta-species object as value.
        symbol: Usually . or _dot_, connects the elements from the list using the
            specified symbol.
    """
    if characteristics == STD_CHAR or isinstance(characteristics, str):
        char_set: set[str] = set()
    else:
        char_set = characteristics
    spe_object = spe_or_reactiong_object.get_spe_object()

    ordered_references_list = spe_object.get_ordered_references()

    objects_to_characteristic = characteristics_dictionary(
        char_set, characteristics_to_object
    )

    species_char_list: list[Any] = [spe_object]

    for obj in ordered_references_list:
        if obj in objects_to_characteristic:
            species_char_list.append(objects_to_characteristic[obj])
        else:
            species_char_list.append(obj.first_characteristic)

    if symbol is not None:
        result: str = (
            symbol.join([spe_object.get_name(), *species_char_list[1:]])
            if len(species_char_list) > 1
            else spe_object.get_name()
        )
        return result

    return species_char_list


def construct_all_combinations(
    spe_or_reactiong_object: Any,
    characteristics: set[str] | str,
    characteristics_to_object: dict[str, Any],
    symbol: str | None = None,
) -> list[Any]:
    """
    This function constructs all possible list in the format
    ['species_name', 'char1', 'char2', ...] using all
    combinations of characteristics from the vector
    coordinates not used in the characteristics specified
    in the function argument. If a symbol is given it
    generates a string from the list using the symbol
    to join it.

    Args:
        spe_object: Meta-species object to be used.
        characteristics: Of characteristics given.
        characteristics_to_object: With characteristics as keys and their respective
            base meta-species object as value.
        symbol: Usually . or _dot_, connects the elements from the list using the
            specified symbol.
    """

    if characteristics == STD_CHAR or isinstance(characteristics, str):
        char_set2: set[str] = set()
    else:
        char_set2 = characteristics
    spe_object = spe_or_reactiong_object.get_spe_object()

    spe_object.order_references()
    ordered_references_list = spe_object.get_ordered_references()

    objects_to_characteristic = characteristics_dictionary(
        char_set2, characteristics_to_object
    )

    list_of_all_possibilities: list[list[Any]] = [[spe_object]]
    for obj in ordered_references_list:
        if obj in objects_to_characteristic:
            list_of_all_possibilities.append([objects_to_characteristic[obj]])
        else:
            list_of_all_possibilities.append(list(obj.get_characteristics()))

    to_return: list[Any] = []
    for i in itertools_product(*list_of_all_possibilities):
        if symbol is not None:
            to_return.append(
                symbol.join([spe_object.get_name(), *list(i)[1:]])
                if len(i) > 1
                else spe_object.get_name()
            )
        else:
            to_return.append(list(i))

    return to_return


# ------------------------------------------------------------------
# Structured variants returning ConcreteSpeciesId
# ------------------------------------------------------------------


def construct_species_id(
    spe_or_reacting_object: Any,
    characteristics: set[str] | str,
    characteristics_to_object: dict[str, Any],
) -> ConcreteSpeciesId:
    """Build a single ConcreteSpeciesId for the default state.

    Structured alternative to ``construct_species_char_list``
    with ``symbol=DOT_SEPARATOR``.
    """
    char_list = construct_species_char_list(
        spe_or_reacting_object, characteristics, characteristics_to_object
    )
    if isinstance(char_list, str):
        return ConcreteSpeciesId.from_sbml_id(char_list)
    spe_obj = char_list[0]
    return ConcreteSpeciesId(
        base=spe_obj.get_name(),
        characteristics=tuple(str(c) for c in char_list[1:]),
    )


def construct_all_species_ids(
    spe_or_reacting_object: Any,
    characteristics: set[str] | str,
    characteristics_to_object: dict[str, Any],
) -> list[ConcreteSpeciesId]:
    """Build all ConcreteSpeciesId combinations.

    Structured alternative to ``construct_all_combinations``
    with ``symbol=DOT_SEPARATOR``.
    """
    raw = construct_all_combinations(
        spe_or_reacting_object, characteristics, characteristics_to_object
    )
    result: list[ConcreteSpeciesId] = []
    for entry in raw:
        if isinstance(entry, list):
            spe_obj = entry[0]
            result.append(
                ConcreteSpeciesId(
                    base=spe_obj.get_name(),
                    characteristics=tuple(str(c) for c in entry[1:]),
                )
            )
        else:
            result.append(ConcreteSpeciesId.from_sbml_id(str(entry)))
    return result
