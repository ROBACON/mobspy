"""Assign initial species counts from user calls, resolving inheritance and units."""

from __future__ import annotations

import sys
from types import FrameType
from typing import Any

from numpy import floating as np_float_
from numpy import integer as np_int_
from pint import Quantity

from mobspy.constants import ALL_CHAR
from mobspy.exceptions import ValidationError
from mobspy.modules.list_species import List_Species
from mobspy.modules.mobspy_parameters import (
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)
from mobspy.modules.reactions import Reacting_Species
from mobspy.modules.species import Species


def set_counts(count_dic: dict[Any, Any]) -> List_Species:
    """Adds counts to meta-species using a given dictionary.

    Keys from this dictionary can be either meta-species objects or
    strings. Items must be the assigned counts to species.

    Args:
        count_dic: Dictionary where keys can be either meta-species or strings.


    Returns:
        List_Species object. All meta-species that had a count assigned in this
        dictionary will be returned as a List_Species which can be passed as a model to
        the simulation object.
    """
    count_dic = _validate_count_values(count_dic)

    all_found_species: set[Species] = set()
    for key in count_dic:
        if isinstance(key, str):
            all_found_species = _find_species_in_stack()
            break

    model: set[Any] = set()
    for key, item in count_dic.items():
        if isinstance(key, str):
            _resolve_string_key(key, item, all_found_species, model)
        else:
            _resolve_species_key(key, item, model)

    return List_Species(model)


def _validate_count_values(count_dic: dict[Any, Any]) -> dict[Any, Any]:
    """Validate and normalize count values."""
    new_count_dict: dict[Any, Any] = {}
    for key, item in count_dic.items():
        if isinstance(item, (int, float, Quantity, mp_Mobspy_Parameter)):
            new_count_dict[key] = item
        elif isinstance(item, (np_int_, np_float_)):
            new_count_dict[key] = float(item)
        else:
            raise ValidationError(
                "Reactant_species count assignment does not"
                f" support the type {type(item)}"
            )
    return new_count_dict


def _find_species_in_stack() -> set[Species]:
    """Walk the call stack to collect all Species instances."""
    found_species: set[Species] = set()
    frame: FrameType | None = sys._getframe(1)
    while frame is not None:
        for ns in (frame.f_locals, frame.f_globals):
            for obj in ns.values():
                try:
                    if isinstance(obj, Species) and type(obj) != type:  # noqa: E721
                        found_species.add(obj)
                except AttributeError:
                    pass
        frame = frame.f_back
    return found_species


def _resolve_string_key(
    key: str,
    item: Any,
    all_found_species: set[Species],
    model: set[Any],
) -> None:
    """Resolve a string key to a species and assign count."""
    already_found = False
    str_name = key.split(".", maxsplit=1)[0]
    str_characteristics = set(key.split(".")[1:])
    for spe in all_found_species:
        if spe.get_name() == str_name and not already_found:
            already_found = True
            temp_set = set(str_characteristics)
            temp_set.discard(ALL_CHAR)
            if temp_set.issubset(spe.get_all_characteristics()):
                spe.add_quantities(str_characteristics, item)
            else:
                raise ValidationError(
                    "Characteristics not found in species with equal name"
                )
            model.add(spe)
        elif spe.get_name() == str_name and already_found:
            raise ValidationError(
                "There are two different meta-species with"
                " the same name. Set_counts cannot resolve"
            )
    if not already_found:
        raise ValidationError(f"Meta-species with the following name {key} not found")


def _resolve_species_key(key: Any, item: Any, model: set[Any]) -> None:
    """Resolve a Species or Reacting_Species key and assign count."""
    try:
        if isinstance(key, (Species, Reacting_Species)):
            if not isinstance(key, Species):
                if len(key.list_of_reactants) != 1:
                    raise ValidationError(
                        "Assignment used incorrectly. Only one species at a time"
                    )
                model.add(key.list_of_reactants[0]["object"])
            if isinstance(key, Species):
                model.add(key)
            key(item)
    except AttributeError as e:
        raise ValidationError("Keys must be either meta-species or strings") from e
