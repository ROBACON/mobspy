"""Factory functions for creating meta-species: BaseSpecies, New, Zero."""

from __future__ import annotations

import linecache
import sys
from typing import TYPE_CHECKING

from mobspy.exceptions import ReactionError, ValidationError
from mobspy.modules.assignments_implementation import (
    Assign as asgi_Assign,
)
from mobspy.modules.list_species import List_Species
from mobspy.modules.reactions import _Last_rate_storage
from mobspy.modules.species import Species, clean_species_name

if TYPE_CHECKING:
    pass


def compile_species_number_line(code_line: str) -> tuple[int, list[str]]:
    """Compile code line for BaseSpecies and New.

    :param code_line: Line of code where BaseSpecies or New was called
    :return: (number of variables, list of name strings)
    """
    before_eq, after_eq = code_line.split("=")[0], code_line.split("=")[1]
    if after_eq.count("BaseSpecies") > 1:
        raise ReactionError(
            f"At {after_eq}: \n" + "BaseSpecies can only be called once at a time"
        )

    code_line = before_eq
    n = code_line.count(",") + 1
    code_line = code_line.replace(" ", "")
    names = code_line.split(",")

    new_names = []
    for name in names:
        new_names.append(clean_species_name(name))

    return n, new_names


def _Create_Species(
    species: Species | None,
    code_line: str,
    number_or_names: int | list[str] | None = None,
) -> Species | tuple[Species, ...]:
    if number_or_names is not None:
        if isinstance(number_or_names, int):
            if number_or_names < 1:
                raise ValidationError(
                    "Please use strictly positive integers for number of properties"
                )
        elif isinstance(number_or_names, list):
            pass
        else:
            raise ValidationError("Only numbers or lists of strings accepted")

    if number_or_names is None or isinstance(number_or_names, int):
        number_of_properties, compiled_names = compile_species_number_line(code_line)
        names = compiled_names
        if number_or_names is not None:  # noqa: SIM102
            if number_of_properties != number_or_names:
                raise ValidationError(
                    "The number of properties is not equal to the number of variables"
                )
    elif isinstance(number_or_names, list):
        number_of_properties = len(number_or_names)
        names = number_or_names

    to_return = []
    for i in range(number_of_properties):
        _Last_rate_storage.entity_counter += 1
        if names is None:
            name = "N$" + str(_Last_rate_storage.entity_counter)
        else:
            name = names[i]
        if species is None:
            to_return.append(Species(name))
        else:
            temp = One * species
            temp.name(name)
            to_return.append(temp)

    if len(to_return) == 1:
        return to_return[0]
    else:
        return tuple(to_return)


def BaseSpecies(
    number_or_names: int | list[str] | None = None,
) -> Species | tuple[Species, ...]:
    """Return base species with no inheritance.

    :param number_or_names: (int) number of base species
        or (list) list of names
    :return: Species objects
    """
    asgi_Assign.reset_context()

    frame = sys._getframe(1)
    code_line = linecache.getline(frame.f_code.co_filename, frame.f_lineno).rstrip("\n")
    return _Create_Species(None, code_line, number_or_names)


def New(
    species: Species,
    number_or_names: int | list[str] | None = None,
) -> Species | tuple[Species, ...]:
    """Return meta-species that inherit from the supplied species.

    :param species: Species to inherit from
    :param number_or_names: (int) number of meta-species
        or (list) list of names
    :return: Species objects
    """
    frame = sys._getframe(1)
    code_line = linecache.getline(frame.f_code.co_filename, frame.f_lineno).rstrip("\n")
    return _Create_Species(species, code_line, number_or_names)


def ListSpecies(
    number_of_elements: int,
    inherits_from: Species | None = None,
) -> List_Species:
    frame = sys._getframe(1)
    code = linecache.getline(frame.f_code.co_filename, frame.f_lineno).rstrip("\n")
    before_eq, after_eq = code.split("=")
    before_eq = before_eq.replace(" ", "")
    if "," in before_eq:
        raise ReactionError(
            f"At: {code} \n"
            "No comas are allowed in the right-side of "
            "the equality during the creation of "
            "a ListSpecies"
        )
    name = before_eq.split("=")[0]

    temp_list = []
    for i in range(number_of_elements):
        temp_name = name + "_" + str(i + 1)

        if inherits_from is None:
            temp_list.append(BaseSpecies([temp_name]))
        else:
            temp_list.append(New(inherits_from, [temp_name]))

    return List_Species(temp_list)


# Module-level singletons created at import time
Zero, One, EndFlagSpecies = BaseSpecies(3)
Zero._bypass_name("_S0")
One._bypass_name("_S1")
EndFlagSpecies._bypass_name("_End_Flag_MetaSpecies")
