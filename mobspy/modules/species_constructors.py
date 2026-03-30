"""Factory functions for creating meta-species: BaseSpecies, New, Zero."""

from __future__ import annotations

import linecache
import sys

from mobspy.constants import (
    END_FLAG_SPECIES_NAME,
    ONE_SPECIES_NAME,
    RATE_NAME_PREFIX,
    ZERO_SPECIES_NAME,
)
from mobspy.exceptions import ReactionError, ValidationError
from mobspy.modules.assignments_implementation import (
    Assign as asgi_Assign,
)
from mobspy.modules.list_species import List_Species
from mobspy.modules.reactions import _Last_rate_storage
from mobspy.modules.species import Species, clean_species_name


def _read_caller_source_line(stack_depth: int = 2) -> str:
    """Read the source line of the caller for name inference.

    Returns:
        The source line, or an empty string if unavailable.
    """
    try:
        frame = sys._getframe(stack_depth)
        line = linecache.getline(frame.f_code.co_filename, frame.f_lineno).rstrip("\n")
    except (ValueError, OSError):
        return ""
    return line


def compile_species_number_line(code_line: str) -> tuple[int, list[str]]:
    """Parse variable names from an assignment source line.

    Args:
        code_line: Line of code where BaseSpecies or New was called.

    Returns:
        Tuple of (count, names) extracted from the assignment.

    Raises:
        ValidationError: If names cannot be inferred from the source.
    """
    if "=" not in code_line:
        raise ValidationError(
            "Could not infer species names from source.\n"
            "Pass names explicitly: "
            "A, B = BaseSpecies(['A', 'B'])"
        )

    before_eq, after_eq = (
        code_line.split("=", maxsplit=1)[0],
        code_line.split("=")[1],
    )
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
    """Instantiate one or more Species from the caller's line.

    Infers names from the assignment source.
    """
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
        _Last_rate_storage.increment_entity_counter()
        if names is None:
            name = RATE_NAME_PREFIX + str(_Last_rate_storage.get_entity_counter())
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

    Args:
        number_or_names: Number of base species or (list) list of names.


    Returns:
        Species objects.

    Examples:
        >>> from mobspy import *
        >>> A, B = BaseSpecies(['A', 'B'])
        >>> A.get_name()
        'A'
        >>> A.is_species()
        True
    """
    asgi_Assign.reset_context()
    code_line = _read_caller_source_line(stack_depth=2)
    return _Create_Species(None, code_line, number_or_names)


def New(
    species: Species,
    number_or_names: int | list[str] | None = None,
) -> Species | tuple[Species, ...]:
    """Return meta-species that inherit from the supplied species.

    Args:
        species: Species to inherit from.
        number_or_names: Number of meta-species or (list) list of names.


    Returns:
        Species objects.

    Examples:
        >>> from mobspy import *
        >>> A = BaseSpecies(['A'])
        >>> B = New(A, ['B'])
        >>> B.get_name()
        'B'
        >>> B.is_species()
        True
    """
    code_line = _read_caller_source_line(stack_depth=2)
    return _Create_Species(species, code_line, number_or_names)


def ListSpecies(
    number_of_elements: int,
    inherits_from: Species | None = None,
) -> List_Species:
    """Create a numbered list of species, optionally inheriting from a parent."""
    code = _read_caller_source_line(stack_depth=2)
    if "=" not in code:
        raise ValidationError(
            "Could not infer ListSpecies name from source.\n"
            "Assign to a variable: my_list = ListSpecies(3)"
        )
    before_eq, _after_eq = code.split("=")
    before_eq = before_eq.replace(" ", "")
    if "," in before_eq:
        raise ReactionError(
            f"At: {code} \n"
            "No commas are allowed in the left-side of "
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

    return List_Species(temp_list)  # type: ignore[arg-type]


# Module-level singletons created at import time
Zero, One, EndFlagSpecies = BaseSpecies(3)
Zero._bypass_name(ZERO_SPECIES_NAME)
One._bypass_name(ONE_SPECIES_NAME)
EndFlagSpecies._bypass_name(END_FLAG_SPECIES_NAME)
