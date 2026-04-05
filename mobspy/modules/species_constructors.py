"""Factory functions for creating meta-species: BaseSpecies, New, Zero."""

from __future__ import annotations

import linecache
import sys
from typing import overload

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
        hint = "Pass names explicitly: A, B = BaseSpecies(['A', 'B'])"
        if hasattr(sys, "ps1") or "ipykernel" in sys.modules:
            hint += (
                "\nNote: In interactive sessions (REPL/notebook), "
                "name inference from source is not always available."
            )
        raise ValidationError("Could not infer species names from source.\n" + hint)

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

    new_names = [clean_species_name(name) for name in names]

    return n, new_names


def _validate_number_or_names(number_or_names: int | list[str] | None) -> None:
    """Validate the number_or_names argument for species creation."""
    if number_or_names is None:
        return
    if isinstance(number_or_names, int):
        if number_or_names < 1:
            raise ValidationError(
                "Please use strictly positive integers for number of properties"
            )
    elif not isinstance(number_or_names, list):
        raise ValidationError("Only numbers or lists of strings accepted")


def _resolve_names(
    code_line: str,
    number_or_names: int | list[str] | None,
) -> tuple[int, list[str]]:
    """Resolve the number of species and their names from the caller's context."""
    if isinstance(number_or_names, list):
        return len(number_or_names), number_or_names

    number_of_properties, compiled_names = compile_species_number_line(code_line)
    if number_or_names is not None and number_of_properties != number_or_names:
        raise ValidationError(
            f"The number of properties ({number_or_names}) is not equal to "
            f"the number of variables ({number_of_properties}): "
            f"inferred names {compiled_names}"
        )
    return number_of_properties, compiled_names


def _Create_Species(  # noqa: N802
    species: Species | None,
    code_line: str,
    number_or_names: int | list[str] | None = None,
) -> Species | tuple[Species, ...]:
    """Instantiate one or more Species from the caller's line.

    Infers names from the assignment source.
    """
    _validate_number_or_names(number_or_names)
    number_of_properties, names = _resolve_names(code_line, number_or_names)

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
            temp = One * species  # type: ignore[operator]
            temp.name(name)  # type: ignore[union-attr,operator]
            to_return.append(temp)  # type: ignore[arg-type]

    if len(to_return) == 1:
        return to_return[0]
    return tuple(to_return)


@overload
def BaseSpecies() -> Species: ...


@overload
def BaseSpecies(number_or_names: list[str]) -> tuple[Species, ...]: ...


@overload
def BaseSpecies(number_or_names: int) -> tuple[Species, ...]: ...


def BaseSpecies(  # noqa: N802
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


@overload
def New(species: Species) -> Species: ...


@overload
def New(species: Species, number_or_names: list[str]) -> tuple[Species, ...]: ...


@overload
def New(species: Species, number_or_names: int) -> tuple[Species, ...]: ...


def New(  # noqa: N802
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


def ListSpecies(  # noqa: N802
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
