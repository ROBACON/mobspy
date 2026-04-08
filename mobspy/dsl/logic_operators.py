"""Comparison and logical operators for building event triggers on species counts."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeAlias

from mobspy.compiler.species_strings import (
    construct_all_combinations as ssg_construct_all_combinations,
)
from mobspy.constants import DOT_SEPARATOR
from mobspy.exceptions import EventError

if TYPE_CHECKING:
    from pint import Quantity

    from mobspy.simulation import Simulation

_SpeciesDict: TypeAlias = dict[str, Any]
_OpElem: TypeAlias = "_SpeciesDict | int | float | str | Quantity"
_CompNum: TypeAlias = "int | float | SpeciesComparator | Quantity"
_ScalarNum: TypeAlias = "int | float | Quantity"


class SpeciesComparator:
    """This class implements the comparisons necessary
    for events and conditional durations for Species.
    Ex: (A > 5).

    Args:
        _simulation_context: Current simulation under context.
    """

    def __init__(self) -> None:
        self._simulation_context: Simulation | None = None

    def add_operation_and_number(
        self,
        symbol: str,
        number: _CompNum,
    ) -> MetaSpeciesLogicResolver:
        """Creates a MetaSpeciesLogicResolver from the
        comparison of a meta-species with a value or
        another meta-species.

        Args:
            symbol: Comparison symbol '>=', '<=', '>' or '<'.
            number: If compared to a value, (Species or ReactingSpecies) if compared to
                the objects.
        """
        if self.__dict__.get("_from_mul", False):
            raise EventError(
                "Species created by the * operator (inheritance) "
                "cannot be used directly in event comparisons.\n"
                "The * operator creates inheritance, not multiplication."
            )
        reformatted = self.reformat_number_and_species(number)
        if isinstance(reformatted, list):
            operation: list[_OpElem] = [
                self.logical_add_species(self),
                symbol,
                *reformatted,
            ]
        else:
            operation = [
                self.logical_add_species(self),
                symbol,
                reformatted,
            ]
        return MetaSpeciesLogicResolver(operation, self._simulation_context)

    def reformat_number_and_species(
        self,
        number: _CompNum,
    ) -> int | float | _SpeciesDict | list[_OpElem] | Quantity:
        """Discovers if the number is a Species,
        ReactingSpecies or integer and prepares the
        output accordingly.

        Args:
            number: If compared to a value, (Species or ReactingSpecies) if compared to
                the objects.
        """
        if isinstance(number, SpeciesComparator):
            if number.is_species():  # type: ignore[attr-defined]  # mixin attr
                return self.logical_add_species(number)
            return self.logical_add_reacting_species(number)
        return number

    @classmethod
    def logical_add_species(cls, species: SpeciesComparator) -> _SpeciesDict:
        """Adds a species object to an event trigger
        operation.

        Args:
            species: Species object.
        """
        return {
            "object": species,
            "characteristics": set(),
        }

    @classmethod
    def logical_add_reacting_species(
        cls, react_spe: SpeciesComparator
    ) -> list[_OpElem]:
        """Adds a reacting species object to an event
        trigger operation by summing over all the
        indicated characteristics. It also accepts sums
        of species multiplied by integers.

        Args:
            react_spe: Object.
        """
        operation: list[_OpElem] = []

        react_spe.check_context()  # type: ignore[attr-defined]  # mixin attr
        for i, react_dict in enumerate(react_spe.list_of_reactants):  # type: ignore[attr-defined]  # mixin attr
            if i > 0:
                dl: list[_OpElem] = ["+"]
            else:
                dl = []
            dl = [*dl, react_dict["stoichiometry"], "*"]
            dl = [
                *dl,
                {
                    "object": react_dict["object"],
                    "characteristics": react_dict["characteristics"],
                },
            ]
            operation = operation + dl
        return operation

    def __lt__(self, number: _CompNum) -> MetaSpeciesLogicResolver:
        return self.add_operation_and_number("<", number)

    def __le__(self, number: _CompNum) -> MetaSpeciesLogicResolver:
        return self.add_operation_and_number("<=", number)

    def __gt__(self, number: _CompNum) -> MetaSpeciesLogicResolver:
        return self.add_operation_and_number(">", number)

    def __ge__(self, number: _CompNum) -> MetaSpeciesLogicResolver:
        return self.add_operation_and_number(">=", number)

    def __eq__(self, other: object) -> bool:
        if self._simulation_context is not None:
            raise EventError(
                "Equality assignment not allowed for "
                "event condition in MobsPy.\n"
                "Please if necessary "
                "use ( >= ) & ( =< )"
            )
        return id(self) == id(other)

    def __ne__(self, other: object) -> bool:
        return id(self) != id(other)

    def __hash__(self) -> int:
        return hash(id(self))


class ReactingSpeciesComparator(SpeciesComparator):
    """This class implements the comparisons necessary
    for events and conditional durations for Reacting
    Species. Ex: (A.a1 > 5).

    Args:
        _simulation_context: Current simulation under context.
    """

    def __init__(self) -> None:
        super().__init__()
        self._simulation_context: Simulation | None = None

    def check_context(self) -> None:
        """Checks if a reacting species is under
        context and adds it to the attribute
        self._simulation_context.
        """
        for react_dict in self.list_of_reactants:  # type: ignore[attr-defined]
            if react_dict["object"]._simulation_context is not None:
                self._simulation_context = react_dict["object"]._simulation_context
                break
        else:
            self._simulation_context = None

    def add_operation_and_number(
        self,
        symbol: str,
        number: _CompNum,
    ) -> MetaSpeciesLogicResolver:
        """Adds the symbol and number to create a
        MetaSpeciesLogicResolver object.

        Args:
            symbol: Comparative symbols '<=', '<', '>=', or '>'.
            number: If compared to a value, (Species or ReactingSpecies) if compared to
                the objects.


        Returns:
            MetaSpeciesLogicResolver object containing the comparison.
        """
        reformatted = self.reformat_number_and_species(number)
        if isinstance(reformatted, list):
            operation: list[_OpElem] = [
                *self.logical_add_reacting_species(self),
                symbol,
                *reformatted,
            ]
        else:
            operation = [
                *self.logical_add_reacting_species(self),
                symbol,
                reformatted,
            ]
        return MetaSpeciesLogicResolver(operation, self._simulation_context)


class MetaSpeciesLogicResolver:
    """This object stores the logical and comparison
    operations that become the event triggers for
    conditions or the conditional duration for
    simulations.

    Args:
        operation: The logical operation that will be transformed in a string for
            copasi. Format ex: [{object: ..., characteristics:...}, '<=', 10].
        simulation_context: Simulation under context.
    """

    def __init__(
        self,
        operation: list[_OpElem],
        simulation_context: Simulation | None = None,
    ) -> None:
        self.operation = operation
        self.simulation_context = simulation_context

    def __bool__(self) -> bool:
        raise EventError(
            "MetaSpeciesLogicResolver cannot be used in boolean context.\n"
            "This typically happens with Python chained comparisons "
            "(e.g., 10 >= A >= 10).\n"
            "Use logical operators instead: (A >= 10) & (A <= 20)"
        )

    def __and__(self, other: MetaSpeciesLogicResolver) -> MetaSpeciesLogicResolver:
        return self._join(other, "&&")

    def __or__(self, other: MetaSpeciesLogicResolver) -> MetaSpeciesLogicResolver:
        return self._join(other, "||")

    def _join(
        self,
        other: MetaSpeciesLogicResolver,
        symbol: str,
    ) -> MetaSpeciesLogicResolver:
        """Combine two logic expressions with the given boolean operator symbol."""
        if not isinstance(other, MetaSpeciesLogicResolver):
            raise EventError(
                "Logic operations require MetaSpeciesLogicResolver operands"
            )

        new_operation: list[_OpElem] = [
            "(",
            "(",
            *self.operation,
            ")",
            symbol,
            "(",
            *other.operation,
            ")",
            ")",
        ]
        return MetaSpeciesLogicResolver(new_operation, self.simulation_context)

    def __lt__(self, number: _ScalarNum) -> MetaSpeciesLogicResolver:
        return self._add_comparison("<", number)

    def __le__(self, number: _ScalarNum) -> MetaSpeciesLogicResolver:
        return self._add_comparison("<=", number)

    def __gt__(self, number: _ScalarNum) -> MetaSpeciesLogicResolver:
        return self._add_comparison(">", number)

    def __ge__(self, number: _ScalarNum) -> MetaSpeciesLogicResolver:
        return self._add_comparison(">=", number)

    def _add_comparison(
        self, symbol: str, number: _ScalarNum
    ) -> MetaSpeciesLogicResolver:
        """Create a new resolver with a comparison operator prepended.

        Raises:
            EventError: If a comparison operator already exists.
        """
        comparison_symbols = {"<", "<=", ">", ">="}
        has_cmp = any(
            op in comparison_symbols for op in self.operation if isinstance(op, str)
        )
        if has_cmp:
            raise EventError(
                "Chained comparisons are not supported in MobsPy events.\n"
                "Use logical operators instead: (A >= 10) & (A <= 20)"
            )
        new_operation = [number, symbol, *list(self.operation)]
        return MetaSpeciesLogicResolver(new_operation, self.simulation_context)

    @classmethod
    def find_all_species_strings(
        cls,
        species: SpeciesComparator,
        characteristics: set[str],
        species_for_sbml: dict[str, Any],
    ) -> list[str]:
        """Find all SBML species names matching the given characteristics."""
        reference_set = set(characteristics)
        reference_set.add(str(species))
        return [
            x
            for x in species_for_sbml
            if reference_set.issubset(set(x.split(DOT_SEPARATOR)))
        ]

    def generate_string(
        self,
        characteristics_to_object: dict[str, Any],
        to_sort: bool = False,
    ) -> str:
        """When a meta-species is used in a logic
        expression, this function transforms a
        meta-species in the sum of all individual
        states.

        Args:
            characteristics_to_object: Orthogonal characteristic space.
            to_sort: Sort strings or not - so the sum will always appear in the same
                order.
        """
        copasi_str = ""
        for _i, e in enumerate(self.operation):
            if isinstance(e, (int, float, str)):
                copasi_str = copasi_str + str(e) + " "
            else:
                copasi_str = copasi_str + "("
                ite = ssg_construct_all_combinations(
                    e["object"],
                    e["characteristics"],  # pyright: ignore[reportArgumentType]
                    characteristics_to_object,
                    DOT_SEPARATOR,
                )
                if to_sort:
                    ite = sorted(ite)

                for j, species_str in enumerate(ite):
                    if j == 0:
                        copasi_str = copasi_str + f"{species_str}"
                    else:
                        copasi_str = copasi_str + f" + {species_str}"
                copasi_str = copasi_str + ")" + " "

        return copasi_str
