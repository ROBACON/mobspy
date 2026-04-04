"""Dot-notation species operators and boolean overrides for rate expression queries."""

from __future__ import annotations

from contextvars import ContextVar
from typing import TYPE_CHECKING

from mobspy.constants import DOT_SEPARATOR, NULL_SPECIES
from mobspy.exceptions import CompilationError

if TYPE_CHECKING:
    from mobspy.modules.species import Species

__all__ = [
    "Bool_Override",
    "Specific_Species_Operator",
    "_ms_active_ctx",
]

# Global context variable for expression mode.
# When True, arithmetic on ExpressionDefiner subclasses builds expression
# trees instead of performing plain numeric/unit operations.
_ms_active_ctx: ContextVar[bool] = ContextVar("_ms_active_ctx", default=False)


class Bool_Override:  # noqa: N801
    """
    Just a base class for implementing the . operation in the rate
    function arguments through boolean overriding. It is responsible
    for returning true when the reactant has the specified
    characteristic when using the dot notation

    Args:
        _stocked_characteristics: Stocks the characteristics of the queries performed by
            the user.
        species_string: String value from an individual species in MobsPY format.
    """

    species_string: str
    _stocked_characteristics: set[str]

    def __bool__(self) -> bool:
        """
        The implementation of the .dot operation for rate function arguments
        Returns true when the argument possesses the characteristics
        Bool is called after the .dot operations

        Parameters:
            self

        Returns:
            True if the object contains all the characteristics queried
            False otherwise
        """
        if self.species_string == NULL_SPECIES:
            return False

        species_string_split = self.species_string.split(DOT_SEPARATOR)[1:]
        if all(char in species_string_split for char in self._stocked_characteristics):
            to_return_boolean = True
        else:
            to_return_boolean = False

        self._stocked_characteristics = set()
        return to_return_boolean


class Specific_Species_Operator(Bool_Override):  # noqa: N801
    """
    Creates objects from species strings from the meta-species to
    pass them to rate functions as arguments. Uses Bool_Override to
    return true or false to the .dot operation inside rate functions.

    Args:
        _stocked_characteristics: Stocks the characteristics of the queries performed by
            the user.
        species_string: String value from an individual species in MobsPY format.
        species_object: Meta-species object which originated the meta-species str.
    """

    def __init__(self, species_string: str, species_object: Species | None) -> None:
        """
        Constructs the object from the species strings from the
        meta-species to pass them to rate functions as arguments.

        Args:
            species_string: A string from MobsPy meta-species format.
            species_object: Meta-species set for which the species_string is contained
                in.
        """
        self.species_string = species_string
        self._stocked_characteristics: set[str] = set()
        self._species_object = species_object

    def __getattr__(self, characteristic: str) -> Specific_Species_Operator:
        """
        Stores the characteristics for the boolean query inside the
        rate function by adding them to the set.

        Args:
            characteristic: Characteristic being queried.
        """
        self._stocked_characteristics.add(characteristic)
        return self

    def __str__(self) -> str:
        """
        Returns the species_string from the MobsPy meta-species
        used in the object construction.
        """
        return self.species_string

    def is_a(self, reference: Species) -> bool | None:
        """
        Checks if the meta-species the species_string belongs to has
        inherited from the parameter reference (reminder: every
        meta-species inherits from itself).

        Args:
            reference: Meta-species object.

        Returns:
            True if the meta-species in Specific_Species_Operator has inherited from the
            reference, False otherwise.
        """
        if not self._stocked_characteristics:
            if self._species_object is None:
                raise CompilationError("Species object is not set for this operator")
            return reference in self._species_object.get_references()
        raise CompilationError(
            "Cannot chain is_a() with dot-notation characteristic queries. "
            "Use them in separate conditions: "
            "'r.is_a(X) and r.alive' instead of chaining."
        )

    def add(self, characteristic: str) -> None:
        """
        Adds characteristic to the set of characteristics when checking for all
        referred characteristics by the user

        Args:
            characteristic: Characteristic to add to the set.
        """
        self._stocked_characteristics.add(characteristic)

    def get_name(self) -> str:
        """
        Returns: The name of the species the reactant is in string format
        """
        if self._species_object is None:
            raise CompilationError("Species object is not set for this operator")
        return self._species_object.get_name()

    def get_characteristics(self) -> set[str]:
        """
        Returns: The characteristics of the species in this given state
        """
        return set(self.species_string.split(DOT_SEPARATOR)[1:])

    def get_state(self) -> str:
        """
        Returns: the string of a state with dots instead of _dot_
        """
        return self.species_string.replace(DOT_SEPARATOR, ".")
