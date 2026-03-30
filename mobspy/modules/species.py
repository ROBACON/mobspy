"""Species class and related helpers.

The core meta-species object that users interact with to define
models in MobsPy's DSL.
"""

from __future__ import annotations

import linecache
import re
import sys
from contextvars import ContextVar
from typing import TYPE_CHECKING, Any, Self

from numpy import floating as np_float_
from numpy import integer as np_int_
from pint import Quantity

from mobspy.constants import DOT_SEPARATOR, NOT_CHAR, STD_CHAR
from mobspy.exceptions import ReactionError, ValidationError
from mobspy.mobspy_logging import get_logger
from mobspy.modules.assignments_implementation import (
    Asg as asgi_Asg,
)
from mobspy.modules.assignments_implementation import (
    Assign as asgi_Assign,
)
from mobspy.modules.logic_operators import (
    SpeciesComparator as lop_SpeciesComparator,
)
from mobspy.modules.mobspy_expressions import (
    Specific_Species_Operator as me_Specific_Species_Operator,
)
from mobspy.modules.mobspy_parameters import (
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)
from mobspy.modules.reactions import (
    Assignment_Opp_Imp,
    Reacting_Species,
    Reactions,
    _Last_rate_storage,
    _methods_Reacting_Species,
)
from mobspy.modules.species_string_generator import (
    construct_all_combinations as ssg_construct_all_combinations,
)
from mobspy.modules.species_utils import (
    check_orthogonality_between_references as mcu_check_orthogonality_between_references,  # noqa: E501
)
from mobspy.modules.species_utils import (
    combine_references as mcu_combine_references,
)
from mobspy.modules.species_utils import (
    unite_characteristics as mcu_unite_characteristics,
)

if TYPE_CHECKING:
    from collections.abc import Generator

    from mobspy.modules.list_species import List_Species

_logger = get_logger(__name__)

_simulation_context_cv: ContextVar[Any] = ContextVar(
    "_simulation_context_cv", default=None
)
_meta_any_ctx_cv: ContextVar[set[str]] = ContextVar("_meta_any_ctx_cv")


def _get_multiline_code_context(filename: str, lineno: int) -> str:
    """Extract multi-line code context for reaction parsing.

    Handles cases where formatters split >> operators across lines.
    """
    line = linecache.getline(filename, lineno)
    if not line:
        return ""

    return line.rstrip("\n")


# Easter Egg: I finished the first version on a sunday at the
# BnF in Paris. If anyone is reading this, I highly recommend
# you study there, it is quite a nice place


class Species(lop_SpeciesComparator, Assignment_Opp_Imp):
    """Fundamental class - The meta-species object.

    Contains the characteristics, the species name,
    the reactions it is involved in and finally the
    other species it references. Objects store all the
    basic information necessary to create an SBML file
    and construct a model.

    Args:
        _name: Name of the species.
        _characteristics: Set of characteristics DIRECTLY added to a species.
        _references: Set of meta-species a meta-species has inherited from.
        first_characteristic: First characteristic added to the species.
        _reactions: Every species stores all reactions it is involved in.
        _species_counts: Counts listed for the species.

    Examples:
        >>> from mobspy.modules.species import Species
        >>> A = Species("A")
        >>> A.get_name()
        'A'
        >>> A.is_species()
        True
    """

    def __init__(self, name: str) -> None:
        """Object constructor - We recommend using BaseSpecies instead.

        Args:
            name: Name of the species.
        """
        super().__init__()
        self.name(name)
        self._characteristics: set[str] = set()
        self._references: set[Species] = {self}
        self._ordered_references: list[Species] = []
        self._reference_index_dictionary: dict[Species, int] = {}
        self._unit: str = ""
        self._assignments: dict[str, Any] = {}
        self._linked_species: set[Species] = set()

        self.first_characteristic: str | None = None

        self._reactions: set[Reactions] = set()

        self._species_counts: list[dict[str, Any]] = []

    @classmethod
    def check_if_valid_characteristic(cls, affected_object: Any, char: str) -> bool:
        """Check if the characteristic name is valid.

        Args:
            affected_object: The object being checked.
            char: Name of the characteristic to be added.

        Raises:
            ValidationError: If the characteristic name is not allowed.


        Returns:
            True if characteristic is allowed.
        """
        black_list = {"list_of_reactants", "first_characteristic"}
        system_attrs = {"_pytestfixturefunction", "__sphinx_mock__"}
        if (
            char[0] == "_"
            and char not in system_attrs
            and not (char.startswith("__") and char.endswith("__"))
        ):
            raise ValidationError(
                f"Characteristic name {char} in object "
                f"{affected_object} is not allowed."
                " Please pick another name"
            )

        if (
            char in _methods_Reacting_Species
            or char in _methods_Species
            or char in black_list
        ):
            raise ValidationError(
                f"Characteristic name {char} in object "
                f"{affected_object} is not allowed. "
                "Please pick another name"
            )
            return False
        return True

    @classmethod
    def str_under_context(cls, species_object: Species, characteristics: Any) -> str:
        """Return the str representation of a species under context.

        Args:
            species_object: Meta-species object.
            characteristics: Characteristics to filter.


        Returns:
            String in format (A_dot_a1 + A_dot_a2 + ....).
        """
        ref_char_to_spe_obj = (
            Species.get_simulation_context().orthogonal_vector_structure
        )
        all_strings = sorted(
            ssg_construct_all_combinations(
                species_object, characteristics, ref_char_to_spe_obj, DOT_SEPARATOR
            )
        )
        to_str = all_strings[0]
        for i, e in enumerate(all_strings):
            if i == 0:
                continue
            to_str = to_str + " + " + e
        to_str = "(" + to_str + ")"
        return to_str  # type: ignore[no-any-return]

    def __str__(self) -> str:
        """String representation, returns the species name."""
        if Species.get_simulation_context() is None:
            return self._name
        else:
            return Species.str_under_context(self, STD_CHAR)

    def c(self, item: Any) -> Reacting_Species:
        """c query implementation, queries by value.

        Args:
            item: Value to query over.
        """
        item = str(item)
        Species.check_if_valid_characteristic(self, item)
        return self.__getattr__(item)  # type: ignore[no-any-return]

    def label(self, label: int | float | str) -> Reacting_Species:
        """Label function implementation.

        Args:
            label: Value for the label for matching.


        Returns:
            Reacting_Species object created with the label.
        """
        return Reacting_Species(self, set(), label=label)

    def show_reactions(self) -> None:
        """Print the reactions inside the object."""
        _logger.debug(str(self) + DOT_SEPARATOR)
        for reference in self._references:
            for reaction in reference.get_reactions():
                _logger.debug(str(reaction))

    def show_characteristics(self) -> None:
        """Print the characteristics directly added to this object."""
        _logger.debug(str(self) + " has the following characteristics referenced:")
        for i, reference in enumerate(self.get_references()):  # noqa: B007
            if reference.get_characteristics():
                _logger.debug(str(reference) + ": ")
                reference.show_characteristics()

    def show_references(self) -> None:
        """Print the objects this object has inherited from."""
        _logger.debug(str(self) + DOT_SEPARATOR)
        _logger.debug("{")
        for i, reference in enumerate(self.get_references()):  # noqa: B007
            if reference.get_characteristics():
                _logger.debug(" " + str(reference) + " ")
        _logger.debug("}")

    def show_quantities(self) -> None:
        """Show the species counts stored in this object."""
        _logger.info(str(self._species_counts))

    def __or__(self, other: Species | List_Species) -> List_Species:
        """Create an instance of List_Species using the ``|`` operator.

        Args:
            other: To combine.
        """
        from mobspy.modules.list_species import List_Species

        if isinstance(other, List_Species):
            other.append(self)
            return other
        elif isinstance(other, Species):
            return List_Species([self, other])
        else:
            raise ValidationError("Only Species and List_Species can be concatenated")

    def __iter__(self) -> Generator[Self, None, None]:
        """Iter defined to be consistent with List_Species behavior."""
        yield self

    def remove_repeated_elements(self) -> Self:
        """Defined to be consistent with the model behavior from List_Species."""
        return self

    def __getitem__(self, item: Any) -> Self:
        """Override of __getitem__ for dealing with reaction rates.

        Args:
            item: Reaction rate.
        """
        return _Last_rate_storage.override_get_item(self, item)  # type: ignore[no-any-return]

    def __rmul__(self, stoichiometry: Any) -> Reacting_Species | Any:
        """Multiplication by the stoichiometry.

        Args:
            stoichiometry: Stoichiometry.


        Returns:
            Reacting_Species with stoichiometry.
        """
        if not asgi_Assign.check_context():
            if isinstance(stoichiometry, (int, float)):
                r = Reacting_Species(self, set(), stoichiometry)
            else:
                raise ReactionError(
                    "Stoichiometry can only be an int or "
                    f"float - Received {stoichiometry}"
                )
            return r
        else:
            return asgi_Assign.mul(stoichiometry, self)

    def __add__(
        self,
        other: Species | Reacting_Species | Any,
    ) -> Reacting_Species | Any:
        """Addition for reaction construction.

        Args:
            other: Other object added to construct a reaction.


        Returns:
            Reacting Species from the sum.
        """
        if not asgi_Assign.check_context():
            r1 = Reacting_Species(self, set())
            if isinstance(other, Reacting_Species):
                r2 = other
            else:
                r2 = Reacting_Species(other, set())
            return r1 + r2
        else:
            return asgi_Assign.add(self, other)

    def __radd__(self, other: Any) -> Reacting_Species | Any:
        """Making addition symmetric, see __add__."""
        if not asgi_Assign.check_context():
            return Species.__add__(self, other)
        else:
            return asgi_Assign.add(other, self)

    def __invert__(self) -> Reacting_Species:
        return self.c(NOT_CHAR)

    def __neg__(self) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.mul(-1, self)
        else:
            raise ValidationError(
                "The negative operator was applied to a Species in the wrong context"
            )

    @classmethod
    def _compile_defined_reaction(cls, code_line: str, line_number: int) -> bool | None:
        """Validate that a reaction line ends with a rate.

        Raises on missing bracket-enclosed rate.
        """
        pattern = r"\][\s\n\)\],]*(#.*)?$"
        set_pattern = r"Set\s*\[.*>>.*\]"

        if re.search(set_pattern, code_line):
            return True

        if ">>" not in code_line:
            return True

        if code_line.rstrip().endswith(">>"):
            return True

        if code_line.rstrip().endswith("["):
            return True

        if not bool(re.search(pattern, code_line)):
            raise ReactionError(
                f"At: {code_line} \n"
                + f"Line number: {line_number} \n"
                + "There must be a rate in the end of the reaction. "
                "Avoid comments in the same line as the reaction."
            )
        return None

    def __rshift__(self, other: Species | Reacting_Species) -> Reactions:
        """Reaction definition (``>>`` operator).

        Args:
            other: Reaction products.


        Returns:
            The reaction.
        """
        myself = Reacting_Species(self, set())
        frame = sys._getframe(1)
        code_line = _get_multiline_code_context(
            frame.f_code.co_filename, frame.f_lineno
        )
        line_number = frame.f_lineno
        Species._compile_defined_reaction(code_line, line_number)

        if isinstance(other, Species):
            p = Reacting_Species(other, set())
        elif other == 0:
            raise ReactionError(
                "Use the Zero meta-species for degradation reactions, not the integer 0"
            )
        else:
            p = other

        reaction = Reactions(myself.list_of_reactants, p.list_of_reactants)
        return reaction

    def __getattr__(self, characteristic: str) -> Any:
        """Add characteristics via Species.characteristic.

        Args:
            characteristic: Characteristic to add/query.


        Returns:
            Reacting_Species with the characteristic.
        """
        if characteristic == "_ipython_canary_method_should_not_exist_":
            return 0

        if characteristic == "assign":
            asgi_Assign.set_context()
            return asgi_Asg(self, species_or_reacting=True)

        Species.check_if_valid_characteristic(self, characteristic)

        characteristics_from_references = mcu_unite_characteristics(
            list(self.get_references())
        )
        characteristics = {characteristic}

        if (
            characteristic not in characteristics_from_references
            and "$" not in characteristic
        ):
            if len(self.get_characteristics()) == 0:
                self.first_characteristic = characteristic
            self.add_characteristic(characteristic)
        return Reacting_Species(self, characteristics)

    def __call__(self, quantity: Any) -> Self | str | None:  # type: ignore[return]
        """Handle count assignment and characteristic extraction.

        Args:
            quantity: For count assignment, or Specific_Species_Operator
                for characteristic extraction.

        Returns:
            Self to allow for assigning counts mid-reaction.
        """
        if isinstance(quantity, (np_int_, np_float_)):
            quantity = float(quantity)

        quantity_dict: dict[str, Any] | None = None
        _any_ctx = Species.get_meta_specie_named_any_context()
        if len(_any_ctx) != 0:
            for i in _any_ctx:
                self.c(i)
            quantity_dict = self.add_quantities(_any_ctx.copy(), quantity)

        elif (
            isinstance(quantity, (int, float, Quantity, mp_Mobspy_Parameter))
        ) and not asgi_Assign.check_context():
            quantity_dict = self.add_quantities(STD_CHAR, quantity)

        elif asgi_Assign.check_context():
            self.assign(quantity)
        elif isinstance(quantity, me_Specific_Species_Operator):
            for cha in str(quantity).split(DOT_SEPARATOR)[1:]:
                if cha in self._characteristics:
                    return cha
            raise ReactionError(
                f"{quantity} contains no characteristics from {self._name}"
            )
        elif isinstance(quantity, Reacting_Species):
            raise ReactionError(
                "Assignments of counts using meta-species "
                "are only allowed under events in "
                "simulation context"
            )
        elif Species.get_simulation_context() is None:
            raise ReactionError(
                f"Species count assignment does not support the type {type(quantity)}"
                " if not under a simulation context"
            )

        if self.get_simulation_context() is not None:
            sim_under_context = self.get_simulation_context()

            if isinstance(quantity, str):
                quantity_dict = self.add_quantities(STD_CHAR, quantity)
            try:
                if quantity_dict is None:
                    raise ValidationError(
                        "quantity_dict is None during event context assignment"
                    )
                sim_under_context.current_event_count_data.append(
                    {
                        "species": self,
                        "characteristics": quantity_dict["characteristics"],
                        "quantity": quantity_dict["quantity"],
                    }
                )
            except (AttributeError, KeyError, TypeError, ValueError) as e:
                raise ReactionError(
                    str(e)
                    + "\n Only species count assignments are allowed in a model context"
                ) from e
        else:
            return self

    def add_quantities(
        self,
        characteristics: Any,
        quantity: Any,
    ) -> dict[str, Any] | None:
        """Set the quantity of a specific string of species.

        Args:
            characteristics: Characteristics of the species.
            quantity: Counts.
        """
        if self.get_simulation_context() is None:
            already_in = False
            for e in self._species_counts:
                if characteristics == e["characteristics"]:
                    e["quantity"] = quantity
                    already_in = True
            if not already_in:
                self._species_counts.append(
                    {"characteristics": characteristics, "quantity": quantity}
                )
        else:
            return {"characteristics": characteristics, "quantity": quantity}
        return None

    def reset_quantities(self) -> None:
        """Reset the counts inside a species."""
        self._species_counts = []

    def get_quantities(self) -> list[dict[str, Any]]:
        """Returns the list of species_counts."""
        return self._species_counts

    def __mul__(self, other: Species | Any) -> Species | Any:
        """Multiplication to construct more complex species.

        Args:
            other: For the multiplication.


        Returns:
            Higher order species from the multiplication.
        """
        if asgi_Assign.check_context():
            return asgi_Assign.mul(self, other)

        frame = sys._getframe(1)
        code_line = linecache.getline(frame.f_code.co_filename, frame.f_lineno).rstrip(
            "\n"
        )

        if not isinstance(other, Species):
            raise ReactionError(
                f"At {code_line}: \n"
                + "Meta-Species can only be multiplied by other meta-species \n"
                + f"It was multiplied by the type {type(other)}"
            )

        name = code_line.replace(" ", "").split("=")[0]

        _Last_rate_storage.increment_entity_counter()
        new_entity = Species(name)
        new_entity.set_references(mcu_combine_references(self, other))
        new_entity.add_reference(new_entity)
        new_entity._from_mul = True  # type: ignore[attr-defined]

        mcu_check_orthogonality_between_references(new_entity.get_references())

        return new_entity

    def get_spe_object(self) -> Self:
        """Return the underlying Species object."""
        return self

    def get_query_characteristics(self) -> str:
        """Return the default query characteristic string."""
        return STD_CHAR

    def link_a_species(self, other_species: Species) -> None:
        """Link a species with another.

        Args:
            other_species: Other species to be linked.
        """
        self._linked_species.add(other_species)

    def unit(self, unit: Any) -> None:
        """Set the unit for this species (no-op at base level)."""
        pass

    def name(self, name: str) -> None:
        """Name a species.

        Args:
            name: Name of the species.
        """
        if name[0] == "_":
            raise ValidationError(
                f"In species name {name}: Species cannot "
                "have a name starting with underscore"
            )

        name = clean_species_name(name)
        self._name = name

    def _bypass_name(self, name: str) -> None:
        name = clean_species_name(name)
        self._name = name

    def get_name(self) -> str:
        """Return the species name."""
        return self._name

    def get_characteristics(self) -> set[str]:
        """Return the set of characteristics for this species."""
        return self._characteristics

    def add_characteristic(self, characteristic: str) -> None:
        """Add a characteristic to this species."""
        self._characteristics.add(characteristic)

    def remove_characteristic(self, characteristic: str) -> None:
        """Remove a characteristic from this species."""
        self._characteristics.remove(characteristic)

    def print_characteristics(self) -> None:
        """Log the species characteristics at debug level."""
        _logger.debug(str(self._characteristics))

    def get_references(self) -> set[Species]:
        """Return the set of reference species."""
        return self._references

    def get_all_characteristics(self) -> set[str]:
        """Get all characteristics including from references."""
        all_char: set[str] = set()
        for reference in self._references:
            all_char = all_char.union(reference.get_characteristics())
        return all_char

    def add_reference(self, reference: Species) -> None:
        """Add a reference species to this species."""
        self._references.add(reference)

    def set_references(self, reference_set: set[Species]) -> None:
        """Replace the reference set with the given one."""
        self._references = reference_set

    def reset_references(self) -> None:
        """Reset references to only contain this species."""
        self._references = {self}

    def get_reactions(self) -> set[Reactions]:
        """Return the set of reactions involving this species."""
        return self._reactions

    def set_reactions(self, reactions: set[Reactions]) -> None:
        """Replace the reaction set with the given one."""
        self._reactions = reactions

    def reset_reactions(self) -> None:
        """Clear all reactions from this species."""
        self._reactions = set()

    def add_reaction(self, reaction: Reactions) -> None:
        """Add a reaction to this species."""
        self._reactions.add(reaction)

    def reset_counts(self) -> None:
        """Clear all initial count assignments."""
        self._species_counts = []

    @classmethod
    def set_simulation_context(cls, sim: Any) -> None:
        """Set the active simulation context for all species.

        Raises:
            ValidationError: If a context is already set.
        """
        if _simulation_context_cv.get() is None:
            _simulation_context_cv.set(sim)
        else:
            raise ValidationError(
                "A different Simulation Object was assigned "
                "to a meta-species object under context \n"
                "Please use only one Simulation Object "
                "per context assignment"
            )

    @classmethod
    def update_meta_specie_named_any_context(
        cls,
        meta_specie_named_any_characteristics: set[str],
    ) -> None:
        """Update the meta-species Any context for the current thread.

        Args:
            meta_specie_named_any_characteristics: Set of characteristics of the
                currently active any context.
        """
        _meta_any_ctx_cv.set(meta_specie_named_any_characteristics)

    @classmethod
    def reset_simulation_context(cls) -> None:
        """Clear the active simulation context."""
        _simulation_context_cv.set(None)

    @classmethod
    def get_simulation_context(cls) -> Any:
        """Return the active simulation context, or None."""
        return _simulation_context_cv.get()

    @classmethod
    def get_meta_specie_named_any_context(cls) -> set[str]:
        """Return the current Any context characteristics for this thread."""
        try:
            return _meta_any_ctx_cv.get()
        except LookupError:
            s: set[str] = set()
            _meta_any_ctx_cv.set(s)
            return s

    def order_references(self) -> None:
        """Sort references by characteristics and build an index map."""
        cleaned_references = [
            x for x in self.get_references() if x.get_characteristics() != set()
        ]
        self._ordered_references = sorted(
            cleaned_references, key=lambda x: sorted(x.get_characteristics())
        )
        i = 1
        for reference in self._ordered_references:
            self._reference_index_dictionary[reference] = i
            i = i + 1

    def get_ordered_references(self) -> list[Species]:
        """Return references sorted by characteristics."""
        return self._ordered_references

    def get_index_from_reference_dict(self, reference: Species) -> int:
        """Return the 1-based index of a reference species."""
        return self._reference_index_dictionary[reference]

    @classmethod
    def is_species(cls) -> bool:
        """Return True; this is a Species."""
        return True

    @classmethod
    def is_spe_or_reac(cls) -> bool:
        """Return True; this is a Species or Reactions type."""
        return True


def clean_species_name(species_name: str) -> str:
    """Strip tabs and spaces from a species name."""
    species_name = species_name.replace("\t", "")
    species_name = species_name.replace(" ", "")
    return species_name


_methods_Species = set(dir(Species))
