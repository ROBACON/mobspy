"""Species class and related helpers.

The core meta-species object that users interact with to define
models in MobsPy's DSL.
"""

from __future__ import annotations

import linecache
import sys
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
from mobspy.modules.declarations import CountAssignment, RatedProduct, get_registry
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
    compute_ordered_references as mcu_compute_ordered_references,
)
from mobspy.modules.species_utils import (
    unite_characteristics as mcu_unite_characteristics,
)

if TYPE_CHECKING:
    from collections.abc import Generator

    from mobspy.modules.list_species import List_Species

_logger = get_logger(__name__)


def _create_reaction_from_rated(
    reactant_list: list[dict[str, Any]],
    rated: RatedProduct,
) -> Reactions:
    """Create a Reactions object from a RatedProduct (``@`` syntax).

    Handles reversible reactions when rated.is_reversible is True.
    Validates units early when rate is a Pint Quantity.
    """
    from mobspy.modules.unit_validation import validate_rate_units  # noqa: PLC0415

    validate_rate_units(rated.rate, len(reactant_list))
    _Last_rate_storage.set_last_rate(rated.rate)
    reaction = Reactions(reactant_list, rated.products, rate=rated.rate)
    if rated.is_reversible:
        validate_rate_units(rated.reverse_rate, len(rated.products))
        Reactions(rated.products, reactant_list, rate=rated.reverse_rate)
    return reaction


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
        self._ordered_references: list[Species] | None = None
        self._reference_index_dictionary: dict[Species, int] | None = None
        self._unit: str = ""
        self._assignments: dict[str, Any] = {}
        self._linked_species: set[Species] = set()

        self.first_characteristic: str | None = None

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
        result: str = "(" + to_str + ")"
        return result

    def __str__(self) -> str:
        """String representation, returns the species name."""
        if Species.get_simulation_context() is None:
            return self._name
        return Species.str_under_context(self, STD_CHAR)

    def __repr__(self) -> str:
        return f"Species({self._name!r})"

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
        for _i, reference in enumerate(self.get_references()):
            if reference.get_characteristics():
                _logger.debug(str(reference) + ": ")
                reference.show_characteristics()

    def show_references(self) -> None:
        """Print the objects this object has inherited from."""
        _logger.debug(str(self) + DOT_SEPARATOR)
        _logger.debug("{")
        for _i, reference in enumerate(self.get_references()):
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
        from mobspy.modules.list_species import List_Species  # noqa: PLC0415

        if isinstance(other, List_Species):
            other.append(self)
            return other
        if isinstance(other, Species):
            return List_Species([self, other])
        raise ValidationError("Only Species and List_Species can be concatenated")

    def __iter__(self) -> Generator[Self, None, None]:
        """Iter defined to be consistent with List_Species behavior."""
        yield self

    def remove_repeated_elements(self) -> Self:
        """Defined to be consistent with the model behavior from List_Species."""
        return self

    def __getitem__(self, item: Any) -> Self:
        """Attach a rate via ``[]`` syntax.

        Stores the rate in a ContextVar for ``>>`` to retrieve.
        Unlike ``@``, returns ``self`` to allow stoichiometry:
        ``2 * A[rate]`` works because ``A[rate]`` returns ``A``.

        Args:
            item: Reaction rate.
        """
        return _Last_rate_storage.override_get_item(self, item)  # type: ignore[no-any-return]

    def __matmul__(self, rate: Any) -> RatedProduct:
        """Attach a rate via the ``@`` operator.

        ``B @ rate`` returns a RatedProduct consumed by ``>>``.
        A tuple ``(k_fwd, k_rev)`` creates a reversible reaction.

        Args:
            rate: Reaction rate (number, callable, tuple for reversible).
        """
        products = Reacting_Species(self, set()).list_of_reactants
        if isinstance(rate, tuple) and len(rate) == 2:  # noqa: PLR2004
            return RatedProduct(
                products=products,
                rate=rate[0],
                is_reversible=True,
                reverse_rate=rate[1],
            )
        return RatedProduct(products=products, rate=rate)

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
        return asgi_Assign.mul(stoichiometry, self)

    def __add__(
        self,
        other: Species | Reacting_Species | RatedProduct | Any,
    ) -> Reacting_Species | RatedProduct | Any:
        """Addition for reaction construction.

        When the right-hand side is a ``RatedProduct`` (from ``C @ rate``),
        this species is prepended to the product list so that
        ``B + C @ rate`` works without parentheses.

        Args:
            other: Other object added to construct a reaction.

        Returns:
            Reacting Species from the sum, or RatedProduct if other is rated.
        """
        if not asgi_Assign.check_context():
            r1 = Reacting_Species(self, set())
            if isinstance(other, RatedProduct):
                return r1 + other
            if isinstance(other, Reacting_Species):
                r2 = other
            else:
                r2 = Reacting_Species(other, set())
            return r1 + r2
        return asgi_Assign.add(self, other)

    def __radd__(self, other: Any) -> Reacting_Species | Any:
        """Making addition symmetric, see __add__."""
        if not asgi_Assign.check_context():
            return Species.__add__(self, other)
        return asgi_Assign.add(other, self)

    def __invert__(self) -> Reacting_Species:
        return self.c(NOT_CHAR)

    def __neg__(self) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.mul(-1, self)
        raise ValidationError(
            "The negative operator was applied to a Species in the wrong context"
        )

    def __rshift__(self, other: Species | Reacting_Species | RatedProduct) -> Reactions:
        """Reaction definition (``>>`` operator).

        Accepts Species, Reacting_Species, or RatedProduct (from ``@``).

        Args:
            other: Reaction products (possibly with attached rate).


        Returns:
            The reaction.
        """
        myself = Reacting_Species(self, set())

        if isinstance(other, RatedProduct):
            return _create_reaction_from_rated(myself.list_of_reactants, other)

        if isinstance(other, Species):
            p = Reacting_Species(other, set())
        elif other == 0:
            raise ReactionError(
                "Use the Zero meta-species for degradation reactions, not the integer 0"
            )
        else:
            p = other

        return Reactions(myself.list_of_reactants, p.list_of_reactants)

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

        if isinstance(quantity, me_Specific_Species_Operator):
            return self._extract_characteristic(quantity)

        quantity_dict = self._resolve_quantity(quantity)

        if self.get_simulation_context() is not None:
            self._apply_event_context(quantity, quantity_dict)
            return None
        return self

    def _extract_characteristic(self, quantity: me_Specific_Species_Operator) -> str:
        """Find and return matching characteristic, or raise."""
        for cha in str(quantity).split(DOT_SEPARATOR)[1:]:
            if cha in self._characteristics:
                return cha
        raise ReactionError(f"{quantity} contains no characteristics from {self._name}")

    def _resolve_quantity(self, quantity: Any) -> dict[str, Any] | None:
        """Resolve the quantity into a quantity_dict based on context and type."""
        _any_ctx = Species.get_meta_specie_named_any_context()
        if len(_any_ctx) != 0:
            for i in _any_ctx:
                self.c(i)
            return self.add_quantities(_any_ctx.copy(), quantity)

        if (
            isinstance(quantity, (int, float, Quantity, mp_Mobspy_Parameter))
        ) and not asgi_Assign.check_context():
            return self.add_quantities(STD_CHAR, quantity)

        if asgi_Assign.check_context():
            self.assign(quantity)
            return None
        if isinstance(quantity, Reacting_Species):
            raise ReactionError(
                "Assignments of counts using meta-species "
                "are only allowed under events in "
                "simulation context"
            )
        if Species.get_simulation_context() is None:
            raise ReactionError(
                f"Species count assignment does not support the type {type(quantity)}"
                " if not under a simulation context"
            )
        return None

    def _apply_event_context(
        self,
        quantity: Any,
        quantity_dict: dict[str, Any] | None,
    ) -> None:
        """Apply quantity assignment within a simulation event context."""
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

            frozen_chars = (
                frozenset(characteristics)
                if isinstance(characteristics, set)
                else characteristics
            )
            get_registry().add_count(
                CountAssignment(
                    species=self,
                    characteristics=frozen_chars,
                    quantity=quantity,
                )
            )
        else:
            return {"characteristics": characteristics, "quantity": quantity}
        return None

    def reset_quantities(self) -> None:
        """Reset the counts inside a species.

        Also removes counts from the global ModelRegistry.
        """
        self._species_counts = []
        get_registry().remove_counts_for(self)

    def get_quantities(self) -> list[dict[str, Any]]:
        """Returns the list of species_counts."""
        return self._species_counts

    def __mul__(self, other: Species | Any) -> Species | Any:
        """Multiplication to construct more complex species.

        Attempts to infer the variable name from the source line.
        Falls back to a generated name if inference fails.
        Use ``(A * B).named("C")`` to set the name explicitly.

        Args:
            other: For the multiplication.


        Returns:
            Higher order species from the multiplication.
        """
        if asgi_Assign.check_context():
            return asgi_Assign.mul(self, other)

        if not isinstance(other, Species):
            raise ReactionError(
                "Meta-Species can only be multiplied by other meta-species. "
                f"It was multiplied by the type {type(other)}"
            )

        name = self._infer_mul_name()
        return self._mul_impl(other, name)

    @staticmethod
    def _infer_mul_name() -> str:
        """Try to infer the variable name from the caller's source line."""
        try:
            frame = sys._getframe(2)
            code_line = linecache.getline(
                frame.f_code.co_filename, frame.f_lineno
            ).rstrip("\n")
            candidate = code_line.replace(" ", "").split("=")[0]
            if candidate and candidate[0] != "_" and "(" not in candidate:
                return candidate
        except (ValueError, OSError, IndexError):
            pass
        _Last_rate_storage.increment_entity_counter()
        return f"_Mul_{_Last_rate_storage.get_entity_counter()}"

    def _mul_impl(self, other: Species, name: str) -> Species:
        """Core logic for species multiplication (shared by __mul__ and named).

        Orthogonality is checked eagerly at definition time.
        Ordered references are computed lazily (characteristics
        may be added to parents after multiplication).
        """
        _Last_rate_storage.increment_entity_counter()
        counter = _Last_rate_storage.get_entity_counter()
        new_entity = Species(f"Mul{counter}")
        new_entity._bypass_name(name)

        combined = mcu_combine_references(self, other)
        combined.add(new_entity)
        new_entity._references = combined
        new_entity._from_mul = True  # type: ignore[attr-defined]

        # Fail fast: check orthogonality at definition time
        mcu_check_orthogonality_between_references(combined)

        return new_entity

    def freeze_references(self) -> None:
        """Eagerly compute and freeze ordered references.

        Called during compilation, after all characteristics have
        been added. After this call, the ordered references and
        index dict are immutable.
        """
        ordered, index_dict = mcu_compute_ordered_references(self)
        self._ordered_references = ordered
        self._reference_index_dictionary = index_dict

    def named(self, name: str) -> Species:
        """Set or rename this species.

        Primarily used after ``*`` to avoid frame introspection::

            Thing = (Color * Size).named("Thing")

        Args:
            name: The species name.
        """
        self.name(name)
        return self

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

    def _invalidate_ordered_references(self) -> None:
        """Invalidate the cached ordered references."""
        self._ordered_references = None
        self._reference_index_dictionary = None

    def add_reference(self, reference: Species) -> None:
        """Add a reference species to this species."""
        self._references.add(reference)
        self._invalidate_ordered_references()

    def set_references(self, reference_set: set[Species]) -> None:
        """Replace the reference set with the given one."""
        self._references = reference_set
        self._invalidate_ordered_references()

    def reset_references(self) -> None:
        """Reset references to only contain this species."""
        self._references = {self}
        self._invalidate_ordered_references()

    def get_reactions(self) -> set[Reactions]:
        """Return the set of reactions involving this species.

        Queries the :class:`~mobspy.modules.declarations.ModelRegistry`.
        """
        registry = get_registry()
        return registry.reactions_for_species(frozenset({id(self)}))

    def reset_reactions(self) -> None:
        """Remove all reactions involving this species from the registry.

        Subsequent Simulations will not pick up the removed reactions.
        """
        registry = get_registry()
        registry.remove_reactions_for(self)

    def reset_counts(self) -> None:
        """Clear all initial count assignments."""
        self._species_counts = []

    @classmethod
    def set_simulation_context(cls, sim: Any) -> None:
        """Set the active simulation context for all species.

        Raises:
            ValidationError: If a context is already set.
        """
        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        session = get_session()
        if session.simulation_context is None:
            session.simulation_context = sim
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
        """Update the meta-species Any context for the current thread."""
        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        get_session().meta_any_ctx = meta_specie_named_any_characteristics

    @classmethod
    def reset_simulation_context(cls) -> None:
        """Clear the active simulation context."""
        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        get_session().simulation_context = None

    @classmethod
    def get_simulation_context(cls) -> Any:
        """Return the active simulation context, or None."""
        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        return get_session().simulation_context

    @classmethod
    def get_meta_specie_named_any_context(cls) -> set[str]:
        """Return the current Any context characteristics for this thread."""
        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        return get_session().meta_any_ctx

    def _ensure_ordered_references(self) -> None:
        """Compute ordered references and index map if not already frozen.

        Normally called lazily on first access. After
        ``freeze_references()`` (called during compilation), this
        is a no-op.
        """
        if self._ordered_references is None:
            ordered, index_dict = mcu_compute_ordered_references(self)
            self._ordered_references = ordered
            self._reference_index_dictionary = index_dict

    def order_references(self) -> None:
        """Sort references by characteristics and build an index map.

        Now delegates to the lazy computation. Calling this explicitly
        is no longer required but remains safe.
        """
        self._ordered_references = None
        self._reference_index_dictionary = None
        self._ensure_ordered_references()

    def get_ordered_references(self) -> list[Species]:
        """Return references sorted by characteristics."""
        self._ensure_ordered_references()
        return self._ordered_references  # type: ignore[return-value]

    def get_index_from_reference_dict(self, reference: Species) -> int:
        """Return the 1-based index of a reference species."""
        self._ensure_ordered_references()
        return self._reference_index_dictionary[reference]  # type: ignore[index]

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
    return species_name.replace(" ", "")


_methods_Species = set(dir(Species))  # noqa: N816
