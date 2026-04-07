"""Reaction-related classes: Reactions, Reacting_Species, and helpers.

Contains _Last_rate_storage, Assignment_Opp_Imp, Reactions, and
Reacting_Species. These are tightly coupled with Species but separated
here for organizational clarity.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Any, Self

from numpy import floating as np_float_
from numpy import integer as np_int_
from pint import Quantity

from mobspy.constants import CONTEXT_SPECIES_NAME, NOT_CHAR, ZERO_SPECIES_NAME
from mobspy.exceptions import ReactionError, ValidationError
from mobspy.modules.assignments_implementation import (
    Asg as asgi_Asg,
)
from mobspy.modules.assignments_implementation import (
    Assign as asgi_Assign,
)
from mobspy.modules.declarations import (
    RatedProduct,
    ReactantRef,
    ReactionDecl,
    get_registry,
)
from mobspy.modules.logic_operators import (
    ReactingSpeciesComparator as lop_ReactingSpeciesComparator,
)
from mobspy.modules.mobspy_expressions import (
    ExpressionDefiner as me_ExpressionDefiner,
)
from mobspy.modules.mobspy_expressions import (
    OverrideQuantity as me_OverrideQuantity,
)
from mobspy.modules.mobspy_parameters import (
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)
from mobspy.modules.species_utils import (
    unite_characteristics as mcu_unite_characteristics,
)

if TYPE_CHECKING:
    from mobspy.modules.declarations import RateValue
    from mobspy.modules.mobspy_expressions import MobsPyExpression
    from mobspy.modules.species import Species


class _Last_rate_storage:  # noqa: N801
    """Legacy rate buffering for the ``[]`` bracket syntax.

    When the user writes ``A[rate] >> B``, Python evaluates
    ``A[rate]`` first (via ``__getitem__``), which stores the rate
    in a thread-local ``ContextVar``.  Then ``A >> B`` (via
    ``__rshift__``) retrieves it.

    **This mechanism is not needed for the ``@`` rate syntax.**
    ``A >> B @ rate`` produces a ``RatedProduct`` that carries the
    rate directly, bypassing this class entirely.

    The entity counter is unrelated to rate storage and is used
    to generate unique names for auto-created species.
    """

    @staticmethod
    def get_last_rate() -> RateValue:
        """Return the stored rate for the current thread."""
        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        return get_session().last_rate  # type: ignore[no-any-return]

    @staticmethod
    def set_last_rate(value: RateValue) -> None:
        """Store a rate value for the current thread."""
        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        get_session().last_rate = value

    @staticmethod
    def get_entity_counter() -> int:
        """Return the entity counter for the current thread."""
        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        return get_session().entity_counter

    @staticmethod
    def increment_entity_counter() -> None:
        """Increment the entity counter for the current thread."""
        from mobspy.modules.session_context import get_session  # noqa: PLC0415

        get_session().entity_counter += 1

    @classmethod
    def override_get_item(
        cls,
        object_to_return: Reactions | Reacting_Species,
        item: RateValue,
    ) -> Reactions | Reacting_Species:
        """Store the rate before the reaction is fully defined.

        .. deprecated::
            Use ``A >> B @ rate`` instead of ``A >> B[rate]``.

        Args:
            object_to_return: Returns the object that __getitem__ was performed on.
            item: Stored reaction rate.
        """
        import warnings  # noqa: PLC0415

        warnings.warn(
            "A >> B[rate] is deprecated. Use 'A >> B @ rate' instead.",
            DeprecationWarning,
            stacklevel=3,
        )
        cls.set_last_rate(item)
        return object_to_return

    @classmethod
    def process_rate(cls, rate: RateValue) -> RateValue:
        """Validate and normalize a reaction rate value."""
        from mobspy.modules.rate_builder import RateExpression  # noqa: PLC0415
        from mobspy.modules.species import Species  # noqa: PLC0415

        if isinstance(rate, RateExpression):
            return rate.node

        if isinstance(rate, (np_int_, np_float_)):
            rate = float(rate)

        if isinstance(rate, (Species, Reacting_Species, Reactions)):
            raise ReactionError(
                f"Reaction rate of type {type(cls.get_last_rate())} not valid"
            )

        if not (
            isinstance(rate, (int, float, str))
            or callable(rate)
            or isinstance(
                rate,
                (
                    me_OverrideQuantity,
                    Quantity,
                    mp_Mobspy_Parameter,
                    me_ExpressionDefiner,
                ),
            )
            or rate is None
        ):
            raise ReactionError(
                "Reaction rate of type " + str(type(rate)) + " not valid"
            )

        return rate


def _apply_any_context(
    reactants: list[dict[str, Any]],
    products: list[dict[str, Any]],
) -> None:
    """Apply Any meta-species context characteristics to reactants and products."""
    from mobspy.modules.species import Species  # noqa: PLC0415

    _any_ctx = Species.get_meta_specie_named_any_context()
    if len(_any_ctx) != 0:
        for j in _any_ctx:
            for r in reactants:
                r["object"].c(j)
                r["characteristics"].add(j)
            for p in products:
                p["object"].c(j)
                p["characteristics"].add(j)


class Reactions:
    """Reaction class storing reactants, products, rate and order.

    Reactions are created with the ``>>`` operator and stored
    in all involved objects.

    Args:
        reactants: List of meta-species used as reactants.
        products: List of meta-species used as products.
        order: Reaction order operator (Round-Robin default).
        rate: Reaction rate.
    """

    def __init__(
        self,
        reactants: list[dict[str, Any]],
        products: list[dict[str, Any]],
        rate: RateValue = None,
    ) -> None:
        """Construct a reaction from reactants and products.

        The order and the rate are assigned later by the compiler.

        Args:
            reactants: List of meta-species reactants.
            products: List of meta-species products.
        """
        from mobspy.modules.species import Species  # noqa: PLC0415

        self._validate_reaction_context(products)

        self.rate = self._process_rate_assignment(reactants, products, rate)

        _apply_any_context(reactants, products)

        try:
            _ = reactants[0]["object"]
        except IndexError:
            try:
                _ = products[0]["object"]
            except IndexError as e:
                raise ReactionError("No Meta-Species detected in the reaction") from e

        if Species.get_simulation_context() is not None:
            raise ReactionError(
                "Reactions cannot be defined under event context. Only species counts"
            )

        self.reactants = reactants
        self.products = products

        self.order = None

        self._register_declaration()

    def _register_declaration(self) -> None:
        """Register this reaction as a ReactionDecl in the ModelRegistry."""
        registry = get_registry()
        decl = ReactionDecl(
            reactants=tuple(
                ReactantRef(
                    species=r["object"],
                    characteristics=frozenset(r["characteristics"]),
                    stoichiometry=r["stoichiometry"],
                    label=r.get("label"),
                )
                for r in self.reactants
            ),
            products=tuple(
                ReactantRef(
                    species=p["object"],
                    characteristics=frozenset(p["characteristics"]),
                    stoichiometry=p["stoichiometry"],
                    label=p.get("label"),
                )
                for p in self.products
            ),
            rate=self.rate,
        )
        registry.add_reaction(decl, reaction_obj=self)

    @staticmethod
    def _validate_reaction_context(products: list[dict[str, Any]]) -> None:
        """Validate that no forbidden context is active and products are valid."""
        for p in products:
            if p["object"].get_name() == CONTEXT_SPECIES_NAME:
                raise ReactionError("The Any meta-species cannot be used in reactions")

        if asgi_Assign.check_context():
            asgi_Assign.reset_context()
            raise ReactionError(
                "A MobsPy context error has happened. "
                "A reaction was defined with the assignment context activated. "
                "The assignment context was deactivated. "
                "Please try to redefine the model"
            )

    def _process_rate_assignment(
        self,
        reactants: list[dict[str, Any]],
        products: list[dict[str, Any]],
        rate: RateValue,
    ) -> RateValue:
        """Resolve the reaction rate from explicit arg or ContextVar.

        The ``@`` path passes rate explicitly. The ``[]`` path
        stores the rate in ``_Last_rate_storage`` and passes
        ``rate=None`` here. Handles reversible tuple rates from
        ``Rev[...][k1, k2]``.
        """
        if rate is not None:
            return _Last_rate_storage.process_rate(rate)

        stored = _Last_rate_storage.get_last_rate()

        # Reversible tuple from Rev[...][k1, k2]
        import contextlib  # noqa: PLC0415

        is_tuple = False
        with contextlib.suppress(TypeError, AttributeError):
            is_tuple = isinstance(stored, tuple) and len(stored) == 2  # noqa: PLR2004

        if is_tuple and isinstance(stored, tuple):
            fwd_rate = stored[0]
            Reactions(reactants=products, products=reactants, rate=stored[1])
            resolved = _Last_rate_storage.process_rate(fwd_rate)
        else:
            resolved = _Last_rate_storage.process_rate(stored)

        _Last_rate_storage.set_last_rate(None)
        return resolved

    @staticmethod
    def __create_reactants_string(list_of_reactants: list[dict[str, Any]]) -> str:
        """Format a list of reactants/products as a string."""
        reaction_string = ""
        for i, r in enumerate(list_of_reactants):
            if r["stoichiometry"] > 1:
                reaction_string += str(r["stoichiometry"]) + "*" + str(r["object"])
            else:
                reaction_string += str(r["object"])

            if len(r["characteristics"]) > 0:
                reaction_string += "." + ".".join(r["characteristics"])

            if i != len(list_of_reactants) - 1:
                reaction_string += " "
                reaction_string += "+"
                reaction_string += " "

        return reaction_string

    def __str__(self) -> str:
        """Print the meta-reaction in the format A + B -> C + D."""
        return (
            self.__create_reactants_string(self.reactants)
            + " -> "
            + self.__create_reactants_string(self.products)
        )

    def __getitem__(self, item: RateValue) -> Self:
        """Attach a rate via ``[]`` syntax.

        Args:
            item: Reaction rate.
        """
        return _Last_rate_storage.override_get_item(self, item)  # type: ignore[return-value]

    def __matmul__(self, rate: RateValue | tuple[RateValue, RateValue]) -> Reactions:
        """Attach or update a rate via the ``@`` operator.

        ``(A >> B) @ rate`` sets the rate on an existing reaction.
        A tuple ``(k_fwd, k_rev)`` creates the reverse reaction too.

        Args:
            rate: Reaction rate (number, callable, tuple for reversible).
        """
        if isinstance(rate, tuple) and len(rate) == 2:  # noqa: PLR2004
            self.rate = _Last_rate_storage.process_rate(rate[0])
            Reactions(self.products, self.reactants, rate=rate[1])
        else:
            self.rate = _Last_rate_storage.process_rate(rate)
        return self

    def set_rate(self, rate: RateValue) -> None:
        """Set the stored reaction rate.

        Args:
            rate: Reaction rate.
        """
        self.rate = rate


class Assignment_Opp_Imp:  # noqa: N801
    """Mixin that routes arithmetic operators to the assignment context when active."""

    @staticmethod
    def _dispatch_assign_op(
        first: Any,
        second: Any,
        op_func: Callable[..., MobsPyExpression],
        error_verb: str,
    ) -> MobsPyExpression:
        """Dispatch an arithmetic operator to the assignment or expression context.

        Args:
            first: left operand (already ordered by caller)
            second: right operand (already ordered by caller)
            op_func: the ``asgi_Assign`` method to call
            error_verb: verb for the error message (e.g. "Addition")
        """
        if asgi_Assign.check_context():
            return op_func(first, second)
        # Outside assignment context: wrap as rate expressions.
        # This enables both lambda replay (expression mode) and
        # standalone rate building (e.g. Protein ** 2 / (Km + Protein))
        from mobspy.modules.mobspy_expressions import (  # noqa: PLC0415
            check_if_non_expression_operated,
        )
        from mobspy.modules.species_operators import _ms_active_ctx  # noqa: PLC0415

        if _ms_active_ctx.get():
            # Inside lambda replay: wrap as MobsPyExpression
            wrapped_first = check_if_non_expression_operated(first)
            wrapped_second = check_if_non_expression_operated(second)
            op_map = {
                "Addition": "__add__",
                "Subtraction": "__sub__",
                "Multiplication": "__mul__",
                "Division": "__truediv__",
                "Exponentiation": "__pow__",
            }
            py_op = op_map.get(error_verb, "__mul__")
            return getattr(wrapped_first, py_op)(wrapped_second)  # type: ignore[no-any-return]

        # Outside any context: wrap as RateExpression for rate builder use
        from mobspy.modules.expression_nodes import (  # noqa: PLC0415
            BinaryOpNode,
            LiteralNode,
            SpeciesRefNode,
        )
        from mobspy.modules.rate_builder import RateExpression  # noqa: PLC0415
        from mobspy.modules.species import Species as _Spe  # noqa: PLC0415

        def _to_node(x: Any) -> Any:
            if isinstance(x, _Spe):
                return SpeciesRefNode(x.get_name())
            if isinstance(x, RateExpression):
                return x.node
            if isinstance(x, (int, float)):
                return LiteralNode(x)
            return x

        sbml_op = {
            "Addition": "+",
            "Subtraction": "-",
            "Multiplication": "*",
            "Division": "/",
            "Exponentiation": "^",
        }.get(error_verb, "*")
        return RateExpression(BinaryOpNode(_to_node(first), sbml_op, _to_node(second)))  # type: ignore[return-value]

    def __add__(self, other: object) -> Any:
        return self._dispatch_assign_op(self, other, asgi_Assign.add, "Addition")

    def __radd__(self, other: object) -> Any:
        return self._dispatch_assign_op(other, self, asgi_Assign.add, "Addition")

    def __sub__(self, other: object) -> Any:
        return self._dispatch_assign_op(self, other, asgi_Assign.sub, "Subtraction")

    def __rsub__(self, other: object) -> Any:
        return self._dispatch_assign_op(other, self, asgi_Assign.sub, "Subtraction")

    def __truediv__(self, other: object) -> Any:
        return self._dispatch_assign_op(self, other, asgi_Assign.div, "Division")

    def __rtruediv__(self, other: object) -> Any:
        return self._dispatch_assign_op(other, self, asgi_Assign.div, "Division")

    def __pow__(self, other: object) -> Any:
        return self._dispatch_assign_op(self, other, asgi_Assign.pow, "Exponentiation")

    def __rpow__(self, other: object) -> Any:
        return self._dispatch_assign_op(other, self, asgi_Assign.pow, "Exponentiation")

    def __mul__(self, other: object) -> Any:
        return self._dispatch_assign_op(self, other, asgi_Assign.mul, "Multiplication")

    def __rmul__(self, other: object) -> Any:
        return self._dispatch_assign_op(other, self, asgi_Assign.mul, "Multiplication")


class Reacting_Species(lop_ReactingSpeciesComparator, Assignment_Opp_Imp):  # noqa: N801
    """Intermediary object created when a species enters a reaction.

    Transforms a Species object into a list-compatible format
    for the reaction object. The ``>>`` operator calls the
    Reaction constructor.

    Attributes:
        list_of_reactants: Represents reactants or products.
            Each dict contains 'object', 'characteristics',
            'stoichiometry', and optionally 'label'.
    """

    def __init__(
        self,
        object_reference: Species,
        characteristics: set[str],
        stoichiometry: int | float = 1,
        label: int | float | str | None = None,
    ) -> None:
        """Construct a Reacting_Species.

        Args:
            object_reference: Meta-species object reference.
            characteristics: Characteristics used to query.
            stoichiometry: Stoichiometry in the reaction.
            label: Label value for matching.
        """
        super().__init__()
        is_zero = object_reference.get_name() == ZERO_SPECIES_NAME
        if is_zero and characteristics == set():
            self.list_of_reactants: list[dict[str, Any]] = []
        else:
            self.list_of_reactants = [
                {
                    "object": object_reference,
                    "characteristics": characteristics,
                    "stoichiometry": stoichiometry,
                    "label": label,
                }
            ]

    def __enter__(self) -> int:
        """Enter context manager for characteristics."""
        self.context_initiator_for_reacting_specie()
        return 0

    def __exit__(self, *args: object) -> None:
        """Exit context manager for characteristics."""
        self.context_finish_for_reacting_specie()

    def __str__(self) -> str:
        """String representation of the list of reactants."""
        from mobspy.modules.species import Species  # noqa: PLC0415

        species_object = self.list_of_reactants[0]["object"]
        characteristics = self.list_of_reactants[0]["characteristics"]
        if len(self.list_of_reactants) == 1:
            if Species.get_simulation_context() is None:
                to_return = str(species_object)
                for cha in self.list_of_reactants[0]["characteristics"]:
                    to_return += "." + cha
                return to_return
            return Species.str_under_context(species_object, characteristics)
        if Species.get_simulation_context() is not None:
            raise ReactionError(
                "Please separate the species when using "
                "string based assignments under event "
                "context. Ex: str(A) + str(B)"
            )
        return str(self.list_of_reactants)

    def c(self, item: Any) -> Reacting_Species:
        """Query by value instead of name.

        Calls __getattr__ with the value inside the variable.

        Args:
            item: Value to query over.
        """
        from mobspy.modules.species import Species  # noqa: PLC0415

        item = str(item)
        Species.check_if_valid_characteristic(self, item)
        return self.__getattr__(item)  # type: ignore[no-any-return]

    def label(self, label: int | float | str) -> Self:
        """Assign a label to a meta-species for compiler matching.

        Args:
            label: Value for the label for matching.
        """
        if len(self.list_of_reactants) == 1:
            self.list_of_reactants[0]["label"] = label
        else:
            raise ReactionError(
                "Labels cannot be assigned to multiple "
                "reacting species at the same time."
            )
        return self

    def __getitem__(self, item: RateValue) -> Self:
        """Attach a rate via ``[]`` syntax.

        Stores the rate in a ContextVar for ``>>`` to retrieve.
        Unlike ``@``, returns ``self`` to allow stoichiometry:
        ``2 * A.x[rate]`` works because ``A.x[rate]`` returns ``A.x``.

        Args:
            item: Reaction rate.
        """
        return _Last_rate_storage.override_get_item(self, item)  # type: ignore[return-value]

    def __matmul__(self, rate: RateValue | tuple[RateValue, RateValue]) -> RatedProduct:
        """Attach a rate via the ``@`` operator.

        ``B @ rate`` returns a RatedProduct consumed by ``>>``.
        A tuple ``(k_fwd, k_rev)`` creates a reversible reaction.

        Args:
            rate: Reaction rate (number, callable, tuple for reversible).
        """
        if isinstance(rate, tuple) and len(rate) == 2:  # noqa: PLR2004
            return RatedProduct(
                products=list(self.list_of_reactants),
                rate=rate[0],
                is_reversible=True,
                reverse_rate=rate[1],
            )
        return RatedProduct(products=list(self.list_of_reactants), rate=rate)

    def get_spe_object(self) -> Species:
        """Return the single underlying Species object.

        Raises:
            ReactionError: If more than one reactant is present.
        """
        if len(self.list_of_reactants) != 1:
            raise ReactionError(
                "The internal method get_queried_characteristics can only be used for "
                "Reacting_Species with a single reactant."
            )
        return self.list_of_reactants[0]["object"]  # type: ignore[no-any-return]

    def get_query_characteristics(self) -> set[str]:
        """Return the queried characteristics of the single reactant.

        Raises:
            ReactionError: If more than one reactant is present.
        """
        if len(self.list_of_reactants) != 1:
            raise ReactionError(
                "The internal method get_queried_characteristics can only be used for "
                "Reacting_Species with a single reactant."
            )
        return self.list_of_reactants[0]["characteristics"]  # type: ignore[no-any-return]

    def __rmul__(self, stoichiometry: Any) -> Self | MobsPyExpression:
        """Multiply by stoichiometry for reactions.

        Args:
            stoichiometry: Stoichiometry value.
        """
        if not asgi_Assign.check_context():
            if isinstance(stoichiometry, (int, float)):
                self.list_of_reactants[0]["stoichiometry"] = stoichiometry
            else:
                raise ReactionError(
                    "Stoichiometry can only be an int or "
                    f"float - Received {stoichiometry}"
                )
            return self
        return asgi_Assign.mul(stoichiometry, self)

    def __add__(
        self,
        other: Species | Reacting_Species | RatedProduct | Any,
    ) -> Self | RatedProduct | MobsPyExpression:
        """Addition of meta-species to construct the reaction.

        When the right-hand side is a ``RatedProduct`` (from ``C @ rate``),
        this species is prepended to the product list so that
        ``B + C @ rate`` works without parentheses.

        Args:
            other: Other object being added.
        """
        from mobspy.modules.species import Species  # noqa: PLC0415

        if not asgi_Assign.check_context():
            if isinstance(other, RatedProduct):
                merged = list(self.list_of_reactants) + other.products
                return RatedProduct(
                    products=merged,
                    rate=other.rate,
                    is_reversible=other.is_reversible,
                    reverse_rate=other.reverse_rate,
                )
            if isinstance(other, Species):
                other = Reacting_Species(other, set())
            try:
                self.list_of_reactants += other.list_of_reactants
            except AttributeError as e:
                raise ReactionError(
                    "Addition between meta-species and "
                    f"types {type(other)} is not supported"
                ) from e
            return self
        return asgi_Assign.add(self, other)

    def __radd__(self, other: Any) -> Self | RatedProduct | MobsPyExpression:
        if not asgi_Assign.check_context():
            return Reacting_Species.__add__(self, other)  # pyright: ignore[reportReturnType]
        return asgi_Assign.add(other, self)

    def __invert__(self) -> Reacting_Species:
        return self.c(NOT_CHAR)

    def __neg__(self) -> MobsPyExpression:
        if asgi_Assign.check_context():
            return asgi_Assign.mul(-1, self)
        raise ValidationError(
            "The negative operator was applied to a "
            "Reacting Species in the wrong context"
        )

    def __rshift__(self, other: Species | Reacting_Species | RatedProduct) -> Reactions:
        """The ``>>`` operator for defining reactions.

        Accepts Species, Reacting_Species, or RatedProduct (from ``@``).

        Args:
            other: Product side of the reaction being added.
        """
        if isinstance(other, RatedProduct):
            from mobspy.modules.species import (  # noqa: PLC0415
                _create_reaction_from_rated,
            )

            return _create_reaction_from_rated(self.list_of_reactants, other)

        from mobspy.modules.species import Species  # noqa: PLC0415

        p = Reacting_Species(other, set()) if isinstance(other, Species) else other

        return Reactions(self.list_of_reactants, p.list_of_reactants)

    def __call__(  # type: ignore[return]
        self,
        quantity: int | float | str | Quantity | mp_Mobspy_Parameter,
    ) -> Self | None:
        """Assign counts to species non-default state.

        Args:
            quantity: Count to be assigned.
        """
        from mobspy.modules.species import Species  # noqa: PLC0415

        if isinstance(quantity, (np_int_, np_float_)):
            quantity = float(quantity)

        _any_ctx = Species.get_meta_specie_named_any_context()
        if len(_any_ctx) > 0:
            for i in _any_ctx:
                self = self.c(i)  # type: ignore[assignment]  # noqa: PLW0642

        species_object = self.list_of_reactants[0]["object"]
        characteristics = self.list_of_reactants[0]["characteristics"]
        simulation_under_context = self.list_of_reactants[0][
            "object"
        ].get_simulation_context()
        quantity_dict: dict[str, Any] | None = None
        if (
            isinstance(quantity, (int, float, Quantity, mp_Mobspy_Parameter))
        ) and not asgi_Assign.check_context():
            if len(self.list_of_reactants) != 1:
                raise ReactionError(
                    "Assignment used incorrectly. Only one species at a time"
                )
            quantity_dict = species_object.add_quantities(characteristics, quantity)
        elif asgi_Assign.check_context():
            dummy_rsp = species_object
            dummy_rsp = [dummy_rsp.c(char) for char in characteristics][-1]
            dummy_rsp.assign(quantity)
        elif simulation_under_context is None:
            raise ReactionError(
                "Reactant_Species count assignment does "
                f"not support the type {type(quantity)}"
            )

        if simulation_under_context is not None:
            try:
                if isinstance(quantity, str):
                    quantity_dict = species_object.add_quantities(
                        characteristics, quantity
                    )
                if quantity_dict is None:
                    raise ReactionError(
                        "Could not resolve quantity for event assignment"
                    )
                simulation_under_context.current_event_count_data.append(
                    {
                        "species": species_object,
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
            return self  # pyright: ignore[reportReturnType]

    @property
    def ref(self) -> Any:
        """Return a rate expression referencing this species state's concentration."""
        from mobspy.modules.expression_nodes import SpeciesRefNode  # noqa: PLC0415
        from mobspy.modules.rate_builder import RateExpression  # noqa: PLC0415

        return RateExpression(SpeciesRefNode(str(self)))

    def __getattr__(self, characteristic: str) -> Any:
        """Implementation of the .dot operation.

        Args:
            characteristic: Characteristic for the query.
        """
        from mobspy.modules.species import Species  # noqa: PLC0415

        if characteristic == "_ipython_canary_method_should_not_exist_":
            return 0

        Species.check_if_valid_characteristic(self, characteristic)

        if characteristic == "assign":
            return asgi_Asg(self, species_or_reacting=False)

        for reactant in self.list_of_reactants:
            species_object = reactant["object"]
            characteristics_from_references = mcu_unite_characteristics(
                species_object.get_references()
            )

            if (
                characteristic not in characteristics_from_references
                and "$" not in characteristic
            ):
                if len(species_object.get_characteristics()) == 0:
                    species_object.first_characteristic = characteristic

                species_object.add_characteristic(characteristic)

            reactant["characteristics"].add(characteristic)

        return self

    @classmethod
    def is_species(cls) -> bool:
        """Return False; this is a reacting species, not a Species."""
        return False

    @classmethod
    def is_spe_or_reac(cls) -> bool:
        """Return True; this is a species-or-reaction type."""
        return True

    def context_initiator_for_reacting_specie(self) -> None:
        """Add the current context and update the Cts context in all meta-species."""
        from mobspy.modules.species import Species  # noqa: PLC0415

        if len(self.list_of_reactants) == 1:
            self._old_context = Species.get_meta_specie_named_any_context()
            new_context = Species.get_meta_specie_named_any_context().union(
                self.list_of_reactants[0]["characteristics"]
            )
            Species.update_meta_specie_named_any_context(new_context)
        else:
            raise ReactionError(
                "Contexts can only be used on basic Reacting meta species"
            )

    def context_finish_for_reacting_specie(self) -> None:
        """Remove the ending context and update."""
        from mobspy.modules.species import Species  # noqa: PLC0415

        Species.update_meta_specie_named_any_context(self._old_context)


_methods_Reacting_Species = set(dir(Reacting_Species))  # noqa: N816
