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
from mobspy.dsl.assignments_implementation import (
    Asg as asgi_Asg,
)
from mobspy.dsl.assignments_implementation import (
    Assign as asgi_Assign,
)
from mobspy.dsl.declarations import (
    RatedProduct,
    ReactantRef,
    ReactionDecl,
    get_registry,
)
from mobspy.dsl.logic_operators import (
    ReactingSpeciesComparator as lop_ReactingSpeciesComparator,
)
from mobspy.dsl.mobspy_parameters import (
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)
from mobspy.dsl.species_utils import (
    unite_characteristics as mcu_unite_characteristics,
)
from mobspy.exceptions import ReactionError
from mobspy.expressions.evaluation import (
    ExpressionDefiner as me_ExpressionDefiner,
)
from mobspy.expressions.evaluation import (
    OverrideQuantity as me_OverrideQuantity,
)

if TYPE_CHECKING:
    from mobspy.dsl.declarations import RateValue
    from mobspy.dsl.species import Species
    from mobspy.expressions.evaluation import MobsPyExpression


class _Last_rate_storage:  # noqa: N801  # legacy DSL public API name
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
        from mobspy.dsl.session_context import get_session  # noqa: PLC0415

        return get_session().last_rate  # type: ignore[no-any-return]

    @staticmethod
    def set_last_rate(value: RateValue) -> None:
        """Store a rate value for the current thread."""
        from mobspy.dsl.session_context import get_session  # noqa: PLC0415

        get_session().last_rate = value

    @staticmethod
    def get_entity_counter() -> int:
        """Return the entity counter for the current thread."""
        from mobspy.dsl.session_context import get_session  # noqa: PLC0415

        return get_session().entity_counter

    @staticmethod
    def increment_entity_counter() -> None:
        """Increment the entity counter for the current thread."""
        from mobspy.dsl.session_context import get_session  # noqa: PLC0415

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
        import warnings  # noqa: PLC0415  # circular import avoidance

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
        from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import
        from mobspy.expressions.rate_builder import RateExpression  # noqa: PLC0415

        if isinstance(rate, RateExpression):
            return rate.node

        if isinstance(rate, (np_int_, np_float_)):
            rate = float(rate)

        if isinstance(rate, (Species, Reacting_Species, Reactions)):
            raise ReactionError(f"Reaction rate of type {type(rate)} not valid")

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
    from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import avoidance

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
        from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import

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
        import contextlib  # noqa: PLC0415  # circular import avoidance

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
        return _Last_rate_storage.override_get_item(self, item)  # type: ignore[return-value]  # pint Unit subtype

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


# Forward/reflected dunder pairs per arithmetic verb, used to build rate
# expressions during lambda replay. The reflected member handles the case
# where the left operand is a plain number (``number.__op__(expr)`` returns
# ``NotImplemented``), e.g. ``1 + T`` or ``2 * T``.
_RATE_OP_DUNDERS: dict[str, tuple[str, str]] = {
    "Addition": ("__add__", "__radd__"),
    "Subtraction": ("__sub__", "__rsub__"),
    "Multiplication": ("__mul__", "__rmul__"),
    "Division": ("__truediv__", "__rtruediv__"),
    "Exponentiation": ("__pow__", "__rpow__"),
}


class Assignment_Opp_Imp:  # noqa: N801  # legacy DSL public API name
    """Mixin that routes arithmetic operators to the assignment context when active."""

    @staticmethod
    def _rate_expression_op(
        first: Any, second: Any, error_verb: str
    ) -> MobsPyExpression:
        """Build the rate-expression term for ``first <op> second``.

        Wraps bare species as ``MobsPyExpression`` and applies the operator.
        When the left operand is a plain number, the forward dunder yields
        ``NotImplemented``; the reflected dunder is then used, which preserves
        operand order for non-commutative operators.
        """
        from mobspy.expressions.evaluation import (  # noqa: PLC0415  # circular import
            check_if_non_expression_operated,
        )

        wrapped_first = check_if_non_expression_operated(first)
        wrapped_second = check_if_non_expression_operated(second)
        forward, reflected = _RATE_OP_DUNDERS.get(error_verb, ("__mul__", "__rmul__"))
        result = getattr(wrapped_first, forward)(wrapped_second)
        if result is NotImplemented:
            result = getattr(wrapped_second, reflected)(wrapped_first)
        return result  # type: ignore[no-any-return]  # dynamic dispatch

    @staticmethod
    def _rate_lambda_active() -> bool:
        """True when a rate-function lambda is being replayed (expression mode).

        In this mode the species/reaction operator overloads must build rate
        expressions instead of constructing reactions.
        """
        from mobspy.dsl.species_operators import _ms_active_ctx  # noqa: PLC0415

        return bool(_ms_active_ctx.get())

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
        from mobspy.dsl.species_operators import _ms_active_ctx  # noqa: PLC0415

        if _ms_active_ctx.get():
            # Inside lambda replay: wrap as MobsPyExpression
            return Assignment_Opp_Imp._rate_expression_op(first, second, error_verb)

        # Outside any context: wrap as RateExpression for rate builder use
        from mobspy.expressions.nodes import BinaryOpNode  # noqa: PLC0415
        from mobspy.expressions.rate_builder import RateExpression  # noqa: PLC0415

        sbml_op = {
            "Addition": "+",
            "Subtraction": "-",
            "Multiplication": "*",
            "Division": "/",
            "Exponentiation": "^",
        }.get(error_verb, "*")
        return RateExpression(  # type: ignore[return-value]  # pint Unit subtype
            BinaryOpNode(
                Assignment_Opp_Imp._operand_to_node(first),
                sbml_op,
                Assignment_Opp_Imp._operand_to_node(second),
            )
        )

    @staticmethod
    def _operand_to_node(x: Any) -> Any:
        """Convert an operand to a rate-expression AST node for standalone building."""
        from mobspy.dsl.species import Species as _Spe  # noqa: PLC0415
        from mobspy.expressions.nodes import (  # noqa: PLC0415  # circular import
            LiteralNode,
            ParamRefNode,
            SpeciesRefNode,
        )
        from mobspy.expressions.rate_builder import RateExpression  # noqa: PLC0415

        if isinstance(x, _Spe):
            return SpeciesRefNode(x.get_name())
        if isinstance(x, Reacting_Species):
            # A queried species (e.g. ``T.x``) renders via its dot string,
            # matching the lambda-replay path's check_if_non_expression_operated.
            return SpeciesRefNode(str(x))
        if isinstance(x, RateExpression):
            return x.node
        if isinstance(x, mp_Mobspy_Parameter):
            return ParamRefNode(x.get_name())
        if isinstance(x, (int, float)):
            return LiteralNode(x)
        return x

    def __neg__(self) -> MobsPyExpression:
        """Unary minus.

        Negation is valid in three contexts: assignment rules, rate-function
        lambdas (expression mode), and standalone rate building. Each routes
        ``-x`` to ``-1 * x`` in the appropriate representation.
        """
        if asgi_Assign.check_context():
            return asgi_Assign.mul(-1, self)

        from mobspy.dsl.species_operators import _ms_active_ctx  # noqa: PLC0415
        from mobspy.expressions.evaluation import (  # noqa: PLC0415  # circular import
            check_if_non_expression_operated,
        )

        if _ms_active_ctx.get():
            # Lambda replay: wrap as a MobsPyExpression and negate it.
            return -check_if_non_expression_operated(self)  # type: ignore[no-any-return]  # dynamic dispatch

        # Standalone rate building: build -1 * self as a RateExpression.
        from mobspy.expressions.nodes import BinaryOpNode, LiteralNode  # noqa: PLC0415
        from mobspy.expressions.rate_builder import RateExpression  # noqa: PLC0415

        return RateExpression(  # type: ignore[return-value]  # pint Unit subtype
            BinaryOpNode(LiteralNode(-1), "*", self._operand_to_node(self))
        )

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
        self._old_context: set[str] | None = None
        from mobspy.expressions.rate_builder import RateExpression  # noqa: PLC0415

        if isinstance(object_reference, RateExpression):
            raise ReactionError(
                "A rate expression (e.g. from negating or doing arithmetic on a "
                "species, like '-A' or '2 - A') cannot be used as a reactant or "
                "product. Such expressions are only valid as reaction rates, "
                "after the '@' operator."
            )
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
        from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import

        if not self.list_of_reactants:
            return "Zero"
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
        from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import

        item = str(item)
        Species.check_if_valid_characteristic(self, item)
        return self.__getattr__(item)  # type: ignore[no-any-return]  # dynamic dispatch

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
        return _Last_rate_storage.override_get_item(self, item)  # type: ignore[return-value]  # pint Unit subtype

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
        return self.list_of_reactants[0]["characteristics"]  # type: ignore[no-any-return]  # dynamic dispatch

    def __rmul__(self, stoichiometry: Any) -> Self | MobsPyExpression:
        """Multiply by stoichiometry for reactions.

        Args:
            stoichiometry: Stoichiometry value.
        """
        if asgi_Assign.check_context():
            return asgi_Assign.mul(stoichiometry, self)
        if Assignment_Opp_Imp._rate_lambda_active():
            # Inside a rate lambda: ``number * species`` is rate arithmetic.
            return Assignment_Opp_Imp._rate_expression_op(
                stoichiometry, self, "Multiplication"
            )
        if isinstance(stoichiometry, (int, float)):
            if not self.list_of_reactants:
                return self
            self.list_of_reactants[0]["stoichiometry"] = stoichiometry
        else:
            raise ReactionError(
                f"Stoichiometry can only be an int or float - Received {stoichiometry}"
            )
        return self

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
        from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import

        if asgi_Assign.check_context():
            return asgi_Assign.add(self, other)
        if Assignment_Opp_Imp._rate_lambda_active():
            return Assignment_Opp_Imp._rate_expression_op(self, other, "Addition")
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

    def __radd__(self, other: Any) -> Self | RatedProduct | MobsPyExpression:
        if asgi_Assign.check_context():
            return asgi_Assign.add(other, self)
        if Assignment_Opp_Imp._rate_lambda_active():
            return Assignment_Opp_Imp._rate_expression_op(other, self, "Addition")
        return Reacting_Species.__add__(self, other)  # pyright: ignore[reportReturnType]  # pint Quantity subtype

    def __invert__(self) -> Reacting_Species:
        return self.c(NOT_CHAR)

    def __rshift__(self, other: Species | Reacting_Species | RatedProduct) -> Reactions:
        """The ``>>`` operator for defining reactions.

        Accepts Species, Reacting_Species, or RatedProduct (from ``@``).

        Args:
            other: Product side of the reaction being added.
        """
        if isinstance(other, RatedProduct):
            from mobspy.dsl.species import (  # noqa: PLC0415  # circular import
                _create_reaction_from_rated,
            )

            return _create_reaction_from_rated(self.list_of_reactants, other)

        from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import

        p = Reacting_Species(other, set()) if isinstance(other, Species) else other

        return Reactions(self.list_of_reactants, p.list_of_reactants)

    def __call__(  # type: ignore[return]  # noqa: PLR0912
        self,
        quantity: int | float | str | Quantity | mp_Mobspy_Parameter,
    ) -> Self | None:
        """Assign counts to species non-default state.

        Args:
            quantity: Count to be assigned.
        """
        from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import

        if isinstance(quantity, (np_int_, np_float_)):
            quantity = float(quantity)

        _any_ctx = Species.get_meta_specie_named_any_context()
        if len(_any_ctx) > 0:
            for i in _any_ctx:
                self = self.c(i)  # type: ignore[assignment]  # noqa: PLW0642

        if not self.list_of_reactants:
            raise ReactionError("Cannot assign count to empty (Zero) species")
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
            return self  # pyright: ignore[reportReturnType]  # pint Quantity subtype
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
            return self  # pyright: ignore[reportReturnType]  # pint Quantity subtype

    @property
    def ref(self) -> Any:
        """Return a rate expression referencing this species state's concentration."""
        from mobspy.expressions.nodes import SpeciesRefNode  # noqa: PLC0415
        from mobspy.expressions.rate_builder import RateExpression  # noqa: PLC0415

        return RateExpression(SpeciesRefNode(str(self)))

    def __getattr__(self, characteristic: str) -> Any:
        """Implementation of the .dot operation.

        Args:
            characteristic: Characteristic for the query.
        """
        from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import

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
        from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import

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
        from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import

        if self._old_context is not None:
            Species.update_meta_specie_named_any_context(self._old_context)


_methods_Reacting_Species = set(dir(Reacting_Species))  # noqa: N816  # legacy name
