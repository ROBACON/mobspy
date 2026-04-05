"""Builder API for constructing rate expressions as AST nodes.

Provides a programmatic alternative to the lambda-replay mechanism.
Instead of writing rate lambdas that are replayed with mock species
objects, users can construct expression trees directly::

    from mobspy.modules.rate_builder import RateExpression, species_ref, param_ref

    # Instead of: lambda r1, r2: k * r1 * r2
    rate = param_ref("k") * species_ref(A) * species_ref(B)

    # Use in reaction:
    A + B >> C @ rate

The builder produces ``ExprNode`` objects that the compiler consumes
directly, without ``sys._getframe`` introspection or lambda replay.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mobspy.modules.expression_nodes import (
    BinaryOpNode,
    ConditionalNode,
    ExprNode,
    FunctionCallNode,
    LiteralNode,
    ParamRefNode,
    SpeciesRefNode,
)

if TYPE_CHECKING:
    from mobspy.modules.species import Species


class RateExpression:
    """Wrapper around ExprNode that supports operator overloading.

    Allows building rate expressions via natural Python syntax::

        rate = RateExpression(2.0) * species_ref(A)
        rate = species_ref(A) ** 2 / param_ref("K")

    The resulting ``.node`` attribute is an ``ExprNode`` ready for
    the compiler.
    """

    __slots__ = ("node",)

    def __init__(self, value: ExprNode | int | float | str) -> None:
        if isinstance(value, ExprNode):
            self.node = value
        elif isinstance(value, (int, float)):
            self.node = LiteralNode(value)
        elif isinstance(value, str):
            self.node = ParamRefNode(value)
        else:
            msg = f"Cannot create RateExpression from {type(value)}"
            raise TypeError(msg)

    def _wrap(self, other: RateExpression | ExprNode | int | float | str) -> ExprNode:
        if isinstance(other, RateExpression):
            return other.node
        if isinstance(other, ExprNode):
            return other
        if isinstance(other, (int, float)):
            return LiteralNode(other)
        if isinstance(other, str):
            return ParamRefNode(other)
        msg = f"Unsupported operand type: {type(other)}"
        raise TypeError(msg)

    def __mul__(
        self, other: RateExpression | ExprNode | int | float | str
    ) -> RateExpression:
        return RateExpression(BinaryOpNode(self.node, "*", self._wrap(other)))

    def __rmul__(
        self, other: RateExpression | ExprNode | int | float | str
    ) -> RateExpression:
        return RateExpression(BinaryOpNode(self._wrap(other), "*", self.node))

    def __truediv__(
        self, other: RateExpression | ExprNode | int | float | str
    ) -> RateExpression:
        return RateExpression(BinaryOpNode(self.node, "/", self._wrap(other)))

    def __rtruediv__(
        self, other: RateExpression | ExprNode | int | float | str
    ) -> RateExpression:
        return RateExpression(BinaryOpNode(self._wrap(other), "/", self.node))

    def __add__(
        self, other: RateExpression | ExprNode | int | float | str
    ) -> RateExpression:
        return RateExpression(BinaryOpNode(self.node, "+", self._wrap(other)))

    def __radd__(
        self, other: RateExpression | ExprNode | int | float | str
    ) -> RateExpression:
        return RateExpression(BinaryOpNode(self._wrap(other), "+", self.node))

    def __sub__(
        self, other: RateExpression | ExprNode | int | float | str
    ) -> RateExpression:
        return RateExpression(BinaryOpNode(self.node, "-", self._wrap(other)))

    def __rsub__(
        self, other: RateExpression | ExprNode | int | float | str
    ) -> RateExpression:
        return RateExpression(BinaryOpNode(self._wrap(other), "-", self.node))

    def __pow__(
        self, other: RateExpression | ExprNode | int | float | str
    ) -> RateExpression:
        return RateExpression(BinaryOpNode(self.node, "^", self._wrap(other)))

    def __rpow__(
        self, other: RateExpression | ExprNode | int | float | str
    ) -> RateExpression:
        return RateExpression(BinaryOpNode(self._wrap(other), "^", self.node))

    def __neg__(self) -> RateExpression:
        return RateExpression(BinaryOpNode(LiteralNode(-1), "*", self.node))

    def where(
        self,
        condition: str,
        if_false: RateExpression | int | float,
    ) -> RateExpression:
        """Conditional expression: ``self`` if condition else ``if_false``.

        Equivalent to ``lambda r: self if r.is_a(X) else if_false``
        in the lambda syntax.

        Args:
            condition: SBML-compatible boolean expression string.
            if_false: Value when condition is false.

        Examples:
            >>> from mobspy.modules.rate_builder import literal
            >>> rate = literal(0.5).where("A_dot_alive > 0", 1.0)
            >>> "piecewise" in str(rate)
            True
        """
        return RateExpression(
            ConditionalNode(
                condition=condition,
                if_true=self.node,
                if_false=self._wrap(if_false),
            )
        )

    def __str__(self) -> str:
        return self.node.render()

    def __repr__(self) -> str:
        return f"RateExpression({self.node.render()})"


# ---------------------------------------------------------------
# Factory functions
# ---------------------------------------------------------------


def species_ref(
    species: Species | str,
    mode: str = "default",
) -> RateExpression:
    """Create a species reference in a rate expression.

    Args:
        species: Species object or SBML species name string.
        mode: ``"default"``, ``"count"``, or ``"concentration"``.

    Examples:
        >>> from mobspy.modules.rate_builder import species_ref
        >>> ref = species_ref("A")
        >>> str(ref)
        'A'
    """
    name = species if isinstance(species, str) else species.get_name()
    return RateExpression(SpeciesRefNode(name, mode=mode))


def param_ref(name: str) -> RateExpression:
    """Create a parameter reference in a rate expression.

    Args:
        name: Parameter name (must match a ModelParameters name).

    Examples:
        >>> from mobspy.modules.rate_builder import param_ref
        >>> ref = param_ref("k1")
        >>> str(ref)
        'k1'
    """
    return RateExpression(ParamRefNode(name))


def literal(value: int | float) -> RateExpression:
    """Create a numeric literal in a rate expression.

    Examples:
        >>> from mobspy.modules.rate_builder import literal
        >>> str(literal(2.5))
        '2.5'
    """
    return RateExpression(LiteralNode(value))


def where(
    condition: str,
    if_true: RateExpression | int | float,
    if_false: RateExpression | int | float,
) -> RateExpression:
    """Piecewise conditional expression.

    Equivalent to ``if_true if condition else if_false``.

    Args:
        condition: SBML-compatible boolean expression string.
        if_true: Value when condition is true.
        if_false: Value when condition is false.

    Examples:
        >>> from mobspy.modules.rate_builder import where, literal
        >>> rate = where("A_dot_alive > 0", literal(0.5), literal(1.0))
        >>> "piecewise" in str(rate)
        True
    """
    true_node = (
        if_true.node if isinstance(if_true, RateExpression) else LiteralNode(if_true)
    )
    false_node = (
        if_false.node if isinstance(if_false, RateExpression) else LiteralNode(if_false)
    )
    return RateExpression(ConditionalNode(condition, true_node, false_node))


def rate_log(expr: RateExpression) -> RateExpression:
    """Natural logarithm in a rate expression."""
    return RateExpression(FunctionCallNode("log", expr.node))


def rate_exp(expr: RateExpression) -> RateExpression:
    """Exponential function in a rate expression."""
    return RateExpression(FunctionCallNode("exp", expr.node))


def rate_sqrt(expr: RateExpression) -> RateExpression:
    """Square root in a rate expression."""
    return RateExpression(FunctionCallNode("sqrt", expr.node))


def hill(
    species: Species | str,
    vmax: int | float | str,
    km: int | float | str,
    n: int | float = 1,
) -> RateExpression:
    """Hill function: vmax * S^n / (km^n + S^n).

    Common rate expression in biochemical kinetics.

    Args:
        species: Substrate species.
        vmax: Maximum rate.
        km: Half-saturation constant.
        n: Hill coefficient.
    """
    s = species_ref(species)
    v = literal(vmax) if isinstance(vmax, (int, float)) else param_ref(vmax)
    s_n = s**n
    k_n = literal(km) ** n if isinstance(km, (int, float)) else param_ref(km) ** n
    return v * s_n / (k_n + s_n)
