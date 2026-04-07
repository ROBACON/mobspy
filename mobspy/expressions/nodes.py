"""AST node types for representing and rendering MobsPy rate expressions."""

from __future__ import annotations

from typing import Any

from mobspy.constants import ASSIGNMENT_PREFIX, CONCENTRATION_PREFIX, COUNT_PREFIX
from mobspy.types import RenderContext

__all__ = [
    "BinaryOpNode",
    "ConditionalNode",
    "ExprNode",
    "FunctionCallNode",
    "LiteralNode",
    "ParamRefNode",
    "SpeciesRefNode",
    "_render_resolved",
    "_to_expr_node",
]


class ExprNode:
    """Base class for expression AST nodes."""

    def render(self) -> str:
        """Render this AST node as a string expression for SBML."""
        raise NotImplementedError

    def __str__(self) -> str:
        return self.render()

    def __repr__(self) -> str:
        return self.render()

    def walk_species(self) -> list[SpeciesRefNode]:
        """Return all SpeciesRefNode leaves in this subtree."""
        return []


class LiteralNode(ExprNode):
    """A numeric or string literal."""

    __slots__ = ("value",)

    def __init__(self, value: int | float | str) -> None:
        self.value = value

    def render(self) -> str:
        return str(self.value)

    def walk_species(self) -> list[SpeciesRefNode]:
        return []


class SpeciesRefNode(ExprNode):
    """Reference to a species variable in an expression.

    mode controls the string prefix:
      'default'       -> species_string  (bare name)
      'count'         -> $count$species_string
      'concentration' -> $concentration$species_string
      'assignment'    -> ($asg_species_string)
    """

    __slots__ = ("mode", "species_string")

    def __init__(self, species_string: str, mode: str = "default") -> None:
        self.species_string = species_string
        self.mode = mode

    def render(self) -> str:
        if self.mode == "count":
            return COUNT_PREFIX + self.species_string
        if self.mode == "concentration":
            return CONCENTRATION_PREFIX + self.species_string
        if self.mode == "assignment":
            return "(" + ASSIGNMENT_PREFIX + self.species_string + ")"
        return self.species_string

    def walk_species(self) -> list[SpeciesRefNode]:
        return [self]


class ParamRefNode(ExprNode):
    """Reference to a named parameter."""

    __slots__ = ("name",)

    def __init__(self, name: str) -> None:
        self.name = name

    def render(self) -> str:
        return str(self.name)

    def walk_species(self) -> list[SpeciesRefNode]:
        return []


class BinaryOpNode(ExprNode):
    """A binary operation: (left op right)."""

    __slots__ = ("left", "op", "right")

    def __init__(self, left: ExprNode, op: str, right: ExprNode) -> None:
        self.left = left
        self.op = op
        self.right = right

    def render(self) -> str:
        return "(" + self.left.render() + self.op + self.right.render() + ")"

    def walk_species(self) -> list[SpeciesRefNode]:
        return self.left.walk_species() + self.right.walk_species()


class FunctionCallNode(ExprNode):
    """A function call: name(arg)."""

    __slots__ = ("arg", "name")

    def __init__(self, name: str, arg: ExprNode) -> None:
        self.name = name
        self.arg = arg

    def render(self) -> str:
        return self.name + "(" + self.arg.render() + ")"

    def walk_species(self) -> list[SpeciesRefNode]:
        return self.arg.walk_species()


class ConditionalNode(ExprNode):
    """Piecewise conditional: if condition then if_true else if_false.

    Renders as ``piecewise(if_true, condition, if_false)`` for SBML.
    """

    __slots__ = ("condition", "if_false", "if_true")

    def __init__(self, condition: str, if_true: ExprNode, if_false: ExprNode) -> None:
        self.condition = condition
        self.if_true = if_true
        self.if_false = if_false

    def render(self) -> str:
        return (
            f"piecewise({self.if_true.render()}, "
            f"{self.condition}, "
            f"{self.if_false.render()})"
        )

    def walk_species(self) -> list[SpeciesRefNode]:
        return self.if_true.walk_species() + self.if_false.walk_species()


def _to_expr_node(value: ExprNode | int | float | str) -> ExprNode:
    """Wrap a raw value into an ExprNode if it isn't one already."""
    if isinstance(value, ExprNode):
        return value
    return LiteralNode(value)


def _render_resolved(
    node: ExprNode | int | float,
    expression_vars: set[Any],
    render_ctx: RenderContext,
) -> str:
    """Render an AST with species references resolved for count/concentration context.

    Handles the mode on each SpeciesRefNode:
      - 'count' in count_in_model -> bare name;
        in concentration_in_model -> (name*c1)
      - 'concentration' in count_in_model -> (name/c1);
        in concentration_in_model -> bare
      - 'default' resolved based on count_in_expression/concentration_in_expression
      - 'assignment' -> rendered as-is (for ODE references)
    """
    var_names = {v.species_string for v in expression_vars}
    count_in_model = render_ctx.count_in_model
    concentration_in_model = render_ctx.concentration_in_model
    count_in_expression = render_ctx.count_in_expression
    concentration_in_expression = render_ctx.concentration_in_expression

    def _resolve(n: ExprNode) -> str:  # noqa: PLR0911  # many returns
        """Resolve an expression node to its SBML string.

        Applies count/concentration conversion as needed.
        """
        if isinstance(n, SpeciesRefNode) and n.species_string in var_names:
            name = n.species_string
            if n.mode == "assignment":
                return n.render()
            if n.mode == "count":
                if count_in_model:
                    return name
                if concentration_in_model:
                    return "(" + name + "*c1)"
                return name
            if n.mode == "concentration":
                if count_in_model:
                    return "(" + name + "/c1)"
                if concentration_in_model:
                    return name
                return name
            # default mode - resolved by expression context
            # Count takes priority when both flags are set
            if (
                count_in_model
                and concentration_in_expression
                and not count_in_expression
            ):
                return "(" + name + "/c1)"
            if (
                concentration_in_model
                and count_in_expression
                and not concentration_in_expression
            ):
                return "(" + name + "*c1)"
            return name

        if isinstance(n, BinaryOpNode):
            return "(" + _resolve(n.left) + n.op + _resolve(n.right) + ")"
        if isinstance(n, FunctionCallNode):
            return n.name + "(" + _resolve(n.arg) + ")"
        return n.render()

    return _resolve(_to_expr_node(node))
