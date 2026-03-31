"""Mathematical function wrappers for use inside MobsPy rate expressions."""

from __future__ import annotations

from typing import Any

from mobspy.constants import NULL_SPECIES
from mobspy.exceptions import ValidationError
from mobspy.modules.assignments_implementation import Assign
from mobspy.modules.meta_class import Reacting_Species, Species
from mobspy.modules.mobspy_expressions import (
    ExprNode,
    FunctionCallNode,
    MobsPyExpression,
    _to_expr_node,
)


class MathFunctionWrapper:
    """Wrapper for mathematical functions that work with MobsPy expressions

    Note for future developers. To add units to these functions,
    use the unit compilation to check if dealing with concentrations or counts

    The units must be compiled inside a function as the input is non-dimensional.
    Therefore, the unit compilation must be performed here at this step

    Currently, for the ODE application, I don't need units,
    so I need to leave this for future needs
    """

    def __init__(self, name: str) -> None:
        self.name = name  # COPASI function name: 'exp', 'sin', 'cos', etc.

    def _create_expression(
        self, expression: MobsPyExpression, new_operation: ExprNode | str
    ) -> MobsPyExpression:
        """Create new MobsPyExpression with this function applied."""
        return MobsPyExpression(
            species_string=NULL_SPECIES,
            species_object=None,
            operation=new_operation,
            unit_count_op=1,
            unit_conc_op=1,
            dimension=expression._dimension,
            expression_variables=set(expression._expression_variables),
            parameter_set=set(expression._parameter_set),
            count_in_model=expression._count_in_model,
            concentration_in_model=expression._concentration_in_model,
            count_in_expression=expression._count_in_expression,
            concentration_in_expression=expression._concentration_in_expression,
            has_units=expression._has_units,
            species_list_operation_order=list(expression.species_list_operation_order),
        )

    def __call__(self, expression: Any) -> MobsPyExpression | None:
        if not Assign.check_context():
            raise ValidationError(
                f"ms_{self.name}() can only be called inside a rate expression "
                "(reaction rate lambda or ODE assignment)"
            )

        # MobsPy Expressions
        if isinstance(expression, MobsPyExpression):
            if expression._has_units:
                raise ValidationError(
                    f"ms_{self.name}() does not support unit-bearing expressions. "
                    "Extract the numeric value before applying the function."
                )

            new_operation = FunctionCallNode(
                self.name, _to_expr_node(expression._operation)
            )
            return self._create_expression(expression, new_operation)

        # Species passed
        if (
            isinstance(expression, (Species, Reacting_Species))
        ) and Assign.check_context():
            if (
                isinstance(expression, Reacting_Species)
                and len(expression.list_of_reactants) > 1
            ):
                raise ValidationError(
                    "Reacting species with multiple"
                    " reactants should not be applied"
                    " to a function"
                )

            expression = Assign.mul(1, expression)
            new_operation = FunctionCallNode(
                self.name, _to_expr_node(expression._operation)
            )
            return self._create_expression(expression, new_operation)

        raise ValidationError(
            f"ms_{self.name}() received an unsupported "
            f"argument type: {type(expression).__name__}."
            " Expected a species, MobsPyExpression,"
            " or numeric value."
        )
        return None


# Create all the COPASI-compatible math functions
# Please add new functions with ms to avoid conflict with Python's existing modules
ms_exp = MathFunctionWrapper("exp")
ms_logn = MathFunctionWrapper("log")
ms_log10 = MathFunctionWrapper("log10")
ms_sin = MathFunctionWrapper("sin")
ms_cos = MathFunctionWrapper("cos")
ms_tan = MathFunctionWrapper("tan")
ms_asin = MathFunctionWrapper("asin")
ms_acos = MathFunctionWrapper("acos")
ms_atan = MathFunctionWrapper("atan")
ms_sinh = MathFunctionWrapper("sinh")
ms_cosh = MathFunctionWrapper("cosh")
ms_tanh = MathFunctionWrapper("tanh")
ms_floor = MathFunctionWrapper("floor")
ms_ceil = MathFunctionWrapper("ceil")
ms_abs = MathFunctionWrapper("abs")
ms_sqrt = MathFunctionWrapper("sqrt")
