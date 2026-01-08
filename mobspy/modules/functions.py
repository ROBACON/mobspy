from mobspy.modules.meta_class import Species, Reacting_Species
from mobspy.modules.mobspy_expressions import MobsPyExpression, OverrideQuantity
from mobspy.modules.assignments_implementation import Assign
from mobspy.simulation_logging.log_scripts import error as simlog_error
import math

class MathFunctionWrapper:
    """Wrapper for mathematical functions that work with MobsPy expressions

        Note for future developers. To add units to these functions,
        use the unit compilation to check if dealing with concentrations or counts

        The units must be compiled inside a function as the input is non-dimensional.
        Therefore, the unit compilation must be performed here at this step

        Currently, for the ODE application, I don't need units,
        so I need to leave this for future needs
    """

    def __init__(self, name):
        self.name = name  # COPASI function name: 'exp', 'sin', 'cos', etc.

    def _create_expression(self, expression, new_operation):
        """Create new MobsPyExpression with this function applied."""
        return MobsPyExpression(
            species_string="$Null",
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
            species_list_operation_order=list(expression.species_list_operation_order)
        )

    def __call__(self, expression):

        if not Assign.check_context():
            simlog_error(
                "The expression functions must only be called "
            )

        # MobsPy Expressions
        if isinstance(expression, MobsPyExpression):

            if expression._has_units == 'T':
                simlog_error('At this current version, MobsPy functions do not support ')

            new_operation = f"{self.name}({expression._operation})"
            return self._create_expression(expression, new_operation)

        # Species passed
        elif ((isinstance(expression, Species) or isinstance(expression, Reacting_Species))
              and Assign.check_context()):

            if isinstance(expression, Reacting_Species):
                if len(expression.list_of_reactants) > 1:
                    simlog_error(
                        message="Reacting species with multiple reactants should not be applied to a function",
                        full_exception_log=True
                    )

            expression = Assign.mul(1, expression)
            new_operation = f"{self.name}({expression._operation})"
            return self._create_expression(expression, new_operation)

        else:
            simlog_error(
                message="MobsPy functions were called on a non-valid context",
                full_exception_log=True
            )



# Create all the COPASI-compatible math functions
# Please add new functions with ms to avoid conflict with Python's existing modules
ms_exp = MathFunctionWrapper('exp')
ms_logn = MathFunctionWrapper('log')
ms_log10 = MathFunctionWrapper('log10')
ms_sin = MathFunctionWrapper('sin')
ms_cos = MathFunctionWrapper('cos')
ms_tan = MathFunctionWrapper('tan')
ms_asin = MathFunctionWrapper('asin')
ms_acos = MathFunctionWrapper('acos')
ms_atan = MathFunctionWrapper('atan')
ms_sinh = MathFunctionWrapper('sinh')
ms_cosh = MathFunctionWrapper('cosh')
ms_tanh = MathFunctionWrapper('tanh')
ms_floor = MathFunctionWrapper('floor')
ms_ceil = MathFunctionWrapper('ceil')
ms_abs = MathFunctionWrapper('abs')
ms_sqrt = MathFunctionWrapper('sqrt')