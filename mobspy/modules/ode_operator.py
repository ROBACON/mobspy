from __future__ import annotations

import inspect
import re
from typing import Any

from mobspy.exceptions import ValidationError
from mobspy.modules.assignments_implementation import Assign
from mobspy.modules.meta_class import Reacting_Species, Species



def generate_ODE_reaction_rate(list_of_used_species: list[Any], expression: Any) -> Any:
    """Generate a rate function from an ODE expression using a closure.

    Replaces species placeholders in the expression string with
    the positional arguments passed at call time.
    """
    expr_template = str(expression._operation)

    # Convert $asg_X to $pos_N based on position in list
    for i, spe in enumerate(list_of_used_species):
        spe_name = str(spe)
        expr_template = expr_template.replace(f"($asg_{spe_name})", f"$_pos_{i}")

    n = len(list_of_used_species)
    param_names = [f"r{i + 1}" for i in range(n)]

    def rate_fn(**kwargs: Any) -> str:
        result = expr_template
        for i, name in enumerate(param_names):
            result = result.replace(f"$_pos_{i}", str(kwargs[name]))
        return result

    # Set proper signature so inspect.signature() returns (r1, r2, ...)
    params = [
        inspect.Parameter(name, inspect.Parameter.POSITIONAL_OR_KEYWORD)
        for name in param_names
    ]
    rate_fn.__signature__ = inspect.Signature(params)  # type: ignore[attr-defined]

    return rate_fn


class ODEBinding:
    """Intermediate object returned by dt[A] that waits for += or -= expression."""

    def __init__(self, state_variable: Species | Reacting_Species) -> None:
        self.state_variable = state_variable

    def _process_ode_expression(
        self,
        expression: Any,
        is_birth: bool = True,
    ) -> ODEBinding:
        """Common logic for processing ODE expressions.

        Args:
            expression: The rate expression
            is_birth: True for += (birth), False for -= (death)

        Returns:
            self for method chaining
        """
        operator = "+=" if is_birth else "-="

        # Validation
        if isinstance(expression, Reacting_Species):  # noqa: SIM102
            if len(expression.list_of_reactants) > 1:
                raise ValidationError((
                        f"ODE expressions must be built within"
                        f" the dt[...] {operator} context.\n"
                        f"Expressions like 'C = A + B' followed"
                        f" by 'dt[X] {operator} C' are not"
                        f" valid.\n"
                        f"Use: dt[X] {operator} A + B"
                    )
                )

        if isinstance(expression, (Species, Reacting_Species)):
            expression = Assign.mul(1, expression)

        Assign.reset_context()

        species_list_operation_order = expression.species_list_operation_order
        rate_fn = generate_ODE_reaction_rate(
            expression.species_list_operation_order, expression
        )

        reactants = None
        for spe in species_list_operation_order:
            if reactants is None:  # noqa: SIM108
                reactants = spe
            else:
                reactants = reactants + spe

        # Create reaction based on type
        if is_birth:
            reactants >> self.state_variable + reactants[rate_fn]
        else:
            reactants + self.state_variable >> reactants[rate_fn]

        return self

    def __iadd__(self, expression: Any) -> ODEBinding:
        """Handles dt[A] += expression (birth/production reactions)."""
        return self._process_ode_expression(expression, is_birth=True)

    def __isub__(self, expression: Any) -> ODEBinding:
        """Handles dt[A] -= expression (death/degradation reactions)."""
        return self._process_ode_expression(expression, is_birth=False)


class DifferentialOperator:
    """Differential operator for ODE syntax: dt[A] += expression."""

    @staticmethod
    def _compile_ode_syntax(code_line: str, line_number: int) -> None:
        """Validate that ODE syntax uses += or -=."""
        if not re.search(r"dt\s*\[.*\]\s*(\+\=|\-\=)", code_line):
            raise ValidationError(
                f"At: {code_line}\n"
                f"Line number: {line_number}\n"
                "ODE syntax requires '+=' or '-=' operator"
                " right after dt[Species] in the same line\n"
                "Use: dt[Species] += expression (for birth)\n"
                "Use: dt[Species] -= expression (for death)"
            )

    def __setitem__(self, key: Any, value: Any) -> None:
        stack_frame = inspect.stack()[1]
        code_line = stack_frame.code_context[0] if stack_frame.code_context else ""

        if re.search(r"dt\s*\[.*\]\s*(\+\=|\-\=)", code_line):
            return  # Valid += or -= syntax, nothing to do

        raise ValidationError((
                "ODE syntax requires '+=' or '-=' operator,"
                " not '=', right after dt[Species] in the"
                " same line\n"
                "Use: dt[Species] += expression (for birth)\n"
                "Use: dt[Species] -= expression (for death)"
            )
        )

    def __getitem__(self, item: Species | Reacting_Species) -> ODEBinding | None:
        if isinstance(item, (Species, Reacting_Species)):
            stack_frame = inspect.stack()[1]
            code_line = stack_frame.code_context[0] if stack_frame.code_context else ""
            line_number = stack_frame.lineno

            self._compile_ode_syntax(code_line, line_number)

            Assign.set_context()  # Turn ON before expression is evaluated
            return ODEBinding(item)
        else:
            raise ValidationError("MobsPy ODE object must only be applied on a species")
            return None


dt = DifferentialOperator()
