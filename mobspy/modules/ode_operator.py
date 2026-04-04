"""Convert ODE-style differential equations into MobsPy reaction rate closures."""

from __future__ import annotations

import inspect
from typing import Any

from mobspy.constants import ASSIGNMENT_PREFIX, POSITION_PREFIX
from mobspy.exceptions import ValidationError
from mobspy.modules.assignments_implementation import Assign
from mobspy.modules.reactions import Reacting_Species
from mobspy.modules.species import Species


def generate_ODE_reaction_rate(list_of_used_species: list[Any], expression: Any) -> Any:  # noqa: N802
    """Generate a rate function from an ODE expression using a closure.

    Replaces species placeholders in the expression string with
    the positional arguments passed at call time.
    """
    expr_template = str(expression._operation)

    # Convert $asg_X to $pos_N based on position in list
    for i, spe in enumerate(list_of_used_species):
        spe_name = str(spe)
        old = f"({ASSIGNMENT_PREFIX}{spe_name})"
        new = f"{POSITION_PREFIX}{i}"
        expr_template = expr_template.replace(old, new)

    n = len(list_of_used_species)
    param_names = [f"r{i + 1}" for i in range(n)]

    def rate_fn(**kwargs: Any) -> str:
        """Evaluate the ODE rate expression with given species values."""
        result = expr_template
        for i, name in enumerate(param_names):
            result = result.replace(f"{POSITION_PREFIX}{i}", str(kwargs[name]))
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
        if (
            isinstance(expression, Reacting_Species)
            and len(expression.list_of_reactants) > 1
        ):
            raise ValidationError(
                "ODE expressions must be built within"
                f" the dt[...] {operator} context.\n"
                "Expressions like 'C = A + B' followed"
                f" by 'dt[X] {operator} C' are not"
                " valid.\n"
                f"Use: dt[X] {operator} A + B"
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
            reactants = spe if reactants is None else reactants + spe

        # Create reaction based on type (return value unused; side effect registers it)
        if is_birth:
            _ = reactants >> (self.state_variable + reactants) @ rate_fn
        else:
            _ = (reactants + self.state_variable) >> reactants @ rate_fn

        return self

    def __iadd__(self, expression: Any) -> ODEBinding:
        """Handles dt[A] += expression (birth/production reactions)."""
        return self._process_ode_expression(expression, is_birth=True)

    def __isub__(self, expression: Any) -> ODEBinding:
        """Handles dt[A] -= expression (death/degradation reactions)."""
        return self._process_ode_expression(expression, is_birth=False)


class DifferentialOperator:
    """Differential operator for ODE syntax: dt[A] += expression."""

    def __setitem__(self, key: Any, value: Any) -> None:
        # With += / -=, Python calls __getitem__ -> __iadd__/__isub__ -> __setitem__.
        # ODEBinding.__iadd__/__isub__ return self, so value will be an ODEBinding.
        # Plain `dt[X] = expr` passes the raw expression, which is invalid.
        if isinstance(value, ODEBinding):
            return

        raise ValidationError(
            "ODE syntax requires '+=' or '-=' operator,"
            " not '=', right after dt[Species] in the"
            " same line\n"
            "Use: dt[Species] += expression (for birth)\n"
            "Use: dt[Species] -= expression (for death)"
        )

    def __getitem__(self, item: Species | Reacting_Species) -> ODEBinding | None:
        if isinstance(item, (Species, Reacting_Species)):
            Assign.set_context()  # Turn ON before expression is evaluated
            return ODEBinding(item)
        raise ValidationError("MobsPy ODE object must only be applied on a species")


dt = DifferentialOperator()
