from mobspy.modules.assignments_implementation import Assign
from mobspy.modules.meta_class import Species, Reacting_Species
from mobspy.simulation_logging.log_scripts import error as simlog_error
import re
from inspect import stack as inspect_stack


def generate_ODE_reaction_rate(list_of_used_species, expression):
    """Generates a rate function from the ODE expression."""
    expr_string = str(expression._operation)

    # Convert $asg_X to $pos_N based on position in list
    for i, spe in enumerate(list_of_used_species):
        spe_name = str(spe)
        expr_string = expr_string.replace(f"($asg_{spe_name})", f"$_pos_{i}")

    # Create the parameter argument r1, r2, r3
    n = len(list_of_used_species)
    param_names = [f"r{i + 1}" for i in range(n)]
    param_str = ", ".join(param_names)

    # Build replacement logic
    func_code = f"def rate_fn({param_str}):\n\t"
    func_code += f"result = {repr(expr_string)}\n\t"

    replace_lines = ""
    for i in range(n):
        replace_lines += f'result = result.replace("$_pos_{i}", str({param_names[i]}))\n\t'
    func_code += replace_lines
    func_code += "return result"

    local_vars = {}
    exec(func_code, {}, local_vars)
    return local_vars["rate_fn"]


class ODEBinding:
    """Intermediate object returned by dt[A] that waits for >> expression."""

    def __init__(self, state_variable):
        self.state_variable = state_variable

    def _process_ode_expression(self, expression, is_birth=True):
        """
        Common logic for processing ODE expressions.

        Args:
            expression: The rate expression
            is_birth: True for += (birth), False for -= (death)

        Returns:
            self for method chaining
        """
        operator = "+=" if is_birth else "-="

        # Validation
        if isinstance(expression, Reacting_Species):
            if len(expression.list_of_reactants) > 1:
                simlog_error(
                    message=f"ODE expressions must be built within the dt[...] {operator} context.\n"
                            f"Expressions like 'C = A + B' followed by 'dt[X] {operator} C' are not valid.\n"
                            f"Use: dt[X] {operator} A + B",
                    full_exception_log=True
                )

        if isinstance(expression, Species) or isinstance(expression, Reacting_Species):
            expression = Assign.mul(1, expression)

        Assign.reset_context()

        species_list_operation_order = expression.species_list_operation_order
        rate_fn = generate_ODE_reaction_rate(
            expression.species_list_operation_order, expression
        )

        reactants = None
        for spe in species_list_operation_order:
            if reactants is None:
                reactants = spe
            else:
                reactants = reactants + spe

        # Create reaction based on type
        if is_birth:
            # Birth reaction: reactants >> state_variable + reactants
            reactants >> self.state_variable + reactants[rate_fn]
        else:
            # Death reaction: reactants + state_variable >> reactants
            reactants + self.state_variable >> reactants[rate_fn]

        return self

    def __iadd__(self, expression):
        """Handles dt[A] += expression (birth/production reactions)."""
        return self._process_ode_expression(expression, is_birth=True)

    def __isub__(self, expression):
        """Handles dt[A] -= expression (death/degradation reactions)."""
        return self._process_ode_expression(expression, is_birth=False)



class DifferentialOperator:
    """Differential operator for ODE syntax: dt[A] >> expression."""


    @staticmethod
    def _compile_ode_syntax(code_line, line_number):
        """Validate that ODE syntax uses += or -="""
        # Check that dt[...] is followed by += or -=
        if not re.search(r'dt\s*\[.*\]\s*(\+\=|\-\=)', code_line):
            simlog_error(
                f"At: {code_line}\n"
                f"Line number: {line_number}\n"
                "ODE syntax requires '+=' or '-=' operator right after dt[Species] in the same line\n"
                "Use: dt[Species] += expression (for birth)\n"
                "Use: dt[Species] -= expression (for death)"
            )

    def __setitem__(self, key, value):
        stack_frame = inspect_stack()[1]
        code_line = stack_frame.code_context[0] if stack_frame.code_context else ""

        if re.search(r'dt\s*\[.*\]\s*(\+\=|\-\=)', code_line):
            return  # Valid += or -= syntax, nothing to do

        simlog_error(
            message="ODE syntax requires '+=' or '-=' operator, not '=', right after dt[Species] in the same line\n"
                    "Use: dt[Species] += expression (for birth)\n"
                    "Use: dt[Species] -= expression (for death)",
            full_exception_log=True
        )

    def __getitem__(self, item):
        if isinstance(item, Species) or isinstance(item, Reacting_Species):
            stack_frame = inspect_stack()[1]
            code_line = stack_frame.code_context[0] if stack_frame.code_context else ""
            line_number = stack_frame.lineno

            self._compile_ode_syntax(code_line, line_number)

            Assign.set_context()  # Turn ON before expression is evaluated
            return ODEBinding(item)
        else:
            simlog_error("MobsPy ODE object must only be applied on a species")


dt = DifferentialOperator()