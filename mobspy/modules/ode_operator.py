from mobspy.modules.assignments_implementation import Assign
from mobspy.modules.meta_class import Species
from mobspy.simulation_logging.log_scripts import error as simlog_error


def generate_ODE_reaction_rate(list_of_used_species, expression):
    """Generates a rate function from the ODE expression."""
    expr_string = str(expression._operation)

    # Convert $asg_X to $pos_N based on position in list
    for i, spe in enumerate(list_of_used_species):
        spe_name = spe.get_name()
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

    def __rshift__(self, expression):
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

        reactants >> self.state_variable + reactants [rate_fn]
        # Returns the reaction, even though it doesn't matter
        return self


class DifferentialOperator:
    """Differential operator for ODE syntax: dt[A] >> expression."""

    def __getitem__(self, item):
        if isinstance(item, Species):
            Assign.set_context()  # Turn ON before expression is evaluated
            return ODEBinding(item)
        else:
            simlog_error("MobsPy ODE object must only be applied on a species")


dt = DifferentialOperator()