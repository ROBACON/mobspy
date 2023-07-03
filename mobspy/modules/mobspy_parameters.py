import inspect
import mobspy.simulation_logging.log_scripts as simlog
from mobspy.modules.mobspy_expressions import *
from pint import Quantity, UnitRegistry


class Mobspy_Parameter(ExpressionDefiner, QuantityConverter):

    # convert_received_unit

    parameter_stack = {}

    def __init__(self, name, value):
        self._generate_necessary_attributes()

        temp_set = set()
        temp_set.add(self)
        self.name = name
        self.original_value = value
        self.parameter_stack[name] = self

        self._ms_active = True

        self._operation = str(name)
        self._parameter_set.add(self)

        if isinstance(value, Quantity):
            # We convert into MobsPy units already during the definition of a parameter
            self.value = self.convert_received_unit(value)

            self._unit_count_op = value
            self._unit_conc_op = value
            self._unit_operation = value
            self._has_units = 'T'
        else:
            self.value = value

            self._unit_count_op = 1
            self._unit_conc_op = 1
            self._has_units = False

    def rename(self, new_name):
        del self.parameter_stack[self.name]
        self.parameter_stack[new_name] = self
        self.name = new_name

    def __str__(self):
        return str(self._operation)


def ModelParameters(*args):

    code_line = inspect.stack()[1].code_context[0][:-1]
    separated_line = code_line.split('=')[-2].replace(" ", "")
    parameter_variable_names = separated_line.split(',')

    if len(args) != len(parameter_variable_names):
        simlog.error('The number of parameters provided does not match the number of variables', stack_index=1)

    if len(parameter_variable_names) > 1:
        parameters_to_return = [Mobspy_Parameter(p, v) for p, v in zip(parameter_variable_names, args)]
    else:
        parameters_to_return = Mobspy_Parameter(parameter_variable_names[0], args[0])

    return parameters_to_return


if __name__ == '__main__':
    u = UnitRegistry()
    a, b, c = ModelParameters(1, [3, 4, 5], 2)
    r1 = (a + b + c)/5
    print(r1._operation)
    # print(type(r1._parameter_set))
    # print(Mobspy_Parameter.parameter_stack)

