from inspect import stack as inspect_stack
import mobspy.simulation_logging.log_scripts as simlog
from mobspy.modules.mobspy_expressions import ExpressionDefiner as me_ExpressionDefiner, \
    QuantityConverter as me_QuantityConverter
from pint import Quantity, UnitRegistry


class Mobspy_Parameter(me_ExpressionDefiner, me_QuantityConverter):
    """
        This is the constructor that is called by ModelParameters to create a model parameter
        (not a simulation parameter). The user is not supposed to create parameters using this object.
    """
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

        def unit_process(value):
            # We convert into MobsPy units already during the definition of a parameter
            self.value = self.convert_received_unit(value).magnitude

            self.original_magnitude = value.magnitude
            self.conversion_factor = self.value / self.original_magnitude
            self.original_unit = value.units

            # For future developers, the T is there to avoid potential bugs with the . query
            self._unit_count_op = value
            self._unit_conc_op = value
            self._unit_operation = value
            self._has_units = 'T'

            return self.value, self.original_unit

        if isinstance(value, Quantity):
            unit_process(value)
        elif type(value) == list or type(value) == tuple:
            new_list = []
            first_unit = None
            for i, val in enumerate(value):
                if isinstance(val, Quantity) and i == 0:
                    new_value, first_unit = unit_process(val)
                elif isinstance(val, Quantity) and i > 0 and first_unit is None:
                    simlog.error("MobsPy parameters must all be the same unit", stack_index=1)
                elif isinstance(val, Quantity) and i > 0 and first_unit is not None:
                    new_value, unit = unit_process(val)
                    if unit != first_unit:
                        simlog.error("MobsPy parameters must all be the same unit", stack_index=1)
                else:
                    new_value = val

                new_list.append(new_value)

            self.value = new_list
        else:
            self.value = value

            self._unit_count_op = 1
            self._unit_conc_op = 1
            self.conversion_factor = 1
            self._has_units = False

    def convert_to_original_unit(self):
        """
            Converts the parameter from the MobsPy standard unit to the original unit of the parameter
        """
        if self.has_units():
            self.set_value(self.value/self.conversion_factor*self.original_unit)

    def rename(self, new_name):
        """
            Renames a parameter. Checks to see if name is available beforehand. It uses the parameter stack to do so
        """
        if new_name in self.parameter_stack:
            simlog.warning(" MobsPy uses a parameter dictionary with parameter names as keys and the respective object"
                           " as value to keep track of created parameters. As, there is a parameter with this name"
                           " already in the stack. The old will be deleted and replaced by this one.")

        del self.parameter_stack[self.name]
        self.parameter_stack[new_name] = self
        self.name = new_name

    def set_value(self, new_value):
        """
            Sets value of parameter
        """
        self.value = new_value
        return self

    def has_units(self):
        """
            Check if is a unit based parameter or not
        """
        if self._has_units == 'T':
            return True
        else:
            return False

    def __str__(self):
        return str(self._operation)


def ModelParameters(*args):
    """
        Creates ModelParameters. Like meta-species, it uses the variable names as parameter names
    """
    code_line = inspect_stack()[1].code_context[0][:-1]
    separated_line = code_line.split('=')[-2].replace(" ", "")
    parameter_variable_names = separated_line.split(',')

    if len(args) != len(parameter_variable_names):
        simlog.error('You must provide an initial value for every parameter variable declared', stack_index=2)

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

