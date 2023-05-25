import inspect
import mobspy.simulation_logging.log_scripts as simlog


class Parameter_Operations:

    def __init__(self, operation, parameter_set):
        self.operation = operation
        self.parameter_set = parameter_set

    def process_other_operator(self, other):
        if isinstance(other, Parameter_Operations):
            self.parameter_set = self.parameter_set.union(other.parameter_set)

    def __add__(self, other):
        self.process_other_operator(other)
        return Parameter_Operations('(' + self.operation + ' + ' + str(other) + ')', self.parameter_set)

    def __radd__(self, other):
        self.process_other_operator(other)
        return Parameter_Operations('(' + self.operation + ' + ' + str(other) + ')', self.parameter_set)

    def __sub__(self, other):
        self.process_other_operator(other)
        return Parameter_Operations('(' + self.operation + ' - ' + str(other) + ')', self.parameter_set)

    def __rsub__(self, other):
        self.process_other_operator(other)
        return Parameter_Operations('(' + str(other) + ' - ' + self.operation + ')', self.parameter_set)

    def __mul__(self, other):
        self.process_other_operator(other)
        return Parameter_Operations('(' + self.operation + '*' + str(other) + ')', self.parameter_set)

    def __rmul__(self, other):
        self.process_other_operator(other)
        return Parameter_Operations('(' + str(other) + '*' + self.operation + ')', self.parameter_set)

    def __truediv__(self, other):
        self.process_other_operator(other)
        return Parameter_Operations('(' + self.operation + '/' + str(other) + ')', self.parameter_set)

    def __rtruediv__(self, other):
        self.process_other_operator(other)
        return Parameter_Operations('(' + str(other) + '/' + self.operation + ')', self.parameter_set)

    def __str__(self):
        return str(self.operation)


class Mobspy_Parameter(Parameter_Operations):

    parameter_stack = {}

    def __init__(self, name, value):
        super().__init__(name, self)
        self.name = name
        self.value = value
        self.parameter_stack[name] = self

    def rename(self, new_name):
        del self.parameter_stack[self.name]
        self.parameter_stack[new_name] = self


def MSParameters(*args):

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
    a, b, c = MSParameters(1, [3, 4, 5], 2)
    r1 = a + b
    print(r1.parameter_set)
    # r1 = (a + b + c)/5
    # print(r1.operation)
    # print(Mobspy_Parameter.parameter_stack)

