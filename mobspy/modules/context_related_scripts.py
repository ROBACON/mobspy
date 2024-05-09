from mobspy.modules.mobspy_expressions import ExpressionDefiner as me_ExpressionDefiner, u
import inspect

def expression_compilation_initiation():
    """
        First find all ExpressionDefiner objects in the stack.
        Then, sets _ms_active to True to change the behavior of quantities and units in expressions
    """
    u._ms_active = True

    expressions_in_stack = []

    for i in range(len(inspect.stack())):
        local_names = inspect.stack()[i][0].f_locals
        global_names = inspect.stack()[i][0].f_globals
        for key, item in global_names.items():
            if isinstance(item, me_ExpressionDefiner):
                expressions_in_stack.append(item)
        for key, item in local_names.items():
            if isinstance(item, me_ExpressionDefiner):
                expressions_in_stack.append(item)

    for expression in expressions_in_stack:
        expression._ms_active = True

    return expressions_in_stack


def expression_compilation_finish(expressions):
    """
        Sets _ms_active to false for all expressions found in the stack by expression_compilation_initiation()

        :param expressions: all expressions found in stack at the moment of reaction compilation
    """
    u._ms_active = False

    for expression in expressions:
        expression._ms_active = False


class Unit_Context_Setter:
    def __enter__(self):
        self.expressions_in_stack = expression_compilation_initiation()

    def __exit__(self, exc_type, exc_value, traceback):
        expression_compilation_finish(self.expressions_in_stack)
