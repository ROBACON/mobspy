from pint import Quantity, UnitRegistry
import mobspy.simulation_logging.log_scripts as simlog
import re
from scipy.constants import N_A
from copy import deepcopy
from numpy import int_ as np_int_, float_ as np_float_, add as np_add, subtract as np_subtract, \
    multiply as np_multiply, divide as np_divide
import time

# @TODO change _ms_active to context based name
class Bool_Override:
    """
            Just a base class for implementing the . operation in the rate function arguments
            through boolean overriding. It is responsible for returning true when the reactant
            has the specified characteristic when using the dot notation

            :param _stocked_characteristics: (str) stocks the characteristics of the queries performed by the user
            :param species_string: (str) string value from an individual species in MobsPY format
    """

    def __bool__(self):
        """
            The implementation of the .dot operation for rate function arguments
            Returns true when the argument possesses the characteristics
            Bool is called after the .dot operations

            Parameters:
                self

            Returns:
                True if the object contains all the characteristics queried
                False otherwise
        """
        if self.species_string == '$Null':
            return False

        species_string_split = self.species_string.split('_dot_')[1:]
        if all([char in species_string_split for char in self._stocked_characteristics]):
            to_return_boolean = True
        else:
            to_return_boolean = False

        self._stocked_characteristics = set()
        return to_return_boolean


class Specific_Species_Operator(Bool_Override):
    """
        This class creates objects from the species strings from the meta-species to pass them to the rate functions
        as arguments. It uses the Bool_Override class to return true or false to the .dot operation inside the rate
        functions


        :param _stocked_characteristics: (str) stocks the characteristics of the queries performed by the user
        :param species_string: (str) string value from an individual species in MobsPY format
        :param species_object: (Species) Meta-species object which originated the meta-species str
    """

    def __init__(self, species_string, species_object):
        """
            Constructs the object from the species strings from the meta-species to pass them to the rate functions
            as arguments.

            :param species_string: (str) A string from MobsPy meta-species format
            :param species_object: (Species) = Meta-species set for which the species_string is contained in
        """
        self.species_string = species_string
        self._stocked_characteristics = set()
        self._species_object = species_object

    def __getattr__(self, characteristic):
        """
            Stores the characteristics for the boolean query inside the rate function by adding them to the set

            :param characteristic: (str) characteristic being queried
        """
        self._stocked_characteristics.add(characteristic)
        return self

    def __str__(self):
        """
            Returns the species_string from the MobsPy meta-species used in the object construction
        """
        return self.species_string

    def is_a(self, reference):
        """
            This function checks to see if the meta-species the species_string belong to has inherited from the
            parameter reference (reminder: every meta-species inherits from itself)

            :param reference: (Species) Meta-species object
            :return: (bool) True if the meta-species in Specific_Species_Operator has inherited from the reference
            False otherwise
        """
        if not self._stocked_characteristics:
            if reference in self._species_object.get_references():
                return True
            else:
                return False
        else:
            simlog.error('Concatenation of is_a and dot operator not supported. Please use them separately',
                         stack_index=2)

    def add(self, characteristic):
        """
            Adds characteristic to the set of characteristics when checking for all
            referred characteristics by the user

            :param characteristic: (str) characteristic to add to the set
        """
        self._stocked_characteristics.add(characteristic)

    def get_name(self):
        """
            Returns: The name of the species the reactant is in string format
        """
        return self._species_object.get_name()

    def get_characteristics(self):
        """
            Returns: The characteristics of the species in this given state
        """
        return set(self.species_string.split('_dot_')[1:])

    def get_state(self):
        """
            Returns: the string of a state with dots instead of _dot_
        """
        return self.species_string.replace('_dot_', '.')


class ExpressionDefiner:
    """
        This class defines a MobsPy expression. It encompasses what is necessary for something to act as an expression.
        For instance, the functionalities of storing the operations an object has gone through
    """
    @classmethod
    def execute_op(cls, first, second, operation):
        """
            Executes a unit operation using the Quantity class

            :param first: first number
            :param second: second number
            :param operation: operation to be executed
        """
        # TODO: FIX HERE NOW! FABRICIO DONT FORGET THIS DUMBASS
        # CHANGED HERE
        #if isinstance(first, OverrideQuantity) or isinstance(second, OverrideQuantity):
        #    raise Exception('Override quantity is not supposed to be in execute op')

        q_object = None
        if isinstance(first, Quantity) and isinstance(second, Quantity):
            if operation == '__add__':
                q_object = Quantity.__add__(first, second)
            elif operation == '__radd__':
                q_object = Quantity.__radd__(first, second)
            elif operation == '__sub__':
                q_object = Quantity.__sub__(first, second)
            elif operation == '__rsub__':
                q_object = Quantity.__rsub__(first, second)
            elif operation == '__mul__':
                q_object = Quantity.__mul__(first, second)
            elif operation == '__rmul__':
                q_object = Quantity.__rmul__(first, second)
            elif operation == '__truediv__':
                q_object = Quantity.__truediv__(first, second)
            elif operation == '__rtruediv__':
                q_object = Quantity.__rtruediv__(first, second)
        else:
            if operation == '__add__':
                q_object = first + second
            elif operation == '__radd__':
                q_object = second + first
            elif operation == '__sub__':
                q_object = first - second
            elif operation == '__rsub__':
                q_object = second - first
            elif operation == '__mul__':
                q_object = first*second
            elif operation == '__rmul__':
                q_object = second*first
            elif operation == '__truediv__':
                q_object = first/second
            elif operation == '__rtruediv__':
                q_object = second/first
            elif operation == '__pow__':
                q_object = first**second
            elif operation == '__rpow__':
                q_object = second**first

        if q_object is not None:
            return q_object
        else:
            raise TypeError('Non Valid Operation resulted in no q_object creation')

    def execute_quantity_op(self, other, operation):
        """
            Executes the unit based operation. It is used to verify the unit of a MobsPy expression

            :param other: other number (or expression) to execute the operation on
            :param operation: string symbol of the operation to be executed
        """
        try:
            if isinstance(other, ExpressionDefiner):
                count_op = self.execute_op(self._unit_count_op, other._unit_count_op, operation)
            else:
                count_op = self.execute_op(self._unit_count_op, other, operation)
        except Exception as e:
            # raise e - For debug
            count_op = e

        try:
            if isinstance(other, ExpressionDefiner):
                conc_op = self.execute_op(self._unit_conc_op, other._unit_conc_op, operation)
            else:
                conc_op = self.execute_op(self._unit_conc_op, other, operation)
        except Exception as e:
            # raise e - For debug
            conc_op = e

        return count_op, conc_op

    def __array_ufunc__(self, ufunc, _, *inputs):
        """
            This function handles the numpy compatibility with MobsPy expressions. ufunc is the operation,
            and inputs are both the numpy element (arrays not yet supported) and the quantity object.
            You have to reverse the direction, for Numpy the numpy object is the mail call, for us it is in reverse!
            That's why the r in every operation

            :param ufunc: numpy operation
            :param _: method, not used, please don't erase or it becomes an input - will break function
            :param inputs: numpy element and the quantity object
        """
        # Implement all numpy operations
        if ufunc == np_add:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return self.__radd__(inputs[0])
            else:
                simlog.error('MobsPy does not yet support array-wise numpy operations, only element-wise')
        elif ufunc == np_subtract:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return self.__rsub__(inputs[0])
            else:
                simlog.error('MobsPy does not yet support array-wise numpy operations, only element-wise')
        elif ufunc == np_multiply:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return self.__rmul__(inputs[0])
            else:
                simlog.error('MobsPy does not yet support array-wise numpy operations, only element-wise')
        elif ufunc == np_divide:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return self.__rtruediv__(inputs[0])
            else:
                simlog.error('MobsPy does not yet support array-wise numpy operations, only element-wise')
        else:
            simlog.error('Numpy operation not yet supported by MobsPy')

    # T is here to avoid problems with __getattr__ from units and quantities that were not overridden
    def __add__(self, other):
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, '__add__')
            return self.create_from_new_operation(other, '+', count_op, conc_op, True)
        else:
            return self.non_expression_add(other)

    def __radd__(self, other):
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, '__radd__')
            return self.create_from_new_operation(other, '+', count_op, conc_op, False)
        else:
            return self.non_expression_radd(other)

    def __sub__(self, other):
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, '__sub__')
            return self.create_from_new_operation(other, '-', count_op, conc_op, True)
        else:
            return self.non_expression_sub(other)

    def __rsub__(self, other):
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, '__rsub__')
            return self.create_from_new_operation(other, '-', count_op, conc_op, False)
        else:
            return self.non_expression_rsub(other)

    def __mul__(self, other):
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, '__mul__')
            return self.create_from_new_operation(other, '*', count_op, conc_op, True)
        else:
            return self.non_expression_mul(other)

    def __rmul__(self, other):
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, '__rmul__')
            return self.create_from_new_operation(other, '*', count_op, conc_op, False)
        else:
            return self.non_expression_rmul(other)

    def __truediv__(self, other):
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, '__truediv__')
            return self.create_from_new_operation(other, '/', count_op, conc_op, True)
        else:
            return self.non_expression_truediv(other)

    def __rtruediv__(self, other):
        if self._ms_active:
            count_op, conc_op = self.execute_quantity_op(other, '__rtruediv__')
            return self.create_from_new_operation(other, '/', count_op, conc_op, False)
        else:
            return self.non_expression_rtruediv(other)

    def __pow__(self, other):
        # if type(other) != int and type(other) != float:
        #    raise TypeError('Power must only be int or float')
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, '__pow__')
            return self.create_from_new_operation(other, '^', count_op, conc_op)
        else:
            return self.non_expression_pow(other)

    def __rpow__(self, other):
        # if type(other) != int and type(other) != float:
        #    raise TypeError('Power must only be int or float')
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, '__rpow__')
            return self.create_from_new_operation(other, '^', count_op, conc_op, False)
        else:
            return self.non_expression_rpow(other)

    def combine_binary_attributes(self, other, attribute):
        """
            Or gates to binary True or False attributes from self and other

            :param other: other expression to execute operation
            :param attribute: attributed to be combined
        """
        to_return = False
        try:
            if self.__dict__[attribute]:
                to_return = True
        except KeyError:
            pass
        except AttributeError:
            pass

        try:
            if other.__dict__[attribute]:
                to_return = True
        except KeyError:
            pass
        except AttributeError:
            pass

        return to_return

    def _generate_necessary_attributes(self):
        """
            This function was implemented as a replacement for innit for classes that inherit from multiple
            Pint objects. This gives the object all the necessary attributes to execute create_from_new_operation.
            I've done this instead of overrinding the __init__ because it had some compatibility issue with Pint
            __init__ at the time of writting this
        """
        # Operation variables
        self._operation = None
        self._unit_count_op = 1
        self._unit_conc_op = 1

        # MobsPy active - Behavior change
        self._ms_active = False

        # Parameter Variables
        self._parameter_set = set()
        self._expression_variables = set()

        # Presence of units
        self._has_units = False

        # Concentration and counts
        self._count_in_model = False
        self._concentration_in_model = False
        self._count_in_expression = False
        self._concentration_in_expression = False

        # Dimension
        self._dimension = None

    def create_from_new_operation(self, other, symbol, count_op, conc_op, direct_sense=True,
                                  operation=None):
        """
            Executes the storing based operation. It also executes a unit operation to store the verify the unit
            of the given expression

            :param other: other number (or expression) to execute the operation on
            :param symbol: string symbol of the operation to be executed
            :param count_op: unit of the expression if the arguments are considered dimentionless
            :param count_op: unit of the expression if the arguments are considered 1/v
            :param direct_sense: sense of the operation
            :param operation: current operation in the stack
        """
        _has_units = False
        try:
            # Returns string True not boolean to avoid risk __getattr__ bugs
            if self._has_units == 'T':
                _has_units = 'T'
            if other._has_units == 'T':
                _has_units = 'T'
        except AttributeError:
            pass

        try:
            c1 = self._has_units == 'T'
        except AttributeError:
            c1 = False
        try:
            c2 = other._has_units == 'T'
        except AttributeError:
            c2 = False
        if c1 or c2:
            if isinstance(count_op, Exception) and isinstance(conc_op, Exception):
                simlog.error("The units cannot be resolved. "
                             "Ex: This error is cased by summing two number with different units "
                             "(1/u.s) + 1*(u.l/u.s), \nor other impossible unit operations", stack_index=3)

        if isinstance(self, Quantity) or isinstance(self, OverrideQuantity):
            self = QuantityConverter.convert_received_unit(self)
        if isinstance(other, Quantity) or isinstance(other, OverrideQuantity):
            other = QuantityConverter.convert_received_unit(other)

        try:
            if isinstance(count_op, Quantity):
                count_op = QuantityConverter.convert_received_unit(count_op)
        except Exception as e:
            count_op = e

        try:
            if isinstance(conc_op, Quantity):
                conc_op = QuantityConverter.convert_received_unit(conc_op)
        except Exception as e:
            conc_op = e

        if type(self._operation) == str or (isinstance(other, ExpressionDefiner) and
                                            type(other._operation)) == str:
            int_op_self = str(self.magnitude) if isinstance(self, Quantity) or isinstance(self, OverrideQuantity)\
                else str(self)
            int_op_other = str(other.magnitude) if isinstance(other, Quantity) or isinstance(other, OverrideQuantity)\
                else str(other)

            if direct_sense:
                op1, op2 = int_op_self, int_op_other
            else:
                op2, op1 = int_op_self, int_op_other

            if operation is None:
                operation = '(' + op1 + symbol + op2 + ')'
        else:
            operation = count_op.magnitude

        new_parameter_set = self._parameter_set
        new_expression_variables = self._expression_variables

        try:
            new_parameter_set = new_parameter_set.union(other._parameter_set)
        except AttributeError:
            pass
        try:
            new_expression_variables = new_expression_variables.union(other._expression_variables)
        except AttributeError:
            pass

        if isinstance(other, ExpressionDefiner):
            new_parameter_set = self._parameter_set.union(other._parameter_set)
        elif isinstance(other, MobsPyExpression):
            self._expression_variables.union(other._expression_variables)

        _count_in_model = self.combine_binary_attributes(other, '_count_in_model')
        _concentration_in_model = self.combine_binary_attributes(other, '_concentration_in_model')
        _count_in_expression = self.combine_binary_attributes(other, '_count_in_expression')
        _concentration_in_expression = self.combine_binary_attributes(other, '_concentration_in_expression')

        if _count_in_model and _concentration_in_model:
            simlog.error('A meta-species in a model cannot be both a count and a concentration')

        if _count_in_expression and _concentration_in_expression:
            simlog.error('A meta-species in an expression cannot be both a count and a concentration')

        for variable in self._expression_variables:
            variable._operation = variable.species_string

        dimension_1 = None
        if isinstance(self, MobsPyExpression):
            dimension_1 = self._dimension

        dimension_2 = None
        if isinstance(other, MobsPyExpression):
            dimension_2 = other._dimension

        if dimension_1 is not None and dimension_2 is None:
            dimension = dimension_1
        elif dimension_1 is None and dimension_2 is not None:
            dimension = dimension_2
        elif dimension_1 is not None and dimension_2 is not None:
            if dimension_1 != dimension_2:
                raise TypeError('Dimensions are inconsistent between different MobsPy expression objects')
            else:
                dimension = dimension_1
        else:
            dimension = None

        return MobsPyExpression(species_string='$Null', species_object=None,
                                operation=operation, unit_count_op=count_op, unit_conc_op=conc_op,
                                dimension=dimension,
                                expression_variables=new_expression_variables, parameter_set=new_parameter_set,
                                count_in_model=_count_in_model, concentration_in_model=_concentration_in_model,
                                count_in_expression=_count_in_expression,
                                concentration_in_expression=_concentration_in_expression,
                                has_units=_has_units)


class OverrideUnitRegistry:
    """
        It was necessary to override Pint's unit registry so it behaves one way under context and
        normally without context
    """
    def __init__(self):
        self.unit_registry_object = UnitRegistry()
        self._ms_active = False

    def __getattr__(self, item):
        # Convert here
        if item == 'h':
            item = 'hour'

        q_object = 1*self.unit_registry_object.__getattr__(item)
        owq_obj = OverrideQuantity(q_object)
        owq_obj._ms_active = self._ms_active
        return owq_obj

# u is defined here
u = OverrideUnitRegistry()
def set_u_context():
    u._ms_active = True

def reset_u_context():
    u._ms_active = False

class QuantityConverter:
    """
        Class that converts an unit from any to MobsPy L-s-counts conversion
    """
    @classmethod
    def convert_received_unit(cls, quantity):
        """
            Converts a received quantity to L-s-counts, standard MobsPy units

            :param quantity: received quantity to convert
        """
        ur = u.unit_registry_object

        is_override = False
        if isinstance(quantity, OverrideQuantity):
            is_override = True
            copied_quantity = deepcopy(quantity.q_object)
        else:
            copied_quantity = deepcopy(quantity)

        to_convert_into = str(copied_quantity.dimensionality)

        if '[length]' in to_convert_into:
            to_convert_into = to_convert_into.replace('[length]', 'dm')
        if '[time]' in to_convert_into:
            to_convert_into = to_convert_into.replace('[time]', 's')
        if '[temperature]' in to_convert_into:
            to_convert_into = to_convert_into.replace('[temperature]', 'K')
        if '[mass]' in to_convert_into:
            to_convert_into = to_convert_into.replace('[mass]', 'kg')
        if '[substance]' in to_convert_into:
            index = to_convert_into.find('[substance]')
            d_string = str(quantity.dimensionality)
            try:
                mol_power = int(d_string[index + 15])
            except ValueError:
                mol_power = 1
            except IndexError:
                mol_power = 1
            except Exception as e:
                raise e

            copied_quantity = copied_quantity*(N_A / (1*ur.mol)) ** mol_power
            if '[substance]' in copied_quantity.dimensionality:
                copied_quantity = copied_quantity * ((1 * ur.mol) / N_A) ** (2*mol_power)
            if '[substance]' in copied_quantity.dimensionality:
                raise TypeError('Could not convert molar quantity')

            to_convert_into = to_convert_into.replace('[substance]', '1')

        copied_quantity.ito(to_convert_into)

        if not is_override:
            return copied_quantity
        else:
            temp = OverrideQuantity(copied_quantity)
            temp.set_ms_active(quantity._ms_active)
            return temp


class OverrideQuantity(ExpressionDefiner, Quantity):

    def __array_ufunc__(self, ufunc, _, *inputs):
        """
            This function handles the numpy compatibility with MobsPy expressions. ufunc is the operation,
            and inputs are both the numpy element (arrays not yet supported) and the quantity object.

            :param ufunc: numpy operation
            :param _: method, not used, please don't erase or it becomes an input - will break function
            :param inputs: numpy element and the quantity object
        """
        # Implement all numpy operations
        if ufunc == np_add:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return OverrideQuantity(float(inputs[0]) + self.q_object)
            else:
                simlog.error('MobsPy does not yet support array-wise numpy operations, only element-wise')
        elif ufunc == np_subtract:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return OverrideQuantity(float(inputs[0]) - self.q_object)
            else:
                simlog.error('MobsPy does not yet support array-wise numpy operations, only element-wise')
        elif ufunc == np_multiply:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return OverrideQuantity(float(inputs[0])*self.q_object)
            else:
                simlog.error('MobsPy does not yet support array-wise numpy operations, only element-wise')
        elif ufunc == np_divide:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return OverrideQuantity(float(inputs[0])/self.q_object)
            else:
                simlog.error('MobsPy does not yet support array-wise numpy operations, only element-wise')
        else:
            simlog.error('Numpy operation not yet supported by MobsPy')

    def non_expression_add(self, other):
        q_object = Quantity.__add__(self.q_object, other.q_object)
        return OverrideQuantity(q_object)

    def non_expression_radd(self, other):
        q_object = Quantity.__radd__(self.q_object, other.q_object)
        return OverrideQuantity(q_object)

    def non_expression_sub(self, other):
        q_object = Quantity.__sub__(self.q_object, other.q_object)
        return OverrideQuantity(q_object)

    def non_expression_rsub(self, other):
        q_object = Quantity.__rsub__(self.q_object, other.q_object)
        return OverrideQuantity(q_object)

    def non_expression_mul(self, other):
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__mul__(self.q_object, other.q_object)
        else:
            q_object = Quantity.__mul__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_rmul(self, other):
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__rmul__(self.q_object, other.q_object)
        else:
            q_object = Quantity.__rmul__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_truediv(self, other):
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__truediv__(self.q_object, other.q_object)
        else:
            q_object = Quantity.__truediv__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_rtruediv(self, other):
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__rtruediv__(self.q_object, other.q_object)
        else:
            q_object = Quantity.__rtruediv__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_pow(self, other):
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__pow__(self.q_object, other.q_object)
        else:
            q_object = Quantity.__pow__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_rpow(self, other):
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__rpow__(self.q_object, other.q_object)
        else:
            q_object = Quantity.__rpow__(self.q_object, other)
        return OverrideQuantity(q_object)

    def __init__(self, quantity_object):

        self._generate_necessary_attributes()
        self.q_object = quantity_object

        self._unit_count_op = quantity_object
        self._unit_conc_op = quantity_object

        for key, item in quantity_object.__dict__.items():
            self.__dict__[key] = item

        self._operation = self.q_object.magnitude
        self._expression_variables = set()
        self._parameter_set = set()
        self._has_units = 'T'

    def __str__(self):
        if self._ms_active:
            return str(self._operation)
        else:
            return str(self.q_object)

    def set_ms_active(self, ms_active):
        """
            ms_active is the boolean responsible for activating the context. This function can be used to set it

            :param ms_active: boolean with the state to set the expression or normal unit context
        """
        self._ms_active = ms_active

    # Had to propose new functions for conversion as Pint is doing something strange and checking the returned object
    def convert(self, unit):
        new_q_object = self.q_object.to(unit)
        temp = OverrideQuantity(new_q_object)
        temp._ms_active = self._ms_active
        return temp

    def convert_into(self, unit):
        self.q_object.ito(unit)


class MobsPyExpression(Specific_Species_Operator, ExpressionDefiner):
    """
        MobsPy expression objects are passed to rate functions to store the expressions they've been through
    """
    def __str__(self):
        return str(self._operation)

    def __init__(self, species_string, species_object, operation=None,
                 unit_count_op=1, unit_conc_op=(1/u.unit_registry_object.liter),
                 dimension=None,
                 expression_variables=None, parameter_set=None,
                 count_in_model=True, concentration_in_model=False, count_in_expression=True,
                 concentration_in_expression=False, has_units=False):

        super().__init__(species_string, species_object)
        self._generate_necessary_attributes()

        self._ms_active = True

        self._dimension = dimension

        if expression_variables is None:
            self._expression_variables = set()
            self._expression_variables.add(self)
        else:
            self._expression_variables = expression_variables

        if operation is None:
            self._operation = self.species_string
        else:
            self._operation = operation

        if parameter_set is None:
            self._parameter_set = set()
        else:
            self._parameter_set = parameter_set

        self._count_in_model = count_in_model
        self._concentration_in_model = concentration_in_model
        self._count_in_expression = count_in_expression
        self._concentration_in_expression = concentration_in_expression

        self._unit_count_op = unit_count_op
        self._unit_conc_op = unit_conc_op

        if self._count_in_model_in_model:
            self._unit_conc_op = None

        if self._concentration_in_model:
            self._unit_conc_op = None

        self._has_units = has_units

    def __getattr__(self, item):
        return super().__getattr__(item)

    # Write string operation return and unit
    def generate_string_operation(self, skip_check=False, reaction_order=None):
        """
            Converts the expression from what has been stored to a string format for the sbml file

            :param skip_check: skip units check or not - always set to False, True only for debugging
            :param reaction_order: order of the reaction - to check if the unit is correct
        """
        ur = u.unit_registry_object

        if self._dimension is None:
            dimension = 3
        else:
            dimension = self._dimension

        operation = str(self._operation)

        if skip_check:
            return operation, True

        if self._has_units != 'T':
            self._count_in_expression = True

        if self._has_units == 'T' and self._expression_variables == set() and reaction_order is not None:
            if self._unit_count_op.units == (1 / ur.second):
                return operation, True
            elif self._unit_conc_op.units == (1*ur.decimeter**(3*reaction_order)/(ur.second*ur.decimeter**dimension)):
                return operation, False

        c1 = self._has_units == 'T'
        c2 = isinstance(self._unit_count_op, Exception)
        c3 = isinstance(self._unit_conc_op,  Exception)
        if c1 and (c2 and c3):
            raise Exception('The unit for the reaction did not resolve due to an illegal operation' + '\n' +
                            'With meta-species as count: ' + str(self._unit_count_op) + '\n' +
                            'With meta-species as concentration: ' + str(self._unit_conc_op) + '\n')

        # If both count and concentration are valid, count takes priority
        if c1 and not c3:
            if self._unit_conc_op.units == \
                    (1 /(u.unit_registry_object.second*u.unit_registry_object.decimeter**dimension)):
                self._concentration_in_expression = True

        if c1 and not c2:
            if self._unit_count_op.units == (1/u.unit_registry_object.second):
                self._count_in_expression = True

        if self._has_units == 'T' and not self._count_in_expression and not self._concentration_in_expression:
            raise TypeError('The automated unit detection resulted was unable to compile the reaction rate. ' +
                            'With meta-species as count: ' + str(self._unit_count_op) + '. ' +
                            'With meta-species as concentration: ' + str(self._unit_conc_op) + '. '
                            'Please check if the resulting units are valid. '
                            'For counts the expression must result in units 1/[time] and for concentration in '
                            '1/([time][volume]).')

        convert_operation = str(self._operation)

        for variable in self._expression_variables:
            count_name = '$count$' + variable.species_string
            concentration_name = '$concentration$' + variable.species_string

            default_name = variable.species_string
            replace_name = variable.species_string

            if self._count_in_model:
                convert_operation = replace_spe_in_expr(convert_operation, count_name, replace_name)
                convert_operation = replace_spe_in_expr(convert_operation, concentration_name,
                                                        '(' + replace_name + '/volume)')

                # Add unit check!!!
                if self._count_in_expression:
                    pass
                elif self._concentration_in_expression:
                    convert_operation = replace_spe_in_expr(convert_operation, default_name,
                                                            '(' + replace_name + '/volume)')
                    convert_operation = '(' + convert_operation + ')' + '*volume'
                else:
                    raise ValueError('The expression did not resolve for lack of concentration/count specifications')
            elif self._concentration_in_model:
                convert_operation = replace_spe_in_expr(convert_operation, concentration_name, replace_name)
                convert_operation = convert_operation.replace(count_name, '(' + replace_name + '*volume)')
                convert_operation = replace_spe_in_expr(convert_operation, count_name,'(' + replace_name + '*volume)')

                if self._count_in_expression:
                    convert_operation = convert_operation.replace(default_name, '(' + replace_name + '*volume)')
                    convert_operation = replace_spe_in_expr(convert_operation, default_name,
                                                            '(' + replace_name + '*volume)')
                    convert_operation = '(' + convert_operation + ')' + '/volume'
                elif self._concentration_in_expression:
                    pass
                else:
                    raise ValueError('The expression did not resolve for lack of concentration/count specifications')

        return convert_operation, self._count_in_expression


def check_if_non_expression_operated(other):
    if not isinstance(other, ExpressionDefiner) and not type(other) == int and not type(other) == float \
            and not isinstance(other, Quantity) and not isinstance(other, (np_int_, np_float_)):
        other = MobsPyExpression('($asg_' + str(other) + ')', None, dimension=None, count_in_model=True,
                                 concentration_in_model=False, count_in_expression=False,
                                 concentration_in_expression=False)
    return other


def replace_spe_in_expr(string, to_replace, replacement):
    pattern = re.compile(re.escape(to_replace) + r'(?![a-zA-Z0-9_])')
    return pattern.sub(replacement, string)

class _Count_Base:

    'self.species_string'
    def __getitem__(self, item):
        try:
            for v in item._expression_variables:
                item._operation = item._operation.replace(v.species_string, '$count$' + v.species_string)
            return item
        except AttributeError:
            simlog.error('Count[] operator can only be used in MobsPy expressions')


Count = _Count_Base()


class _Conc_Base:

    def __getitem__(self, item):
        try:
            for v in item._expression_variables:
                item._operation = item._operation.replace(v.species_string, '$concentration$' + v.species_string)
            return item
        except AttributeError:
            simlog.error('Concentration[] operator can only be used in MobsPy expressions')
        pass


Concentration = _Conc_Base()

if __name__ == '__main__':

    print('Ran local script')

    # u = OverrideUnitRegistry(is_active=True)
    # a = MobsPyExpression('a_dot_alive', None, count_in_model=False, concentration_in_model=True,
    #                      count_in_expression=False, concentration_in_expression=True)

    # Don't create new objects
    # Quantity Object
    # To much unit registry

    x = MobsPyExpression('A', None, dimension=3,
                         count_in_model=True,
                         concentration_in_model=False,
                         count_in_expression=False,
                         concentration_in_expression=False)

    r = x**2.8
    print(r.generate_string_operation())
