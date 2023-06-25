from pint import Quantity, UnitRegistry, errors
import mobspy.simulation_logging.log_scripts as simlog
from scipy.constants import N_A


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


class ExpressionDefiner:

    def __add__(self, other):
        if self._ms_active:
            return self.create_from_new_operation(other, '+', True)
        else:
            return self.non_expression_add(other)

    def __radd__(self, other):
        if self._ms_active:
            return self.create_from_new_operation(other, '+', False)
        else:
            return self.non_expression_radd(other)

    def __sub__(self, other):
        if self._ms_active:
            return self.create_from_new_operation(other, '-', True)
        else:
            return self.non_expression_sub(other)

    def __rsub__(self, other):
        if self._ms_active:
            return self.create_from_new_operation(other, '-', False)
        else:
            return self.non_expression_rsub(other)

    def __mul__(self, other):
        if self._ms_active:
            return self.create_from_new_operation(other, '*', True)
        else:
            return self.non_expression_mul(other)

    def __rmul__(self, other):
        if self._ms_active:
            return self.create_from_new_operation(other, '*', False)
        else:
            return self.non_expression_rmul(other)

    def __truediv__(self, other):
        if self._ms_active:
            return self.create_from_new_operation(other, '/', True)
        else:
            return self.non_expression_truediv(other)

    def __rtruediv__(self, other):
        if self._ms_active:
            return self.create_from_new_operation(other, '/', False)
        else:
            return self.non_expression_rtruediv(other)

    def __pow__(self, other):
        if self._ms_active:
            operation = '(' + self._operation + ')' + f'^{other}'
            if isinstance(other, Quantity):
                simlog.error('Power operation does not allow for units in the exponent')
            unit_operation = '(' + self._unit_operation + ')' + f'^{other}'
            return self.create_from_new_operation(other, '^', True, operation, unit_operation)
        else:
            return self.non_expression_pow(other)

    def combine_binary_attributes(self, other, attribute):
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

    def create_from_new_operation(self, other, symbol, direct_sense=False,
                                  operation=None, unit_operation=None):

        if direct_sense:
            op1, op2 = str(self), str(other)
        else:
            op2, op1 = str(self), str(other)

        if operation is None:
            operation = '(' + op1 + symbol + op2 + ')'

        if unit_operation is None:
            if isinstance(other, MobsPyExpression):
                if direct_sense:
                    unit_operation = '(' + self._unit_operation + symbol + other._unit_operation + ')'
                else:
                    unit_operation = '(' + other._unit_operation + symbol + self._unit_operation  + ')'
            else:
                if direct_sense:
                    unit_operation = '(' + self.__dict__['_unit_operation'] + symbol + '1' + ')'
                else:
                    unit_operation = '(' + '1' + symbol + self.__dict__['_unit_operation'] + ')'

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

        _has_units = False
        try:
            # Returns string True not boolean to avoid risk __getattr__ bugs
            if self._has_units == 'T':
                _has_units = 'T'
            if other._has_units == 'T':
                _has_units = 'T'
        except AttributeError:
            pass

        _count_in_model = self.combine_binary_attributes(other, '_count_in_model')
        _concentration_in_model = self.combine_binary_attributes(other, '_concentration_in_model')
        _count_in_expression = self.combine_binary_attributes(other, '_count_in_expression')
        _concentration_in_expression = self.combine_binary_attributes(other, '_concentration_in_expression')

        if _count_in_model and _concentration_in_model:
            simlog.error('An expression cannot be both a count and a concentration')

        if _count_in_expression and _concentration_in_expression:
            simlog.error('An expression cannot both output a count and a concentration')

        for variable in self._expression_variables:
            variable._operation = '$' + variable.species_string

        return MobsPyExpression('$Null', None, operation, unit_operation,
                                new_expression_variables, new_parameter_set,
                                _count_in_model, _concentration_in_model, _count_in_expression,
                                _concentration_in_expression, _has_units)


class MobsPyExpression(Specific_Species_Operator, ExpressionDefiner):

    def __str__(self):
        return str(self._operation)

    def __init__(self, species_string, species_object, operation=None, unit_based_operation=None,
                 expression_variables=None, parameter_set=None,
                 count_in_model=True, concentration_in_model=False, count_in_expression=True, concentration_in_expression=False,
                 has_units=False):
        super().__init__(species_string, species_object)

        self._ms_active = True

        if expression_variables is None:
            self._expression_variables = set()
            self._expression_variables.add(self)
        else:
            self._expression_variables = expression_variables

        if operation is None:
            self._operation = '$' + self.species_string
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

        if unit_based_operation is None:
            self._unit_operation = '$' + self.species_string
        else:
            self._unit_operation = unit_based_operation

        self._has_units = has_units

    def __getattr__(self, item):
        if item not in ['count', 'concentration']:
            return super().__getattr__(item)

        if item == 'count':
            self._operation = '$' + 'count' + self._operation
        elif item == 'concentration':
            self._operation = '$' + 'concentration' + self._operation

        return self

    @classmethod
    def automatic_unit_deduction(cls, unit_expression, variable_name,
                                 expression_is_count=False, expression_is_concentration=False):
        copied_expression = str(unit_expression)

        try:
            test_conc = copied_expression.replace('$' + variable_name, '(1/decimeter**3)')
            r1 = UnitRegistry().parse_expression(test_conc)
            if r1.units == '1 / decimeter / second':
                expression_is_concentration = True
            elif expression_is_concentration:
                simlog.error('Placeholder Error')
        except errors.DimensionalityError:
            pass

        try:
            test_conc = copied_expression.replace('$' + variable_name, '1')
            r2 = UnitRegistry().parse_expression(test_conc)
            if r2.units == '1 / second':
                expression_is_count = True
            elif expression_is_count:
                simlog.error('Placeholder Error')
        except errors.DimensionalityError:
            pass

        if not expression_is_count and not expression_is_concentration:
            simlog.error('Placeholder Error')
        elif expression_is_count and expression_is_concentration:
            expression_is_count = True
            expression_is_concentration= False

        return expression_is_count, expression_is_concentration
        
    # Write string operation return and unit
    def return_string_operation(self):

        convert_operation = str(self._operation)
        convert_unit_op = str(self._unit_operation)

        for variable in self._expression_variables:
            count_name = '$count$' + variable.species_string
            concentration_name = '$concentration$' + variable.species_string
            default_name = '$' + variable.species_string
            replace_name = variable.species_string

            convert_unit_op = convert_unit_op.replace(count_name, '1')
            convert_unit_op = convert_unit_op.replace(concentration_name, '(1/decimeter**3)')

            # Add correct condition
            self._count_in_expression, self._concentration_in_expression = \
                    self.automatic_unit_deduction(convert_unit_op, variable.species_string)

            if self._count_in_model:
                convert_operation = convert_operation.replace(count_name, replace_name)
                convert_operation = convert_operation.replace(concentration_name, '(' + replace_name + '/volume)')

                # Add unit check!!!
                if self._count_in_expression:
                    convert_operation = convert_operation.replace(default_name, replace_name)
                elif self._concentration_in_expression:
                    convert_operation = convert_operation.replace(default_name, '(' + replace_name + '/volume)')
                    convert_operation = '(' + convert_operation + ')' + '*volume'
                else:
                    simlog.error('Placeholder Error')
            elif self._concentration_in_model:
                convert_operation = convert_operation.replace(concentration_name, replace_name)
                convert_operation = convert_operation.replace(count_name, '(' + replace_name  + '*volume)')

                if self._count_in_expression:
                    convert_operation = convert_operation.replace(default_name, '(' + replace_name + '*volume)')
                    convert_operation = '(' + convert_operation + ')' + '/volume'
                elif self._concentration_in_expression:
                    convert_operation = convert_operation.replace(default_name, replace_name)
                else:
                    simlog.error('Placeholder Error')

        return convert_operation

    def compile_resulting_unit(self):
        pass


class QuantityConverter:

    @classmethod
    def convert_received_unit(cls, quantity):
        to_convert_into = str(quantity.dimensionality)
        mult_by_NA = False
        div_by_NA = False

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
            mult_by_NA = False
            div_by_NA = False
            if to_convert_into[index - 2] == '*':
                mult_by_NA = True
            elif to_convert_into[index - 2] == '/':
                div_by_NA = True

            # In case substance is the first one
            if index - 2 < 0:
                mult_by_NA = True

            to_convert_into = to_convert_into.replace('[substance]', 'mole')

        converted_quantity = quantity.to(to_convert_into)
        mag = converted_quantity.magnitude
        unit = str(converted_quantity.units)

        if mult_by_NA:
            mag = mag * N_A
        if div_by_NA:
            mag = mag / N_A
        if mult_by_NA or div_by_NA:
            unit = str(unit).replace('mole', '1')

        return str(mag), '(' + unit + ')'


class OverrideUnitRegistry(ExpressionDefiner, UnitRegistry, QuantityConverter):

    def non_expression_add(self, other):
        q_object = super(OverrideUnitRegistry, self).__add__(other)
        return OverrideQuantity(q_object)

    def non_expression_radd(self, other):
        q_object = super(OverrideUnitRegistry, self).__radd__(other)
        return OverrideQuantity(q_object)

    def non_expression_sub(self, other):
        q_object = super(OverrideUnitRegistry, self).__sub__(other)
        return OverrideQuantity(q_object)

    def non_expression_rsub(self, other):
        q_object = super(OverrideUnitRegistry, self).__rsub__(other)
        return OverrideQuantity(q_object)

    def non_expression_mul(self, other):
        q_object = super(OverrideUnitRegistry, self).__mul__(other)
        return OverrideQuantity(q_object)

    def non_expression_rmul(self, other):
        q_object = super(OverrideUnitRegistry, self).__rmul__(other)
        return OverrideQuantity(q_object)

    def non_expression_truediv(self, other):
        q_object = super(OverrideUnitRegistry, self).__truediv__(other)
        return OverrideQuantity(q_object)

    def non_expression_rtruediv(self, other):
        q_object = super(OverrideUnitRegistry, self).__rtruediv__(other)
        return OverrideQuantity(q_object)

    def non_expression_pow(self, other):
        q_object = super(OverrideUnitRegistry, self).__pow__(other)
        return OverrideQuantity(q_object)

    def __init__(self, operation=None, unit_operation=None, is_active=False):
        super(OverrideUnitRegistry, self).__init__()
        self._ms_active = is_active

        if operation is None:
            self._operation = '1'
        else:
            self._operation = operation

        if unit_operation is None:
            self._unit_operation = '1'
        else:
            self._unit_operation = unit_operation

        self._expression_variables = set()
        self._parameter_set = set()
        self._has_units = 'T'

    def __getattr__(self, item):
        # Convert here
        if self._ms_active:
            if item == '_has_units':
                return 'T'

            self._ms_active = False
            unit_rec = 1 * self.parse_units(item)
            self._ms_active = True

            mag, unit = self.convert_received_unit(unit_rec)
            return OverrideUnitRegistry(mag, unit, is_active=True)
        else:
            unit_rec = 1 * self.parse_units(item)
            return OverrideQuantity(unit_rec)

    def __str__(self):
        if self._ms_active:
            return str(self._operation)
        else:
            return super(OverrideUnitRegistry, self).__str__()


class OverrideQuantity(ExpressionDefiner, Quantity, QuantityConverter):

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
            q_object = Quantity.__rmul__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_rtruediv(self, other):
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__rtruediv__(self.q_object, other.q_object)
        else:
            q_object = Quantity.__rmul__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_pow(self, other):
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__rtruediv__(self.q_object, other.q_object)
        else:
            q_object = Quantity.__rmul__(self.q_object, other)
        return OverrideQuantity(q_object)

    def __init__(self, quantity_object, ms_active=False):
        self.q_object = quantity_object

        for key, item in quantity_object.__dict__.items():
            self.__dict__[key] = item

        self._ms_active = ms_active

        mag, unit = self.convert_received_unit(self.q_object)

        # Adjust this with magnitude and unit
        self._operation = mag
        self._unit_operation = unit

        self._expression_variables = set()
        self._parameter_set = set()
        self._has_units = 'T'

    def __getattr__(self, item):
        return self.q_object.__getattr__(item)

    def __get__(self, item):
        return self.q_object.__get__(item)

    def __str__(self):
        if self._ms_active:
            return self._operation
        else:
            return str(self.q_object)


if __name__ == '__main__':

    # u = OverrideUnitRegistry(is_active=True)
    # a = MobsPyExpression('a_dot_alive', None, count_in_model=False, concentration_in_model=True,
    #                      count_in_expression=False, concentration_in_expression=True)

    u = OverrideUnitRegistry()

    test_q = 10 * u.molar
    a = MobsPyExpression('a_dot_alive', None, count_in_model=True, concentration_in_model=False,
                         count_in_expression=False, concentration_in_expression=False)
    K = 5*u.molar
    K._ms_active = True
    c = u.h*u.decimeter**3
    c._ms_active = True
    r = 1/c*1/(1 + K/a)
    a = r.return_string_operation()
    print(a)


