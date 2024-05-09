from mobspy.simulation_logging.log_scripts import error as simlog_error
from inspect import stack as inspect_stack
# import mobspy.modules.species_string_generator as ssg
from mobspy.modules.species_string_generator import construct_all_combinations as ssg_construct_all_combinations
from pint import Quantity

"""
    This model implements the logic comparison structure for events and conditional durations.
    The syntax follows:( MetaSpeciesObject Operator Value ) ex: (A > 5)
    For clause combination one has - & and | or.
    All clauses must be inside parenthesis when combined due to Python operator precedence.
"""


class SpeciesComparator:
    """
        This class implements the comparisons necessary for events and conditional durations for Species
        ex (A > 5)

        :param _simulation_context: (Simulation) current simulation under context
    """

    @classmethod
    def check_parenthesis(cls, code_line, line_number, pos, symbol, number_of_comp):
        """
            This function compiles part of a logic expression for meta-species. It checks to see if an operator '<' or
            '>' is properly written inside parenthesis. It only applies if more than one operator is present in the
            code line

        :param code_line: (str) line of code to be compiled
        :param line_number: (str) number of the code line currently being compiled
        :param pos: (int) position of the '>' or '<' operator
        :param symbol: '(' or ')', indicates in which direction the string analysis should proceded
        :param number_of_comp: for distinction the default case where only one operator is present

        :raise simlog.error: if the code line is not properly written isolating the clauses with parenthesis

        :return: (bool) false if the line compiles, raises a error if it does not
        """
        i = pos
        condition_not_satisfied = True
        while condition_not_satisfied:
            try:
                if symbol == ')':
                    i += 1
                elif symbol == '(':
                    i -= 1

                # Check if finished
                if i == 0 or i == len(code_line):
                    if number_of_comp > 1:
                        simlog_error(f'At: {code_line} \n' + f'Line number: {line_number} \n' +
                                     'All clauses must be individually isolated by parenthesis '
                                     '- Ex: (A <= 0) & (B <= 0) \n')
                    else:
                        return False

                char = code_line[i]
                if char == '<' or char == '>' or char == '|' or char == '&':
                    simlog_error(f'At: {code_line} \n' + f'Line number: {line_number} \n' +
                                 'All clauses must be isolated by parenthesis - Ex: (A <= 0) & (B <= 0) \n')
                elif char == ')' and number_of_comp > 1:
                    if symbol == ')':
                        condition_not_satisfied = False
                elif char == '(' and number_of_comp > 1:
                    if symbol == '(':
                        condition_not_satisfied = False
            except IndexError:
                simlog_error(f'Error Compiling the following line {line_number}: {code_line}')
        return condition_not_satisfied

    @classmethod
    def compile_code_line(cls):
        """
            Compiles the code line executing a logical operation

            :raise simlog.error: if the code line is not properly written isolating the clauses with parenthesis
        """
        code_line = inspect_stack()[2].code_context[0][:-1]
        line_number = inspect_stack()[2].lineno
        if 'and' in code_line or 'or' in code_line:
            simlog_error(f'At: {code_line} \n' + f'Line number: {line_number} \n' +
                         'Event notation did not compile, please use & for \'and\' and  | for \'or\' \n'
                         'Please also put the clauses under parentheses: example (A <= 0) & (B <= 0)')

        temp_code_line = code_line.replace(' ', '')
        multiplication_position = [pos for pos, char in enumerate(temp_code_line) if char == '*']
        for pos in multiplication_position:
            if not (temp_code_line[pos - 1].isnumeric() or temp_code_line[pos + 1].isnumeric()):
                simlog_error(f'At: {code_line} \n' + f'Line number: {line_number} \n' +
                             'Multiplication between meta-species under comparison context' +
                             ' not yet supported by MobsPy - \n.')

        number_of_comp = code_line.count('<') + code_line.count('>')
        comparison_position = [pos for pos, char in enumerate(code_line) if char == '<' or char == '>']
        for pos in comparison_position:
            for symbol in ['(', ')']:
                assert not cls.check_parenthesis(code_line, line_number, pos, symbol, number_of_comp)

    def __init__(self):
        self._simulation_context = None

    def add_operation_and_number(self, symbol, number):
        """
            Creates a MetaSpeciesLogicResolver from the comparison of a meta-species with a value or another
            meta-species

            :param symbol: Comparison symbol '>=', '<=', '>' or '<'
            :param number: (int) if compared to a value, (Species or ReactingSpecies) if compared to the objects
        """
        number = self.reformat_number_and_species(number)
        if type(number) == list:
            operation = [self.logical_add_species(self)] + [symbol] + number
        else:
            operation = [self.logical_add_species(self)] + [symbol] + [number]
        logic_re_object = MetaSpeciesLogicResolver(operation, self._simulation_context)
        return logic_re_object

    def reformat_number_and_species(self, number):
        """
            Discovers if the number is a Species, ReactingSpecies or integer and prepares the output accordingly

            :param number: (int) if compared to a value, (Species or ReactingSpecies) if compared to the objects
        """
        if isinstance(number, SpeciesComparator):
            if number.is_species():
                number = self.logical_add_species(number)
            else:
                number = self.logical_add_reacting_species(number)
        return number

    @classmethod
    def logical_add_species(cls, species):
        """
            Adds a species object to an event trigger operation

            :param species: (Species) species object
        """

        return {'object': species, 'characteristics': set()}

    @classmethod
    def logical_add_reacting_species(cls, react_spe):
        """
            Adds a reacting species object to an event trigger operation by summing over all the indicated
            characteristics. It also accepts sums of species multiplied by integers

            :param react_spe: (ReactingSpecies) ReactingSpecies object
        """
        operation = []

        react_spe.check_context()
        for i, react_dict in enumerate(react_spe.list_of_reactants):
            if i > 0:
                dl = ['+']
            else:
                dl = []
            dl = dl + [react_dict['stoichiometry'], '*']
            dl = dl + [{'object': react_dict['object'], 'characteristics': react_dict['characteristics']}]
            operation = operation + dl
        return dl

    def __lt__(self, number):
        self.compile_code_line()
        return self.add_operation_and_number('<', number)

    def __le__(self, number):
        self.compile_code_line()
        return self.add_operation_and_number('<=', number)

    def __gt__(self, number):
        self.compile_code_line()
        return self.add_operation_and_number('>', number)

    def __ge__(self, number):
        self.compile_code_line()
        return self.add_operation_and_number('>=', number)

    def __eq__(self, other):
        if self._simulation_context is not None:
            simlog_error('Equality assignment not allowed for event condition in MobsPy.\n'
                         'Please if necessary use ( >= ) & ( =< )')
        else:
            return id(self) == id(other)

    def __neg__(self, other):
        return id(self) != id(other)

    def __hash__(self):
        return hash(id(self))


class ReactingSpeciesComparator(SpeciesComparator):
    """
        This class implements the comparisons necessary for events and conditional durations for Reacting Species
        ex (A.a1 > 5)

        :param _simulation_context: (Simulation) current simulation under context
    """

    def __init__(self):
        super(ReactingSpeciesComparator, self).__init__()
        self._simulation_context = None

    def check_context(self):
        """
            Checks if a reacting species is under context and add it to the attribute self._simulation_context
        """
        for react_dict in self.list_of_reactants:
            if react_dict['object']._simulation_context is not None:
                self._simulation_context = react_dict['object']._simulation_context
                break
        else:
            self._simulation_context = None

    def add_operation_and_number(self, symbol, number):
        """
            Adds the symbol and number to create a MetaSpeciesLogicResolver object

            :param symbol: Comparative symbols '<=', '<', '>=', or '>'
            :param number: (int) if compared to a value, (Species or ReactingSpecies) if compared to the objects

            :returns: MetaSpeciesLogicResolver object containing the comparison
        """
        number = self.reformat_number_and_species(number)
        if type(number) == list:
            operation = self.logical_add_reacting_species(self) + [symbol] + number
        else:
            operation = self.logical_add_reacting_species(self) + [symbol] + [number]
        logic_re_object = MetaSpeciesLogicResolver(operation, self._simulation_context)
        return logic_re_object


class MetaSpeciesLogicResolver:
    """
        This object stores the logical and comparison operations that become the event triggers for conditions
        or the conditional duration for simulations

        :param operation: (list) the logical operation that will be transformed in a string for copasi. Format ex
        [{object: ..., characteristics:...}, '<=', 10]
        :param simulation_context: (Simulation) simulation under context
    """

    def __init__(self, operation, simulation_context=None):
        self.operation = operation
        self.simulation_context = simulation_context

    def __and__(self, other):
        return self._join(other, '&&')

    def __or__(self, other):
        return self._join(other, '||')

    def _join(self, other, symbol):
        if not isinstance(other, MetaSpeciesLogicResolver):
            print('ERROR')
            exit()

        new_operation = ['('] + ['('] + self.operation + [')'] + [symbol] + ['('] + other.operation + [')'] + [')']
        self.operation = new_operation
        return self

    def __lt__(self, number):
        self.add_double_symbol('<', number)
        return self

    def __le__(self, number):
        self.add_double_symbol('<=', number)
        return self

    def __gt__(self, number):
        self.add_double_symbol('>', number)
        return self

    def __ge__(self, number):
        self.add_double_symbol('>=', number)
        return self

    def add_double_symbol(self, symbol, number):
        if type(number) != int or type(number) != float or not isinstance(number, Quantity):
            self.operation = [number, symbol] + self.operation

    @classmethod
    def find_all_species_strings(cls, species, characteristics, species_for_sbml):
        reference_set = set(characteristics)
        reference_set.add(str(species))
        return [x for x in species_for_sbml.keys() if reference_set.issubset(set(x.split('_dot_')))]

    def generate_string(self, characteristics_to_object, to_sort=False):
        """
            When a meta-species is used in a logic expression, this function transforms a meta-species in the sum of
            all individual states

            :param characteristics_to_object: orthogonal characteristic space
            :param to_sort: sort strings or not - so the sum will always apear in the same order
        """
        copasi_str = ''
        for i, e in enumerate(self.operation):
            if type(e) == int or type(e) == float or type(e) == str:
                copasi_str = copasi_str + str(e) + ' '
            else:
                copasi_str = copasi_str + '('
                ite = ssg_construct_all_combinations(e['object'], e['characteristics'],
                                                     characteristics_to_object, '_dot_')
                if to_sort:
                    ite = sorted(ite)

                for j, species_str in enumerate(ite):
                    if j == 0:
                        copasi_str = copasi_str + f'{species_str}'
                    else:
                        copasi_str = copasi_str + f' + {species_str}'
                copasi_str = copasi_str + ')' + ' '

        return copasi_str
