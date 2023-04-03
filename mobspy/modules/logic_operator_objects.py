import mobspy.simulation_logging.log_scripts as simlog
import inspect
import mobspy.modules.species_string_generator as ssg
from pint import Quantity


class SpeciesComparator:

    @classmethod
    def check_parenthesis(cls, code_line, line_number, pos, symbol, number_of_comp):
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
                        simlog.error('All clauses must be individually isolated by parenthesis '
                                     '- Ex: (A <= 0) & (B <= 0) \n'
                                     + f'Problem in {line_number} : {code_line}')
                    else:
                        return False

                char = code_line[i]
                if char == '<' or char == '>' or char == '|' or char == '&':
                    simlog.error('All clauses must be isolated by parenthesis - Ex: (A <= 0) & (B <= 0) \n' +
                                 f'Problem in {line_number} : {code_line}')
                elif char == ')' and number_of_comp > 1:
                    if symbol == ')':
                        condition_not_satisfied = False
                elif char == '(' and number_of_comp > 1:
                    if symbol == '(':
                        condition_not_satisfied = False
            except IndexError:
                simlog.error(f'Error Compiling the following line {line_number}: {code_line}')
        return condition_not_satisfied

    @classmethod
    def compile_code_line(cls):
        code_line = inspect.stack()[2].code_context[0][:-1]
        line_number = inspect.stack()[2].lineno
        if 'and' in code_line or 'or' in code_line:
            simlog.error('Event notation did not compile, please use & for \'and\' and  | for \'or\' \n'
                         'Please also put the clauses under parentheses: example (A <= 0) & (B <= 0)')

        temp_code_line = code_line.replace(' ', '')
        multiplication_position = [pos for pos, char in enumerate(temp_code_line) if char == '*']
        for pos in multiplication_position:
            if not (temp_code_line[pos - 1].isnumeric() or temp_code_line[pos + 1].isnumeric()):
                simlog.error('Multiplication between meta-species under comparison context'
                             ' not yet supported by MobsPy - \n.' +
                             f'Problem in {line_number} : {code_line}')

        number_of_comp = code_line.count('<') + code_line.count('>')
        comparison_position = [pos for pos, char in enumerate(code_line) if char == '<' or char == '>']
        for pos in comparison_position:
            for symbol in ['(', ')']:
                assert not cls.check_parenthesis(code_line, line_number, pos, symbol, number_of_comp)

    def __init__(self):
        self._simulation_context = None

    def start_check(self):
        if self._simulation_context is not None:
            self._simulation_context.event_context_add()

    def add_operation_and_number(self, symbol, number):
        operation = [{'object': self, 'characteristics': set()}, symbol, number]
        logic_re_object = MetaSpeciesLogicResolver(operation, self._simulation_context)
        return logic_re_object

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
            simlog.error('Equality assignment not allowed for event condition in MobsPy.\n'
                         'Please if necessary use ( >= ) & ( =< )')
        else:
            return id(self) == id(other)

    def __neg__(self, other):
        return id(self) != id(other)

    def __hash__(self):
        return hash(id(self))


class ReactingSpeciesComparator(SpeciesComparator):

    def __init__(self):
        super(ReactingSpeciesComparator, self).__init__()
        self._simulation_context = None

    def start_check(self):
        self.check_context()
        if self._simulation_context is not None:
            self._simulation_context.event_context_add()

    def check_context(self):
        for react_dict in self.list_of_reactants:
            if react_dict['object']._simulation_context is not None:
                self._simulation_context = react_dict['object']._simulation_context
                break
        else:
            self._simulation_context = None

    def add_operation_and_number(self, symbol, number):
        operation = []

        self.check_context()
        for i, react_dict in enumerate(self.list_of_reactants):
            if i > 0:
                dl = ['+']
            else:
                dl = []
            dl = dl + [react_dict['stoichiometry'], '*']
            dl = dl + [{'object': react_dict['object'], 'characteristics': react_dict['characteristics']}]
            operation = operation + dl
        operation = operation + [symbol, number]
        logic_re_object = MetaSpeciesLogicResolver(operation, self._simulation_context)

        return logic_re_object


class MetaSpeciesLogicResolver:

    def __init__(self, operation, simulation_context=None):
        self.operation = operation
        self.simulation_context = simulation_context
        if self.simulation_context is not None:
            simulation_context.trigger_list.append(self)

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
        if self.simulation_context is not None:
            self.simulation_context.trigger_list.append(self.operation)

    @classmethod
    def find_all_species_strings(cls, species, characteristics, species_for_sbml):
        reference_set = set(characteristics)
        reference_set.add(str(species))
        return [x for x in species_for_sbml.keys() if reference_set.issubset(set(x.split('_dot_')))]

    def generate_string(self, characteristics_to_object, to_sort=False):
        copasi_str = ''
        for i, e in enumerate(self.operation):
            if type(e) == int or type(e) == float or type(e) == str:
                copasi_str = copasi_str + str(e) + ' '
            else:
                copasi_str = copasi_str + '('
                ite = ssg.construct_all_combinations(e['object'], e['characteristics'],
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
