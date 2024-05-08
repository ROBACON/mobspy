from mobspy.modules.mobspy_expressions import MobsPyExpression
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.modules.species_string_generator as ssg
from mobspy.modules.mobspy_expressions import u
import re


class Assignment_Operator:
    _asg_context = False
    regex_pattern = r'\(\$arg(?:\.[^\s().]+)?\)(?=[^\s()]|$)'

    @staticmethod
    # To be updated for a regex, it was not working and I need to finish my Ph.D.
    def find_arg_strings(input_string):

        arg_strings = []
        flag_found = False
        stack = ''
        for char in input_string:
            if char == '$':
                flag_found = True

            if not flag_found:
                continue
            else:
                if char != ')':
                    stack += char
                else:
                    flag_found = False
                    arg_strings.append(stack)
                    stack = ''

        return arg_strings

    def __enter__(self):
        self._asg_context = True
        yield 0

    def set_context(self):
        self._asg_context = True

    def __exit__(self, *args):
        self._asg_context = False
        yield 0

    def reset_context(self):
        self._asg_context = False

    def check_context(self):
        return self._asg_context

    @staticmethod
    def check_arguments(first, second):

        if type(first) == int or type(first) == float:
            first = MobsPyExpression(str(first), None, dimension=None, count_in_model=True,
                                     concentration_in_model=False, count_in_expression=False,
                                     concentration_in_expression=False)

        if type(second) == int or type(second) == float:
            second = MobsPyExpression(str(second), None, dimension=None, count_in_model=True,
                                      concentration_in_model=False, count_in_expression=False,
                                      concentration_in_expression=False)

        if not isinstance(first, MobsPyExpression):
            first = MobsPyExpression('($asg_' + str(first) + ')', None, dimension=None, count_in_model=True,
                                     concentration_in_model=False, count_in_expression=False,
                                     concentration_in_expression=False)
        if not isinstance(second, MobsPyExpression):
            second = MobsPyExpression('($asg_' + str(second) + ')', None, dimension=None, count_in_model=True,
                                      concentration_in_model=False, count_in_expression=False,
                                      concentration_in_expression=False)
        return first, second

    @staticmethod
    def add(first, second):
        first, second = Assignment_Operator.check_arguments(first, second)
        return first + second

    @staticmethod
    def sub(first, second):
        first, second = Assignment_Operator.check_arguments(first, second)
        return first - second

    @staticmethod
    def mul(first, second):
        first, second = Assignment_Operator.check_arguments(first, second)
        return first * second

    @staticmethod
    def div(first, second):
        first, second = Assignment_Operator.check_arguments(first, second)
        return first / second

    @staticmethod
    def pow(first, second):
        first, second = Assignment_Operator.check_arguments(first, second)
        return first ** second

    @staticmethod
    def generate_replacement_in_expression(express_spe, ortogonal_vector_structure, meta_species_in_model,
                                           expression_tuple):
        spe_str = express_spe.replace('$asg_', '')
        spe_str = spe_str.split('.')

        # CHECK HERE FOR MISSING SPECIES IN MODEL
        for meta_spe in meta_species_in_model:
            if spe_str[0] == str(meta_spe):
                spe_object = meta_spe
                break
        else:
            spe_name = str(expression_tuple[0][0])
            error_message = f'Assignment {spe_name}, {expression_tuple[0][1]}: ' \
                            f'{expression_tuple[1].replace("$asg_", "")} failed\n' \
                            f'One of the meta-species in the assignment expression was not found in the model'
            simlog.error(error_message)

        if len(spe_str) == 1:
            str_comb = ssg.construct_all_combinations(spe_object, set(), ortogonal_vector_structure, '_dot_')
        else:
            str_comb =ssg.construct_all_combinations(spe_object, set(spe_str[1:]), ortogonal_vector_structure, '_dot_')

        for i, e in enumerate(str_comb):
            if i == 0:
                to_replace = e
            else:
                to_replace += '+' + e
        return to_replace

    @staticmethod
    def compile_assignments_for_sbml(unprocessed_asgns, ortogonal_vector_structure, meta_species_in_model):

        assignments_for_sbml = {}
        for i, asg in enumerate(unprocessed_asgns):
            spe_to_asgn = ssg.construct_species_char_list(asg[0], asg[1],
                                                          ortogonal_vector_structure,
                                                          symbol='_dot_')

            asgn_expression = str(unprocessed_asgns[asg])
            spe_in_expression = Assignment_Operator.find_arg_strings(asgn_expression)

            replacement_dict = {}
            for spe in spe_in_expression:
                replacement_dict[spe] = \
                    Assignment_Operator.generate_replacement_in_expression(spe, ortogonal_vector_structure,
                                                                           meta_species_in_model,
                                                                           (asg, asgn_expression))

            for key, item in replacement_dict.items():
                asgn_expression = asgn_expression.replace(key, item)

            assignments_for_sbml['assignment_' + str(i)] = {'species': spe_to_asgn, 'expression': asgn_expression}

        return assignments_for_sbml


Assign = Assignment_Operator()


class Asg:
    assignments = {}

    def __init__(self, meta_spe, species_or_reacting):
        Assign.set_context()
        self.meta_spe = []
        self.asgn_key = []
        if species_or_reacting:
            self.meta_spe.append(meta_spe)
            self.asgn_key.append((meta_spe, tuple()))
        else:
            for reacting_spe in meta_spe.list_of_reactants:
                self.meta_spe.append(reacting_spe['object'])
                self.asgn_key.append((reacting_spe["object"], tuple(reacting_spe["characteristics"])))
        self.species_or_reacting = species_or_reacting

    def __call__(self, assignment):
        for spe, key in zip(self.meta_spe, self.asgn_key):
            spe._assignments[key] = assignment
        Assign.reset_context()

    def __getattr__(self, item):
        simlog.error("Assignments must be the last query in the stack - Ex: A.young.blue.assign()", stack_index=2)








