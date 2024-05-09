"""
    mobspy.modules.process_result_data.py

    Handles rate_assignment to the reactions:
        The passing of arguments to be executed by the user supplied rate function
        The construction of mass action kinetics
        Handles the different types of rates supplied by the user
"""

from mobspy.modules.meta_class import Zero as mc_Zero
from mobspy.modules.meta_class_utils import count_string_dictionary as mcu_count_string_dictionary

from mobspy.simulation_logging.log_scripts import error as simlog_error, debug as simlog_debug
from mobspy.modules.unit_handler import convert_rate as uh_convert_rate
from mobspy.modules.mobspy_expressions import MobsPyExpression as mbe_MobsPyExpression, \
    ExpressionDefiner as mbe_ExpressionDefiner, Specific_Species_Operator as mbe_Specific_Species_Operator
from pint import Quantity
from mobspy.modules.mobspy_parameters import Mobspy_Parameter as mp_Mobspy_Parameter
from re import split as re_split
import timeit


def extract_reaction_rate(combination_of_reactant_species, reactant_string_list
                          , reaction_rate_function, type_of_model, dimension, function_rate_arguments,
                          parameter_exist, parameters_in_reaction, skip_check):
    """
        This function is responsible for returning the reaction rate string for the model construction. To do this it
        does a different action depending on the type of the reaction_rate_function (we consider constants as functions)
        It passes the rate as a string expression in MobsPy standard units

        :param  function_rate_arguments: list of strings of the function rate argument ex:['r1', 'r2', ....]
        :param  combination_of_reactant_species: (list of Species) Meta-species currently being used in this reaction
        :param  reactant_string_list: (list of strings) list of species strings in order they apear in the reaction
        :param  reaction_rate_function: (float, callable, Quantity) rate stored in the reaction object
        :param  dimension: (int) system dimension (for the rate conversion)

        :return: reaction_rate_string (str) the reaction kinetics as a string for SBML
    """
    is_count = False
    if type(reaction_rate_function) == int or type(reaction_rate_function) == float or \
            isinstance(reaction_rate_function, Quantity):
        # Function is a constant function number int here
        reaction_rate_function, dimension, is_count = uh_convert_rate(reaction_rate_function,
                                                                      len(reactant_string_list),
                                                                      dimension)

        if reaction_rate_function == 0:
            return 0, parameters_in_reaction
        reaction_rate_string = basic_kinetics_string(reactant_string_list,
                                                     reaction_rate_function, type_of_model, is_count)

    elif isinstance(reaction_rate_function, mbe_ExpressionDefiner) and parameter_exist:
        parameters_in_reaction = parameters_in_reaction.union(reaction_rate_function._parameter_set)
        reaction_rate_string = basic_kinetics_string(reactant_string_list,
                                                     reaction_rate_function, type_of_model, is_count)

    elif function_rate_arguments is not None:

        # [''] means that it is a function that takes no arguments (empty signature)
        if function_rate_arguments != ['']:
            arguments = prepare_arguments_for_callable(combination_of_reactant_species,
                                                       reactant_string_list, function_rate_arguments,
                                                       dimension)

            rate = reaction_rate_function(**arguments)

        else:
            rate = reaction_rate_function()

        rate, dimension, is_count = uh_convert_rate(rate, len(reactant_string_list), dimension)

        if rate == 0:
            return 0, parameters_in_reaction

        if type(rate) == int or type(rate) == float:
            reaction_rate_string = basic_kinetics_string(reactant_string_list,
                                                         rate, type_of_model, is_count)
        elif type(rate) == str:
            reaction_rate_string = rate
            # Remove dollar sign symbols from expression object
            reaction_rate_string = reaction_rate_string.replace('$', '')

            if parameter_exist:
                parameters_in_reaction = search_for_parameters_in_str(reaction_rate_string,
                                                                      parameter_exist, parameters_in_reaction)
        elif isinstance(rate, mbe_MobsPyExpression):

            # Having an expression variable implies it is a constructed expression - not mass action
            if len(rate._expression_variables) > 0:
                reaction_rate_string, _ = rate.generate_string_operation(skip_check=skip_check)

            # Having no expression variables implies it is a constant for mass - action kinetics.
            elif len(rate._expression_variables) == 0:
                rate_for_mass_action, is_count = \
                    rate.generate_string_operation(reaction_order=len(reactant_string_list), skip_check=skip_check)
                parameters_in_rate = rate._parameter_set
                parameters_in_reaction = parameters_in_reaction.union(parameters_in_rate)

                reaction_rate_string = basic_kinetics_string(reactant_string_list,
                                                             rate_for_mass_action, type_of_model, is_count)

        elif isinstance(rate, mp_Mobspy_Parameter):
            parameters_in_reaction.add(rate)
            reaction_rate_string = basic_kinetics_string(reactant_string_list,
                                                         str(rate), type_of_model)
        elif rate is None:
            simlog_error('There is a reaction rate missing for the following reactants: \n'
                         + str(reactant_string_list))
        else:
            simlog_error(f'The rate function {reaction_rate_function}, returned a non-valid value. \n' +
                         f'Only int, floats and str are accepted')
    elif reaction_rate_function is None:
        simlog_error('There is a reaction rate missing for the following reactants: \n'
                     + str(reactant_string_list))
    elif type(reaction_rate_function) == str:
        parameters_in_reaction = search_for_parameters_in_str(reaction_rate_function,
                                                              parameter_exist, parameters_in_reaction)
        reaction_rate_string = reaction_rate_function
    else:
        simlog_debug(type(reaction_rate_function))
        simlog_error(f'The {type(reaction_rate_function)}, from {reaction_rate_function} is not supported')

    return reaction_rate_string, parameters_in_reaction


def basic_kinetics_string(reactants, reaction_rate, type_of_model, is_count=False):
    """
        This constructs the bases for mass-action kinetics. Both for stochastic and deterministic depending on the type
        of model

        :params reactants: (list of str) list of reactants in MobsPy str format
        :params reaction_rate: (float) reaction constant
        :params type_of_model: (str) stochastic or deterministic - rate expressions are different depending on each case

        :return:  kinetics_string (str) mass action kinetics expression for the reaction
    """
    counts = mcu_count_string_dictionary(reactants)

    kinetics_string = ""
    for name, number in counts.items():
        if type_of_model.lower() == 'stochastic':
            kinetics_string += stochastic_string(name, number)
        elif type_of_model.lower() == 'deterministic':
            kinetics_string += deterministic_string(name, number)
        kinetics_string += ' * '

    kinetics_string += str(reaction_rate)

    n = 0
    for key, item in counts.items():
        n += item
    n = n - 1

    if n > 0 and not is_count:
        kinetics_string += f' * volume^{-n}'
    # n == -1 is the lowest possible value
    elif n < 0 and not is_count:
        kinetics_string += f' * volume'
    else:
        pass

    return kinetics_string


def stochastic_string(reactant_name, number):
    """
        This function returns the stochastic string expression for mass action kinetics
        For instance the reaction 2A -> 3A would imply A*(A-1)/2
        It only does so for one reactant, so it must be called for all reactants in the reaction

        :params reactant_name: (str) species string involved in the reaction
        :params number: (int) stoichiometry (number of times it appears)

        :return: to_return_string (str) the mass action kinetics string expression for only that species
    """
    to_return_string = ''
    for i in range(number):
        if i == 0:
            to_return_string = reactant_name
        else:
            to_return_string += f' * ({reactant_name} - {i})/{i + 1}'

    return to_return_string


def deterministic_string(reactant_name, number):
    """
        This function returns the deterministic string expression for mass action kinetics
        For instance the reaction 2A -> 3A would imply A*A
        It only does so for one reactant, so it must be called for all reactants in the reaction

        :params reactant_name: (str) species string involved in the reaction
        :params number: (int) stoichiometry (number of times it appears)

        :return: to_return_string (str) the mass action kinetics string expression for only that species
    """
    to_return_string = ''
    for i in range(number):
        if i == 0:
            to_return_string = reactant_name
        else:
            to_return_string += f' * {reactant_name}'
    return to_return_string


def prepare_arguments_for_callable(combination_of_reactant_species, reactant_string_list, rate_function_arguments,
                                   dimension):
    """
        This function prepares the requested arguments to the rate function for a given reaction by creating objects
        of the Specific_Species_Operator class

        :params combination_of_reactant_species: meta-species involved in the reaction
        :params reactant_string_list: species strings involved in the reaction
        :params rate_function_arguments: arguments received by the rate function
        :return: argument_dict - dictionary with arguments for a rate function
    """
    argument_dict = {}

    if rate_function_arguments is not None:
        i = 0
        for i, (species, reactant_string) in enumerate(zip(combination_of_reactant_species, reactant_string_list)):
            try:
                species = species['object']
                argument_dict[rate_function_arguments[i]] = mbe_MobsPyExpression(reactant_string, species,
                                                                                 dimension=dimension,
                                                                                 count_in_model=True,
                                                                                 concentration_in_model=False,
                                                                                 count_in_expression=False,
                                                                                 concentration_in_expression=False)
            except IndexError:
                continue

        if i != 0:
            while len(argument_dict) < len(rate_function_arguments):
                i += 1
                argument_dict[rate_function_arguments[i]] = mbe_Specific_Species_Operator('$Null', mc_Zero)
        elif i == 0:
            while len(argument_dict) < len(rate_function_arguments):
                argument_dict[rate_function_arguments[i]] = mbe_Specific_Species_Operator('$Null', mc_Zero)
                i += 1

    return argument_dict


def search_for_parameters_in_str(reaction_rate_string, parameters_exist, parameters_in_reaction):
    """
        Searches for MobsPy Parameter names in a string reaction rate. Uses the parameters_exit stack.
        If it finds a parameter it adds it to the set parameters_in_reaction

        :param reaction_rate_string: reaction rate in str format
        :param parameters_exist: stack of parameters available
        :param parameters_in_reaction: set of parameters already in reaction
        :return: parameters_in_reaction set of parameters already in reaction
    """
    split_operation = re_split(', |-|!|\*|\+|/|\)|\(| ', reaction_rate_string)
    split_operation = [x.replace(' ', '') for x in split_operation if x.replace(' ', '') != '']

    for name in split_operation:
        if name in parameters_exist:
            parameters_in_reaction.add(parameters_exist[name])

    return parameters_in_reaction


if __name__ == '__main__':
    pass
