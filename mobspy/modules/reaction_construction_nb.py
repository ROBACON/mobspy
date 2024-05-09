"""
    This module is responsible for the construction of all individual reactions from a meta-reaction
"""
from copy import deepcopy
from mobspy.modules.function_rate_code import extract_reaction_rate as fr_extract_reaction_rate
from itertools import product as itertools_product
from mobspy.modules.meta_class import Reactions
from mobspy.modules.species_string_generator import construct_all_combinations as ssg_construct_all_combinations
from inspect import signature as inspect_signature
from mobspy.modules.order_operators import Default
from mobspy.simulation_logging.log_scripts import error as simlog_error
# from mobspy.modules.mobspy_parameters import *
from mobspy.modules.context_related_scripts import Unit_Context_Setter as crs_Unit_Context_Setter
from mobspy.modules.meta_class_utils import count_string_dictionary as mcu_count_string_dictionary


def iterator_for_combinations(list_of_lists):
    """
        Iterate through all the combinations of a list of lists
        [[A,B,C], [D, F, E], [G]] means:
        ADG, AFG, AEG, BDG, BFG, BEG, CDG, CFG, CEG .....

        Parameter:
            list_of_lists (list of lists) = list of all lists to iterate through all combinations
    """
    for i in itertools_product(*list_of_lists):
        yield i


def copy_reaction(reaction):
    """
        Just copies a meta-reaction
        Deepcopy is not working because it calls __getattr__ on the species with a private method

        Parameters:
            reaction (meta-reaction) = meta-reaction to be copied
    """
    reactants = []
    for reactant in reaction.reactants:
        characteristics = deepcopy(reactant['characteristics'])
        species = reactant['object']
        stoichiometry = reactant['stoichiometry']
        reactants.append({'object': species, 'characteristics': characteristics,
                          'stoichiometry': stoichiometry})

    products = []
    for product in reaction.products:
        characteristics = deepcopy(product['characteristics'])
        species = product['object']
        stoichiometry = product['stoichiometry']
        products.append({'object': species, 'characteristics': characteristics,
                         'stoichiometry': stoichiometry})

    reaction_copy = Reactions(reactants, products)
    reaction_copy.set_rate(reaction.rate)
    return reaction_copy


def check_for_invalid_reactions(reactions, ref_characteristics_to_object):
    """
        If query references two independent characteristics from the in meta-species it would result in an empty set
        Like Ecoli.live.dead - (assuming live and dead belong to the same meta-species characteristics set)
        If that ever happens inside a meta-reaction we just pop an error

        :param reactions: (set of meta-reactions) set of meta-reactions inside the model
        :param ref_characteristics_to_object: (dict) dictionary with the characteristics as keys and their respective
        object as values
    """
    for reaction in reactions:
        for reactant in reaction.reactants:

            check_for_duplicates = {}
            for cha in reactant['characteristics']:

                if cha == 'all$':
                    continue

                try:
                    check_for_duplicates[ref_characteristics_to_object[cha]]
                    simlog_error(f'Illegal reaction: {reaction}. \n'
                                 'There is a query with two characteristics '
                                 f'{ref_characteristics_to_object[cha].get_characteristics()} from the same'
                                 f' vector axis resulting in an impossible query. \n'
                                 f'As all those characteristics have been directly added to'
                                 f'{ref_characteristics_to_object[cha]}, they are located in the same vector axis. \n'
                                 f'To solve this either assign the characteristics to new meta-species and use '
                                 f'inheritance or create a new reaction for each desired query.')
                except KeyError:
                    try:
                        check_for_duplicates[ref_characteristics_to_object[cha]] = cha
                    except KeyError:
                        simlog_error(
                            f'A base object for characteristic {cha} was not found in the species supplied to the '
                            f'simulator \n'
                            'Perhaps a species is missing ? ')

        for product in reaction.products:

            check_for_duplicates = {}
            for cha in product['characteristics']:

                try:
                    check_for_duplicates[ref_characteristics_to_object[cha]]
                    simlog_error(f'Illegal reaction: {reaction}. \n'
                                 'There is a transformation with two characteristics '
                                 f'{ref_characteristics_to_object[cha].get_characteristics()} from the same'
                                 f' vector axis resulting in an undefinable transformation. \n'
                                 f'As all those characteristics have been directly added to '
                                 f'{ref_characteristics_to_object[cha]}, they are located in the same vector axis. \n'
                                 f'To solve this either assign the characteristics to new meta-species and use '
                                 f'inheritance or create a new reaction for each desired query.')
                except KeyError:
                    if '$' not in cha:
                        check_for_duplicates[ref_characteristics_to_object[cha]] = cha


def construct_reactant_structures(reactant_species, ref_characteristics_to_object):
    """
        This function finds the corresponding strings according to the reactant meta-species with or without a query
        pack then in a list and return

        :param reactant_species: (meta-species object) species objects of the involved species
        :param ref_characteristics_to_object: (dict) dictionary with characteristics as keys and species objects as
        values
    """
    species_string_combinations = []

    for reactant in reactant_species:
        species_string_combinations.append(ssg_construct_all_combinations(reactant['object'],
                                                                          reactant['characteristics'],
                                                                          ref_characteristics_to_object))

    return species_string_combinations


def construct_order_structure(species_order_list, current_species_string_list):
    """
        Order structure for reaction order operations. Returns the cyclic_dictionary to be used by the order operator.
        The meta-species objects are the keys of this dictionary and a lists of species strings currently being used
        in the reaction are the values - Allowing the product to find it's corresponding species-string in a future
        step

        :param species_order_list: (list of meta-species objects) list of meta-species objects as they appear in the
        meta-reaction
        :param current_species_string_list: (list of strings) list of strings in MobsPy format of the species
        currently in this specific reaction

        :return: cyclic_dict (dict) Dictionary where the keys are meta-species objects and the values are
        lists of species
    """
    cyclic_dict = {}
    for species_object, species_string in zip(species_order_list, current_species_string_list):

        try:
            cyclic_dict[species_object].append(species_string)
        except KeyError:
            cyclic_dict[species_object] = [species_string]

    return cyclic_dict


def construct_product_structure(reaction):
    """
        This function unpacks the products in a meta-reaction

        :param: reaction meta-reaction currently being analysed

        :return: product_list = A list of dictionaries for each product with the meta-species object, the label and the
        characteristics
    """
    product_list = []
    for product in reaction.products:
        for _ in range(product['stoichiometry']):
            product_list.append({'species': product['object'], 'label': product['label'],
                                 'characteristics': product['characteristics']})

    return product_list


def construct_single_reaction_for_sbml(reactant_species_string_list, product_species_string_list, reaction_rate):
    """
        This function constructs the reactions for SBML for the conversion by the model builder script
        It follows the following structure 're':[('stoichmetry', reactantant_string) ....
        The reaction rate must be a string containing the reaction kinetics
        This returns a single reaction to be appended by the reactions_for_sbml dictionary

        :param reactant_species_string_list: (list of strings) list of reactants in MobsPy format
        :param product_species_string_list: (list of strings) list of products in MobsPy format
        :param reaction_rate: (str) reaction rate expression as a string

        :return: to_return (dict) = dictionary that packs the reactants products and rate
    """
    to_return = {'re': [], 'pr': [], 'kin': reaction_rate}

    reactant_count_dict = mcu_count_string_dictionary(reactant_species_string_list)
    product_count_dict = mcu_count_string_dictionary(product_species_string_list)

    for key in reactant_count_dict:
        to_return['re'].append((reactant_count_dict[key], key))

    for key in product_count_dict:
        to_return['pr'].append((product_count_dict[key], key))

    return to_return


def get_involved_species(reaction, meta_species_in_model):
    """
        This extracts all the involved meta-species inside a reaction
        This function is responsible for implementing the inheritance mechanism, by finding within each meta-species
        references set, if they reference the meta-species in the reaction

        :param reaction: (meta-reaction object)
        :param meta_species_in_model: (list) list of meta-species used in the model

        :return: base_species_order (list of meta-species objects) = order that the meta-species appear in the
        meta-reaction, reactant_species_combination_list (list of lists of meta-species) = list of lists of all
        meta-species that have inherited from the meta-species in the meta-reaction
    """
    reactant_species_combination_list = []
    base_species_order = []

    for reactant in reaction.reactants:
        flag_absent_reactant = False
        for _ in range(reactant['stoichiometry']):

            species_for_reactant = []
            base_species_order.append((reactant['object'], reactant['label']))

            for species in meta_species_in_model:
                if reactant['object'] in species.get_references():
                    species_for_reactant.append({'object': species,
                                                 'characteristics': reactant['characteristics'],
                                                 'stoichiometry': reactant['stoichiometry']})
                    flag_absent_reactant = True

            if not flag_absent_reactant:
                simlog_error(f'Species {reactant["object"]} or any inheritors were not found in model \n'
                             f'For reaction {reaction} \n'
                             f'Please add the species or remove the reaction')

            reactant_species_combination_list.append(species_for_reactant)

    return base_species_order, reactant_species_combination_list


def construct_rate_function_arguments(rate_function, reaction):
    rate_function_arguments = str(inspect_signature(rate_function))

    black_list = ['*', '=']
    if any(i in rate_function_arguments for i in black_list):
        simlog_error(f'Rate arguments must not contain = or *. \n'
                     f'Error in reaction {reaction}. \n'
                     f'Error in rate function {rate_function} in signature {str(inspect_signature(rate_function))}')

    rate_function_arguments = str(rate_function_arguments).replace('(', '')
    rate_function_arguments = str(rate_function_arguments).replace(')', '')
    rate_function_arguments = str(rate_function_arguments).replace(' ', '')
    rate_function_arguments = rate_function_arguments.split(',')
    return rate_function_arguments


def create_all_reactions(reactions, meta_species_in_model,
                         ref_characteristics_to_object,
                         type_of_model, dimension, parameter_exist, parameters_in_reaction,
                         skip_check):
    """
        This function creates all reactions
        Returns the reactions_for_sbml and parameters_for_sbml dictionary
        Those will be used by another module to create the SBML file

        :param reactions: (meta-reaction objects) reactions objects constructed by the meta_class module
        :param meta_species_in_model: (list) list of meta-species in model
        :param ref_characteristics_to_object: (dict) Characteristics as keys objects as values
        :param type_of_model: (str) stochastic or deterministic
        :param dimension: (int) model dimension 1D, 2D, 3D, .....

        :returns: reactions_for_sbml (dict) = dictionary with all reactions that will be added to the sbml model file,
        parameters_for_sbml (dict) = parameters for the sbml model file
    """
    reactions_for_sbml = {}

    check_for_invalid_reactions(reactions, ref_characteristics_to_object)


    # Initiate expressions
    with crs_Unit_Context_Setter():
        for reaction in reactions:

            base_species_order, reactant_species_combination_list = get_involved_species(reaction,
                                                                                         meta_species_in_model)

            for combination_of_reactant_species in iterator_for_combinations(reactant_species_combination_list):

                reactant_species_string_combination_list = \
                    construct_reactant_structures(combination_of_reactant_species, ref_characteristics_to_object)

                for reactant_string_list in iterator_for_combinations(reactant_species_string_combination_list):

                    product_object_list = construct_product_structure(reaction)
                    order_structure = construct_order_structure(base_species_order, reactant_string_list)

                    if reaction.order is None:
                        product_species_species_string_combination_list = Default(order_structure, product_object_list,
                                                                                  meta_species_in_model,
                                                                                  ref_characteristics_to_object)
                    else:
                        product_species_species_string_combination_list = reaction.order(order_structure,
                                                                                         product_object_list,
                                                                                         meta_species_in_model,
                                                                                         ref_characteristics_to_object)

                    for product_string_list in \
                            iterator_for_combinations(product_species_species_string_combination_list):

                        reaction_rate_arguments = None
                        if callable(reaction.rate):
                            reaction_rate_arguments = construct_rate_function_arguments(reaction.rate, reaction)

                        reactant_strings = ['_dot_'.join([reactant[0].get_name()] + reactant[1:])
                                            if len(reactant) > 1 else reactant[0].get_name()
                                            for reactant in reactant_string_list]

                        try:
                            rate_string, parameters_in_reaction = \
                                fr_extract_reaction_rate(combination_of_reactant_species,
                                                         reactant_strings,
                                                         reaction.rate, type_of_model,
                                                         dimension,
                                                         reaction_rate_arguments,
                                                         parameter_exist,
                                                         parameters_in_reaction,
                                                         skip_check)
                        except TypeError as e:
                            simlog_error(f'On reaction {reaction} \n' + str(e))

                        if rate_string == 0:
                            continue

                        reactions_for_sbml['reaction_' + str(len(reactions_for_sbml))] = \
                            construct_single_reaction_for_sbml(reactant_strings,
                                                               product_string_list,
                                                               rate_string)

    return reactions_for_sbml, parameters_in_reaction


if __name__ == '__main__':
    pass
