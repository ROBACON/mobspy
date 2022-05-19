"""
    This module is responsible for the construction of all individual reactions from a meta-reaction
"""
from copy import deepcopy
import mobspy.modules.function_rate_code as fr
import itertools
import mobspy.modules.meta_class as mc
import mobspy.simulation_logging.log_scripts as simlog


def iterator_for_combinations(list_of_lists):
    """
        Iterate through all the combinations of a list of lists
        [[A,B,C], [D, F, E], [G]] means:
        ADG, AFG, AEG, BDG, BFG, BEG, CDG, CFG, CEG .....

        Parameter:
            list_of_lists (list of lists) = list of all lists to iterate through all combinations
    """
    for i in itertools.product(*list_of_lists):
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

    reaction_copy = mc.Reactions(reactants, products)
    reaction_copy.set_rate(reaction.rate)
    return reaction_copy


def check_for_invalid_reactions(reactions, ref_characteristics_to_object):
    """
        If query references two independent characteristics from the in meta-species it would result in an empty set
        Like Ecoli.live.dead - (assuming live and dead belong to the same meta-species characteristics set)
        If that ever happens inside a meta-reaction we just pop an error

        Parameters:
            reactions (set of meta-reactions) = set of meta-reactions inside the model
            ref_characteristics_to_object (dict) = dictionary where the characteristics are the keys and the values
            are the meta-species objects that they have been directly added to
    """
    for reaction in reactions:
        for reactant in reaction.reactants:

            check_for_duplicates = {}
            for cha in reactant['characteristics']:

                # Ok I love try catches for checking if something is in a dict. Don't judge me
                try:
                    check_for_duplicates[ref_characteristics_to_object[cha]]
                    simlog.error('Illegal reaction, there is a product with multiple '
                                 'characteristics from the same Species referenced \n'
                                 'Please divide the reaction accordingly \n'
                                 f'Error in {ref_characteristics_to_object[cha]} both {cha} and '
                                 f'{check_for_duplicates[ref_characteristics_to_object[cha]]}'
                                 f' referenced at the same time \n'
                                 f'The characteristics: {ref_characteristics_to_object[cha].get_characteristics()}')
                except KeyError:
                    try:
                        check_for_duplicates[ref_characteristics_to_object[cha]] = cha
                    except KeyError:
                        simlog.error(
                            'The characteristic\'s object was not found in the species supplied to the simulator \n'
                            'Perhaps a species is missing ?')

        for product in reaction.products:

            check_for_duplicates = {}
            for cha in product['characteristics']:

                # Ok I love try catches for checking if something is in a dict. Don't judge me
                try:
                    check_for_duplicates[ref_characteristics_to_object[cha]]
                    simlog.error('Illegal reaction, there is a product with multiple '
                                 'characteristics from the same Species referenced \n'
                                 'Please divide the reaction accordingly \n'
                                 f'Error in {ref_characteristics_to_object[cha]} both {cha} and '
                                 f'{check_for_duplicates[ref_characteristics_to_object[cha]]}'
                                 f' referenced at the same time \n'
                                 f'The characteristics: {ref_characteristics_to_object[cha].get_characteristics()}')
                except KeyError:
                    check_for_duplicates[ref_characteristics_to_object[cha]] = cha


def construct_reactant_structures(reactant_species, species_string_dict):
    """
        This function finds the corresponding strings according to the reactant meta-species with or without a query
        pack then in a list and return

        Parameters:
            reactant_species (meta-species object) = species objects of the involved species
            species_string_dict (dict)=  dictionary with species objects as keys and corresponding species strings
    """
    species_string_combinations = []

    for reactant in reactant_species:
        species_string_combinations.append(extract_species_strings(reactant['object'],
                                                                   reactant['characteristics'], species_string_dict))

    return species_string_combinations


def construct_order_structure(species_order_list, current_species_string_list):
    """
        Order structure for reaction order operations. Returns the cyclic_dictionary to be used by the order operator.
        The meta-species objects are the keys of this dictionary and a lists of species strings currently being used
        in the reaction are the values - Allowing the product to find it's corresponding species-string in a future
        step

        Parameters:
            species_order_list (list of meta-species objects): list of meta-species objects as they appear in the
            meta-reaction
            current_species_string_list (list of strings): list of strings in MobsPy format of the species currently
            in this specific reaction

        Returns:
            cyclic_dict (dict) = Dictionary where the keys are meta-species objects and the values are lists of species
    """
    cyclic_dict = {}
    for species_object, species_string in zip(species_order_list, current_species_string_list):

        try:
            cyclic_dict[species_object].append(species_string)
        except KeyError:
            cyclic_dict[species_object] = [species_string]

    return cyclic_dict


def construct_product_structure(reaction):
    '''
        This function unpacks the products in a meta-reaction

        Parameters:
            reaction = meta-reaction currently being analysed

        Returns:
            product_list: A list of dictionaries for each product with the meta-species object, the label and the
            characteristics
    '''
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

        Parameters:
            reactant_species_string_list (list of strings) = list of reactants in MobsPy format
            product_species_string_list (list of strings) = list of products in MobsPy format
            reaction_rate (str) = reaction rate expression as a string

        Returns:
            to_return (dict) = dictionary that packs the reactants products and rate
    """
    to_return = {'re': [], 'pr': [], 'kin': reaction_rate}

    reactant_count_dict = count_string_dictionary(reactant_species_string_list)
    product_count_dict = count_string_dictionary(product_species_string_list)

    for key in reactant_count_dict:
        to_return['re'].append((reactant_count_dict[key], key))

    for key in product_count_dict:
        to_return['pr'].append((product_count_dict[key], key))

    return to_return


def count_string_dictionary(list_of_strings):
    """
        Count the number of instances in a list and return them in a dictionary where the keys are the strings and
        the value the number of times it appeared in the list
    """
    to_return = {}

    for e in list_of_strings:
        try:
            to_return[e] += 1
        except KeyError:
            to_return[e] = 1

    return to_return


def extract_species_strings(species, characteristics, species_string_dict):
    """
        Extract the species strings from the species_string_dict based on the meta-species and characteristics given

        Parameters:
            species (meta-species object) = Instance of the meta-species to obtain the species strings for
            characteristics (str) = Characteristics to filter through
            species_string_dict (dict) = Dictionary where the meta-species are the keys and the values all their
            species

        Returns:
            species_strings_list (list of str) = species strings filtered from the meta-species
    """
    species_strings_list = []
    species_strings_to_filter = set()

    species_strings_to_filter = species_strings_to_filter.union(species_string_dict[species])

    for species_string in species_strings_to_filter:
        species_string_split = species_string.split('_dot_')[1:]
        if all([char in species_string_split for char in characteristics]):
            species_strings_list.append(species_string)

    return species_strings_list


def get_involved_species(reaction, species_string_dict):
    """
        This extracts all the involved meta-species inside a reaction
        This function is responsible for implementing the inheritance mechanism, by finding within each meta-species
        references set, if they reference the meta-species in the reaction

        Parameters:
            reaction (meta-reaction object)
            species_string_dict (dict) = dictionary where the keys are meta-species and the values are
            lists of species they contain (used only for the keys here)

        Returns:
            base_species_order (list of meta-species objects): order that the meta-species appear in the meta-reaction
            reactant_species_combination_list (list of lists of meta-species): list of lists of all meta-species that
            have inherited from the meta-species in the meta-reaction
    """
    reactant_species_combination_list = []
    base_species_order = []

    for reactant in reaction.reactants:
        flag_absent_reactant = False
        for _ in range(reactant['stoichiometry']):

            species_for_reactant = []
            base_species_order.append((reactant['object'], reactant['label']))

            for species in species_string_dict:
                if reactant['object'] in species.get_references():
                    species_for_reactant.append({'object': species,
                                                 'characteristics': reactant['characteristics']})
                    flag_absent_reactant = True

            if not flag_absent_reactant:
                simlog.error(f'Species {reactant["object"]} was not found in model \n'
                             f'For reaction {reaction} \n'
                             f'Please add the species or remove the reaction')

            reactant_species_combination_list.append(species_for_reactant)

    return base_species_order, reactant_species_combination_list


def create_all_reactions(reactions, species_string_dict,
                         ref_characteristics_to_object,
                         type_of_model, dimension):
    """
        This function creates all reactions
        Returns the reactions_for_sbml and parameters_for_sbml dictionary
        Those will be used by another module to create the SBML file

        Parameters:
            reactions (meta-reaction objects) = reactions objects constructed by the meta_class module
            species_string_dict (dict) = Meta-species object keys and respective species strings values
            ref_characteristics_to_object (dict) =  Characteristics as keys objects as values
            type_of_model (str) = stochastic or deterministic
            dimension (int) = model dimension 1D, 2D, 3D, .....

        Returns:
            reactions_for_sbml (dict) = dictionary with all reactions that will be added to the sbml model file
            parameters_for_sbml (dict) = parameters for the sbml model file
    """

    reactions_for_sbml = {}
    parameters_for_sbml = {}

    check_for_invalid_reactions(reactions, ref_characteristics_to_object)

    for reaction in reactions:

        base_species_order, reactant_species_combination_list = get_involved_species(reaction, species_string_dict)

        for combination_of_reactant_species in iterator_for_combinations(reactant_species_combination_list):

            reactant_species_string_combination_list = \
                construct_reactant_structures(combination_of_reactant_species, species_string_dict)

            for reactant_string_list in iterator_for_combinations(reactant_species_string_combination_list):

                product_object_list = construct_product_structure(reaction)
                order_structure = construct_order_structure(base_species_order, reactant_string_list)

                product_species_species_string_combination_list = reaction.order(order_structure, product_object_list,
                                                                                 species_string_dict,
                                                                                 ref_characteristics_to_object)

                for product_string_list in iterator_for_combinations(product_species_species_string_combination_list):
                    rate_string = fr.extract_reaction_rate(combination_of_reactant_species,
                                                           reactant_string_list
                                                           , reaction.rate, type_of_model,
                                                           dimension)

                    reactions_for_sbml['reaction_' + str(len(reactions_for_sbml))] = \
                        construct_single_reaction_for_sbml(reactant_string_list,
                                                           product_string_list,
                                                           rate_string)

    return reactions_for_sbml, parameters_for_sbml


if __name__ == '__main__':
    pass
