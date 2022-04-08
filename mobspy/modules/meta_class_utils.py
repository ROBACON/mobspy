import copy
import itertools
import mobspy.simulation_logging.log_scripts as simlog


def unite_dictionaries(first_dict, second_dict):
    """
        first_dict = first dictionary
        second_dict = second dictionary
        Just fuse both dictionaries based on the keys. First key value takes precedence
    """
    for key in second_dict:
        try:
            first_dict[key].union(second_dict[key])
        except KeyError:
            first_dict[key] = second_dict[key]


def combine_references(species1, species2):
    """
        Combine the sets of references of two species
    """
    #if species1.get_references().intersection(species2.get_references()):
    #    simlog.warning(f'A product was executed between species with a common meta-species\n')

    return species1.get_references().union(species2.get_references())


def check_orthogonality_between_references(references):
    """
        Check if base species do not have characteristics in common
        The sets of characteristics must be independent for BASE SPECIES
    """
    for i, reference1 in enumerate(references):
        for j, reference2 in enumerate(references):

            if i == j:
                continue

            if len(reference1.get_characteristics().intersection(reference2.get_characteristics())) != 0:
                simlog.error(f'A characteristic must be unique for species construction'
                             f'Repetition in: {reference1}, {reference2}'
                             f'Characteristics: {reference1.get_characteristics()}, {reference2.get_characteristics()}')


def complete_characteristics_with_first_values(spe_object, characteristics, characteristics_to_object):
    """
        This creates a string with the species object name and the set of all first characteristics added
        to it's base species used to construct it
        It allows us to search for the proper string in the model using the result of this function
        It's used to set the quantities using the call method in Species and Reacting_Species
    """
    if characteristics == 'std$':
        characteristics = set()

    vector_elements = {}
    for cha in characteristics:
        vector = characteristics_to_object[cha]
        if vector in vector_elements:
            simlog.error('The assignment refers to multiple strings')
        else:
            vector_elements[vector] = True

    first_characteristics = set()
    for reference in spe_object.get_references() - set(vector_elements.keys()):
        if reference.first_characteristic:
            first_characteristics.add(reference.first_characteristic)

    return {spe_object.get_name()}.union(first_characteristics).union(characteristics)


def extract_characteristics_from_string(species_string):
    """
        Species are named for the SBML as species_name.characteristic1.characteristic2
        So this transforms them into a set
    """
    return set(species_string.split('_dot_'))


def turn_set_into_0_value_dict(set):
    to_return = {}
    for e in set:
        to_return[e] = 0
    return to_return


def unite_characteristics(species):
    '''
        This function checks if there are repeated characteristics in the model
    :param species: All species added to the model
    ;return: characteristics
    '''
    characteristics = set()

    if species is not None:
        for spe in species:
            characteristics = characteristics.union(spe.get_characteristics())

    return characteristics


def extract_characteristics(spe):

    lists_of_characteristics = []
    for reference in spe.get_references():
        lists_of_characteristics.append(reference.get_characteristics())

    return lists_of_characteristics


def create_orthogonal_vector_structure(species):
    '''
    :param species: All species used in the model
    :return: Two hashes one that references objects to characteristics
            and other that references characteristics to objects

            Here we create the basis of the orthogonal vector structure
            Each base property object must contain a set of independent characteristics
            We think of it as properties being unity vectors in a cartesian coordinate system
            And the characteristics as values

            We use the dictionary structure to easily define the reactions as being transformations
            within the same unity vector
    '''
    ref_characteristics_to_object = {}
    for spe in species:
        for prop in spe.get_references():
            for cha in prop.get_characteristics():

                if cha not in ref_characteristics_to_object:
                    ref_characteristics_to_object[cha] = prop
                elif ref_characteristics_to_object[cha] == prop:
                    pass
                else:
                    simlog.error(f'A characteristic must be unique for each species \n'
                                f'Repetition in: {spe}, {ref_characteristics_to_object[cha] } \n'
                                f'Characteristics: {spe.get_characteristics()}, {ref_characteristics_to_object[cha].get_characteristics()} \n')

    return ref_characteristics_to_object


def create_species_strings(spe_object, sets_of_characteristics):
    """
        spe_object : received species object
        sets_off_characteristics : all sets of characteristics from the referenced objects

        This function combines the species name with all the characteristics of it's references.
        This way it creates the strings that will be used for the SBML
        All possible combinations are created, with no intersections between characteristics of
        a same referenced species
    """

    # Remove empty sets from the list and transform sets in lists
    # Sort the lists so the characteristics will appear in the same order e
    lists_of_definitions = [sorted(list(i)) for i in sets_of_characteristics if i != set()]
    set_of_species = set()
    lists_of_definitions = sorted(lists_of_definitions)

    for i in itertools.product(*lists_of_definitions):
        if lists_of_definitions:
            set_of_species.add(spe_object.get_name() + '_dot_' + '_dot_'.join(i))
        else:
            set_of_species.add(spe_object.get_name())

    return set_of_species


if __name__ == '__main__':
    pass
