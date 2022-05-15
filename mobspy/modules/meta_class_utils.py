"""
    This model stores function used by the meta_class.py module
"""
import copy
import itertools
import mobspy.simulation_logging.log_scripts as simlog


def unite_dictionaries(first_dict, second_dict):
    """
        Performs the union of two dictionaries whose values are sets.

        Parameters:
            first_dict (dict)
            second_dict (dict)
    """
    for key in second_dict:
        try:
            first_dict[key].union(second_dict[key])
        except KeyError:
            first_dict[key] = second_dict[key]


def combine_references(species1, species2):
    """
        Combine the sets of references of two species

        Parameters:
            species1 (Meta-species object)
            species2 (Meta-species object)
    """
    #if species1.get_references().intersection(species2.get_references()):
    #    simlog.warning(f'A product was executed between species with a common meta-species\n')

    return species1.get_references().union(species2.get_references())


def check_orthogonality_between_references(references):
    """
        Check if meta-species objects inside a reference do not have characteristics in common
        The sets of characteristics directly added to species must be independent

        Parameter:
            references (set) =  set of meta-species objects to check for independence
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

        Parameters:
            spe_object (meta-species object) = meta-species to generate the species states with the first values
            characteristics (str) = if there are characteristics to use instead of the default
            characteristics_to_object (dict) = dictionary with characteristics as keys and the meta-species which
            they have been added to directly as value (no inheritance or New)

        Returns:
            Set with the species name and characteristics
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
        Species are named for the SBML as species_name_dot_characteristic1_dot_characteristic2
        So this transforms them into a set

        Parameters:
            species_string (str) = species string in MobsPy for SBML format (with _dot_ instead of .)
    """
    return set(species_string.split('_dot_'))


def turn_set_into_0_value_dict(set):
    to_return = {}
    for e in set:
        to_return[e] = 0
    return to_return


def unite_characteristics(species):
    """
        This function unites the characteristics of all the given species

        Parameters:
            species (list of species or ParallelSpecies object)
    """
    characteristics = set()

    if species is not None:
        for spe in species:
            characteristics = characteristics.union(spe.get_characteristics())

    return characteristics


def extract_characteristics(spe):
    """
        Extracts all the characteristics from a species and it's references

        Parameters:
            spe (meta-species object)
    """
    lists_of_characteristics = []
    for reference in spe.get_references():
        lists_of_characteristics.append(reference.get_characteristics())

    return lists_of_characteristics


def create_orthogonal_vector_structure(species):
    """
        This creates the independent state-structure for the model
        It is just a dictionary where the keys are characteristics and the values are meta-species objects that have
        been directly added to that object (no inheritance or New)
        It simplifies the code by allowing to easily keep track of the 'axis' of each characteristic. Allowing
        for easy transformation on the products and others

        Parameters:
            species (meta-species objects) - meta-species objects used in a model

        Returns:
            ref_characteristics_to_object (dict) = a dictionary where the keys are characteristics and the
            values are meta-species objects that have been directly added to that object
    """
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
        This function combines the species name with all the characteristics of it's references.
        This way it creates the strings that will be used for the SBML
        All possible combinations are created, with no intersections between characteristics of
        a same referenced species

        Parameters:
            spe_object (meta-species object) = received species object
            sets_of_characteristics (set of sets (thus not a set)) = all sets of characteristics from the
            referenced objects

        Returns:
            set_of_species (set) = set of species from the meta-species object
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
