"""
    This model stores function used by the meta_class.py module
"""
import mobspy.simulation_logging.log_scripts as simlog

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


def combine_references(species1, species2):
    """
        Combine the sets of references of two species

        :param species1: (Meta-species object)
        :param species2: (Meta-species object)
    """
    return species1.get_references().union(species2.get_references())


def check_orthogonality_between_references(references):
    """
        Check if meta-species objects inside a reference do not have characteristics in common
        The sets of characteristics directly added to species must be independent

        :param references: (set) set of meta-species objects to check for independence
        :raise simlog.error: raises error if there are repeated characteristics in different meta-species
    """
    for i, reference1 in enumerate(references):
        for j, reference2 in enumerate(references):

            if i == j:
                continue

            if len(reference1.get_characteristics().intersection(reference2.get_characteristics())) != 0:
                simlog.error(f'The same characteristic can only be shared through inheritance. ' +
                             f'There are two characteristics directly added to two meta-species \n'
                             f'Repetition in: {reference1}, {reference2}'
                             f'Characteristics: {reference1.get_characteristics()}, {reference2.get_characteristics()}')


def complete_characteristics_with_first_values(spe_object, characteristics, characteristics_to_object):
    """
        This creates a string with the species object name and the set of all first characteristics added
        to it's base species used to construct it
        It allows us to search for the proper string in the model using the result of this function
        It's used to set the quantities using the call method in Species and Reacting_Species

        :param spe_object: (meta-species object) meta-species to generate the species states with the first values
        :param characteristics: (str) if there are characteristics to use instead of the default
        :param characteristics_to_object: (dict) dictionary with characteristics as keys and the meta-species which
        they have been added to directly as value (no inheritance or New)

        :return: Set with the species name and characteristics
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


def unite_characteristics(species):
    """
        This function unites the characteristics of all the given species

        :param species: (list of species or List_Species object)
    """
    characteristics = set()

    if species is not None:
        for spe in species:
            characteristics = characteristics.union(spe.get_characteristics())

    return characteristics


def create_orthogonal_vector_structure(species):
    """
        This creates the independent state-structure for the model
        It is just a dictionary where the keys are characteristics and the values are meta-species objects that have
        been directly added to that object (no inheritance or New)
        It simplifies the code by allowing to easily keep track of the 'axis' of each characteristic. Allowing
        for easy transformation on the products and others

        :param species: (meta-species objects) - meta-species objects used in a model

        :return ref_characteristics_to_object: (dict) a dictionary where the keys are characteristics and the
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
                    simlog.error(f'The same characteristic can only be shared through inheritance. ' +
                                 f'There are two characteristics directly added to two meta-species \n'
                                 f'Repetition in: {spe}, {ref_characteristics_to_object[cha]} \n'
                                 f'Characteristics: {spe.get_characteristics()}, '
                                 f'{ref_characteristics_to_object[cha].get_characteristics()} \n')

    return ref_characteristics_to_object


if __name__ == '__main__':
    pass
