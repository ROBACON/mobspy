from itertools import product as itertools_product

def characteristics_dictionary(characteristics, characteristics_to_object):
    """
        This function constructs a dictionary that leads from the species string characteristics to their respective
        object. This structure allows to easily find where the characteristics in a species string are locates

        :param characteristics: (set) of characteristics
        :param characteristics_to_object: (dict) with characteristics as keys and their respective base meta-species
        object as value
    """
    object_to_characteristic = {}
    for characteristic in characteristics:
        if '$' not in characteristic:
            object_to_characteristic[characteristics_to_object[characteristic]] = characteristic
    return object_to_characteristic


def construct_species_char_list(spe_object, characteristics, characteristics_to_object, symbol=None):
    """
        This function constructs a list in the format ['species_name', 'char1', 'char2', ...]. It generetes this list
        for a given meta-species and the specified characteristics. Values of characteristics not specified in a
        particular position are replaced by their default value. If a symbol is given it generates a string from
        the list using the symbol to join it.

        :param spe_object: meta-species object to be used
        :param characteristics: (set) of characteristics given
        :param characteristics_to_object: (dict) with characteristics as keys and their respective base meta-species
        object as value
        :param symbol: (str) usually . or _dot_, connects the elements from the list using the specified symbol
    """
    if characteristics == 'std$':
        characteristics = set()

    ordered_references_list = spe_object.get_ordered_references()

    objects_to_characteristic = characteristics_dictionary(characteristics, characteristics_to_object)

    species_char_list = [spe_object]
    for obj in ordered_references_list:
        if obj in objects_to_characteristic:
            species_char_list.append(objects_to_characteristic[obj])
        else:
            species_char_list.append(obj.first_characteristic)

    if symbol is not None:
        species_char_list = symbol.join([spe_object.get_name()] + species_char_list [1:]) \
            if len(species_char_list) > 1 else spe_object.get_name()

    return species_char_list


def construct_all_combinations(spe_object, characteristics, characteristics_to_object, symbol=None):
    """
       This function constructs all possible list in the format ['species_name', 'char1', 'char2', ...] using
       all combinations of characteristics from the vector coordinates not used in the characteristics specified
       in the function argument. If a symbol is given it generates a string from the list using the symbol to join it.

       :param spe_object: meta-species object to be used
       :param characteristics: (set) of characteristics given
       :param characteristics_to_object: (dict) with characteristics as keys and their respective base meta-species
       object as value
       :param symbol: (str) usually . or _dot_, connects the elements from the list using the specified symbol
    """

    if characteristics == 'std$':
        characteristics = set()

    spe_object.order_references()
    ordered_references_list = spe_object.get_ordered_references()

    objects_to_characteristic = characteristics_dictionary(characteristics, characteristics_to_object)

    list_of_all_possibilities = [[spe_object]]
    for obj in ordered_references_list:
        if obj in objects_to_characteristic:
            list_of_all_possibilities.append([objects_to_characteristic[obj]])
        else:
            list_of_all_possibilities.append(list(obj.get_characteristics()))

    to_return = []
    for i in itertools_product(*list_of_all_possibilities):
        if symbol is not None:
            to_return.append(symbol.join([spe_object.get_name()] + list(i)[1:])
                             if len(i) > 1 else spe_object.get_name())
        else:
            to_return.append(list(i))

    return to_return



