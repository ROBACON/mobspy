import itertools


def characteristics_dictionary(characteristics, characteristics_to_object):
    object_to_characteristic = {}
    for characteristic in characteristics:
        object_to_characteristic[characteristics_to_object[characteristic]] = characteristic
    return object_to_characteristic


def construct_species_char_list(spe_object, characteristics, characteristics_to_object, symbol=None):

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
            if len(species_char_list ) > 1 else spe_object.get_name()

    return species_char_list


def construct_all_combinations(spe_object, characteristics, characteristics_to_object, symbol=None):

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
    for i in itertools.product(*list_of_all_possibilities):
        if symbol is not None:
            to_return.append(symbol.join([spe_object.get_name()] + list(i)[1:])
                             if len(i) > 1 else spe_object.get_name())
        else:
            to_return.append(list(i))

    return to_return



