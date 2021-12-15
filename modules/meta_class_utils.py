import copy


def unite_dictionaries(first_dict, second_dict):
    for key in second_dict:
        try:
            first_dict[key].union(second_dict[key])
        except KeyError:
            first_dict[key] = second_dict[key]


def combine_references(species1, species2):
    return species1.species_references.union(species2.species_references)


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
            characteristics = characteristics.union(spe.characteristics)

    return characteristics


def extract_characteristics(spe):

    lists_of_characteristics = []
    for reference in spe.species_references:
        lists_of_characteristics.append(reference.characteristics)

    return lists_of_characteristics


def combine_characteristics(spe_object, lists_of_definitions):

    # Remove empty sets from the list
    lists_of_definitions = [i for i in lists_of_definitions if i != set()]

    # Defining needed structures
    set_of_species = set()
    accumulated_list = [spe_object.get_name()]
    species_from_characteristic = {}

    # Recursive exponential combination implementation
    def recursive_combine_properties(i, list_definitions, accumulated_list):

        for j, characteristic in enumerate(list_definitions):

            if j == 0:
                accumulated_list.append(characteristic)
            else:
                accumulated_list[-1] = characteristic

            try:
                recursive_combine_properties(i + 1, lists_of_definitions[i + 1], copy.deepcopy(accumulated_list))
            except IndexError:
                spe = '.'.join(accumulated_list)
                set_of_species.add(spe)

                for e in accumulated_list:
                    # Skip the name of the species
                    if e == spe_object.get_name():
                        continue
                    else:
                        try:
                            species_from_characteristic[e].add(spe)
                        except KeyError:
                            species_from_characteristic[e] = {spe}

    recursive_combine_properties(0, lists_of_definitions[0], accumulated_list)

    species_from_object = {}
    for spe_reference in spe_object.species_references:
        species_from_object[spe_reference] = set_of_species

    return set_of_species, species_from_characteristic, species_from_object



if __name__ == '__main__':

    print(combine_characteristics('Ecoli', [{'young', 'old'},{'green', 'red'}, {'deja', 'vu'}]))

