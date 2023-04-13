import mobspy.simulation_logging.log_scripts as simlog
import inspect
from mobspy.modules.meta_class import List_Species
from pint import Quantity


def set_counts(count_dic):
    """
        Adds counts to meta-species using a given dictionary. Keys from this dictionary can be either meta-species
        objects or strings. Items must be the assign counts to species

        :param count_dic: Dictionary where keys can be either meta-species or strings
        :raise simlog.error: - If a count is assigned to a Reacting_Species with more than one meta-species.
        If there are two species with the same name and a string assignment is performed.
        If the keys are not strings or meta-species.
        If the counts are not Quantities, Floats or Ints.
        If the species was not found in the stack.
        :return: List_Species object. All meta-species that had a count assigned in this dictionary will be returned
        as a List_Species which can be passed as a model to the simulation object
    """
    for key, item in count_dic.items():
        if type(item) == int or type(item) == float or isinstance(item, Quantity):
            continue
        else:
            simlog.error(f'Reactant_species count assignment does not support the type {type(item)}',
                         stack_index=2)

    def find_species():
        found_species = set()

        for i in range(len(inspect.stack())):
            local_names = inspect.stack()[i][0].f_locals
            global_names = inspect.stack()[i][0].f_globals
            for key, item in global_names.items():
                try:
                    if item.is_species() and type(item) != type:
                        found_species.add(item)
                except AttributeError:
                    pass
            for key, item in local_names.items():
                try:
                    if item.is_species() and type(item) != type:
                        found_species.add(item)
                except AttributeError:
                    pass

        return found_species

    for key in count_dic:
        if type(key) == str:
            all_found_species = find_species()

    model = set()
    for key, item in count_dic.items():
        if type(key) == str:
            already_found = False
            str_name = key.split('.')[0]
            str_characteristics = set(key.split('.')[1:])
            for spe in all_found_species:
                if spe.get_name() == str_name and not already_found:
                    already_found = True
                    temp_set = set(str_characteristics)
                    temp_set.discard('all$')
                    if temp_set.issubset(spe.get_all_characteristics()):
                        spe.add_quantities(str_characteristics, item)
                    else:
                        simlog.error('Characteristics not found in species with equal name', stack_index=2)
                    model.add(spe)
                elif spe.get_name() == str_name and already_found:
                    simlog.error(f'There are two different meta-species with the same name. Set_counts cannot resolve',
                                 stack_index=2)
            if not already_found:
                simlog.error(f'Meta-species with the following name {key} not found', stack_index=2)
        else:
            try:
                if key.is_spe_or_reac():
                    if not key.is_species():
                        if len(key.list_of_reactants) != 1:
                            simlog.error('Assignment used incorrectly. Only one species at a time', stack_index=2)
                        model.add(key.list_of_reactants[0]['object'])
                    if key.is_species():
                        model.add(key)
                    key(item)
            except AttributeError:
                simlog.error('Keys must be either meta-species or strings',
                             stack_index=2)

    return List_Species(model)
