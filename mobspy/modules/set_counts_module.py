import mobspy.simulation_logging.log_scripts as simlog
import inspect
from mobspy.modules.meta_class import List_Species


def set_counts(count_dic):
    """
        Adds counts to meta-species using a given dictionary. Keys from this dictionary can be either meta-species
        objects or strings. Items must be the assign counts to species

        :return: List_Species object. All meta-species that had a count assigned in this dictionary will be returned
        as a List_Species which can be passed as a model to the simulation object
    """
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
                        simlog.error(f'At: set_counts with key {key} \n'
                                     + 'Characteristics not found in species with equal name')
                    model.add(spe)
                elif spe.get_name() == str_name and already_found:
                    simlog.error(f'At: set_counts \n'
                                 f'There are two different meta-species with the same name. Set_counts cannot resolve')
        else:
            try:
                if key.is_spe_or_reac():
                    key(item)
                    if key.is_species():
                        model.add(key)
                    else:
                        model.add(key.list_of_reactants[0]['object'])
            except AttributeError:
                simlog.error('Keys must be Meta-Species objects and items the value to assign to them')

    return List_Species(model)
