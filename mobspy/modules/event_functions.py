import mobspy.modules.meta_class_utils as mcu
import mobspy.modules.unit_handler as uh
import mobspy.simulation_logging.log_scripts as simlog


def prepare_event_count_dictionary(meta_species_counts):
    def get_count_from_meta_species(meta_spe):
        if meta_spe not in already_counted:
            for count in meta_spe.get_quantities():
                if count['quantity'] != 0:
                    event_count_list.append({'species': meta_spe,
                                             'characteristics': count['characteristics'],
                                             'quantity': count['quantity']})
            already_counted.add(meta_spe)

    event_count_list = []
    already_counted = set()
    for meta_species in meta_species_counts:

        # Test for species or reacting_species
        is_reacting_species = False
        try:
            meta_species.get_quantities()
        except:
            is_reacting_species = True

        if not is_reacting_species:
            get_count_from_meta_species(meta_species)
        else:
            get_count_from_meta_species(meta_species.list_of_reactants[0]['object'])

    return event_count_list


def extract_results_dictionary(data):
    results_dictionary = {}
    for species in data:
        if species == 'Time':
            continue
        results_dictionary[str(species)] = [run[-1] for run in data[species]['runs']]


def correct_species_count(simulation_dictionary, event_dictionary, results_dictionary, parameters):
    # Event Dictionary is not received in the correct units when it arrives to this point
    # One needs to call the unit conversion function to account for this
    # Not used yet
    for key in event_dictionary:
        event_dictionary[key] = uh.convert_counts(event_dictionary[key], parameters['volume'], parameters['dimension'])
    pass


def join_event_data(total_packed_data, data):
    def sum_list_with_element(e, list):
        to_return_list = []
        for l in list:
            to_return_list.append(e + l)
        return to_return_list

    species_already_added = []
    for key in data:
        if key == 'Time':
            continue
        else:
            species_already_added.append(key)
            if key in total_packed_data:
                for i in len(total_packed_data[key]['runs']):
                    total_packed_data[key]['runs'][i] = total_packed_data[key]['runs'][i] + data[key]['runs'][i]
            else:
                total_packed_data[key] = {'runs': []}
                for run in data[key]['runs']:

                    if 'Time' in total_packed_data.keys():
                        dummy_list = [0 for _ in range(len(total_packed_data['Time']))] + run
                    else:
                        dummy_list = run

                    total_packed_data[key]['runs'].append(dummy_list)

    for key in total_packed_data:
        if key == 'Time':
            continue
        else:
            if key not in species_already_added:
                for i, run in enumerate(total_packed_data[key]['runs']):
                    dummy_list = [total_packed_data[key]['runs'][-1] for _ in range(len(data['Time']))]
                    total_packed_data[key]['runs'][i] = total_packed_data[key]['runs'][i] + dummy_list

    if 'Time' in total_packed_data:
        total_packed_data['Time'] = total_packed_data['Time'] \
                                    + sum_list_with_element(total_packed_data['Time'][-1], data['time'])
    else:
        total_packed_data['Time'] = data['Time']

    return total_packed_data


def format_event_dictionary_for_sbml(species_for_sbml, event_list, characteristics_to_object,
                                     volume, dimension):
    """
    Creates events_for_sbml dictionary for sbml file construction on the sbml_simulator/SBMLWriter.py

    Parameters:
        species_for_sbml = {'species_string': count, ....}
        event_list = [{'species':'meta_species_object', 'characteristics':['list of characteristics'],
            'quantity': number or pint object}
        characteristics_to_object = {'characteristic':'meta_species_object', .....}
        volume = Simulation volume for unit conversion
        dimension = Dimension of the system 2D, 3D, 4D, etc

    Returns:
        event dictionary for the sbml file construction
        events = {'e': \
              { 'trigger': 'true', \
                'delay': '10', \
                'assignments': [('M','1'),], \
              }, \
            }
    """

    reformed_event_list = []
    for ev in event_list:
        if not ev['event_counts']:
            continue
        event_dictionary = {}
        for ec in ev['event_counts']:
            dummy = mcu.complete_characteristics_with_first_values(ec['species'], ec['characteristics'],
                                                                   characteristics_to_object)
            dummy = '_dot_'.join(dummy)
            event_dictionary[dummy] = uh.convert_counts(ec['quantity'], volume, dimension)
        if type(ev['trigger']) == str:
            reformed_event_list.append({'event_time': ev['event_time'], 'event_counts': event_dictionary,
                                        'trigger': ev['trigger']})
        else:
            reformed_event_list.append({'event_time': ev['event_time'], 'event_counts': event_dictionary,
                                        'trigger': ev['trigger'].generate_string_from_vec_space(species_for_sbml)})

    events_for_sbml = {}
    for i, event in enumerate(reformed_event_list):
        assignments = []
        for key_1 in event['event_counts']:
            dummy_set = set(key_1.split('_dot_'))
            for key_2 in species_for_sbml:
                if dummy_set == set(key_2.split('_dot_')):
                    assignments.append((key_2, str(event['event_counts'][key_1])))
                    break
            else:
                simlog.error(f'Species {dummy_set} was not found in the model')

        events_for_sbml['e' + str(i)] = {'trigger': event['trigger'],
                                         'delay': str(event['event_time']),
                                         'assignments': assignments}

    return events_for_sbml


def inspect_event_triggers(self, list_of_logic_resolver_objects):
    pass