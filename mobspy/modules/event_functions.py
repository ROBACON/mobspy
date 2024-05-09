from mobspy.modules.unit_handler import convert_counts as uh_convert_counts
from mobspy.simulation_logging.log_scripts import error as simlog_error
from mobspy.modules.species_string_generator import construct_all_combinations as ssg_construct_all_combinations, \
     construct_species_char_list as ssg_construct_species_char_list
from pint import Quantity
from mobspy.modules.mobspy_parameters import Mobspy_Parameter as mp_Mobspy_Parameter
from mobspy.modules.function_rate_code import search_for_parameters_in_str as frc_search_for_parameters_in_str


def format_event_dictionary_for_sbml(species_for_sbml, event_list, characteristics_to_object,
                                     volume, dimension, meta_species_to_simulate,
                                     parameter_exist,parameters_in_events):
    """
        Creates events_for_sbml dictionary for sbml file construction on the sbml_simulator/SBMLWriter.py

        :param species_for_sbml: (dict) {'species_string': count, ....}
        :param event_list: [{'species':'meta_species_object', 'characteristics':['list of characteristics'],
            'quantity': number or pint object}
        :param characteristics_to_object: {'characteristic':'meta_species_object', .....}
        :param volume: Simulation volume for unit conversion
        :param dimension: Dimension of the system 2D, 3D, 4D, e
        :return: event dictionary for the sbml file construction
        :rtype: events = {'e': { 'trigger': 'true', 'delay': '10', 'assignments': [('M','1'),]}}
    """
    reformed_event_list = []
    species_in_events = set()

    # Convert count from triggers
    for ev in event_list:
        if ev['trigger'] != 'true':
            for i, e in enumerate(ev['trigger'].operation):
                if isinstance(e, Quantity):
                    ev['trigger'].operation[i] = uh_convert_counts(e, volume, dimension)

    for ev in event_list:
        if not ev['event_counts']:
            continue
        event_dictionary = {}

        # All assignments never take priority over specific assignments
        for ec in ev['event_counts']:
            if 'all$' not in ec['characteristics']:
                continue

            temp_char = set(ec['characteristics'])
            temp_char.remove('all$')
            dummy = ssg_construct_all_combinations(ec['species'], temp_char,
                                                   characteristics_to_object, symbol='_dot_')
            for d in dummy:
                if type(ec['quantity']) != str:
                    event_dictionary[d] = uh_convert_counts(ec['quantity'], volume, dimension)
                else:
                    event_dictionary[d] = ec['quantity']

        for ec in ev['event_counts']:
            if 'all$' in ec['characteristics']:
                continue

            dummy = ssg_construct_species_char_list(ec['species'], ec['characteristics'],
                                                    characteristics_to_object, symbol='_dot_')

            if type(ec['quantity']) != str:
                if isinstance(ec['quantity'], mp_Mobspy_Parameter):
                    parameters_in_events.add(ec['quantity'])
                    event_dictionary[dummy] = ec['quantity'].name

                event_dictionary[dummy] = uh_convert_counts(ec['quantity'], volume, dimension)
            else:
                if parameter_exist != {}:
                    frc_search_for_parameters_in_str(ec['quantity'], parameter_exist, parameters_in_events)
                event_dictionary[dummy] = ec['quantity']

        if type(ev['trigger']) == str:
            reformed_event_list.append({'event_time': ev['event_time'], 'event_counts': event_dictionary,
                                        'trigger': ev['trigger']})
        else:
            for e in ev['trigger'].operation:
                if type(e) == dict:
                    if e['object'] not in meta_species_to_simulate:
                        simlog_error(f'Meta species {e["object"]} was used in an event but is not in the model')
            reformed_event_list.append({'event_time': ev['event_time'],
                                        'event_counts': event_dictionary,
                                        'trigger': ev['trigger'].generate_string(characteristics_to_object,
                                                                                 to_sort=True)})

    events_for_sbml = {}
    for i, event in enumerate(reformed_event_list):
        assignments = []
        for key in event['event_counts']:
            if key in species_for_sbml:
                assignments.append((key, str(event['event_counts'][key])))
                species_in_events.add(key)
            else:
                simlog_error(f'Species {key} used in an event assignment but it is not in the model')

        assignments.sort()

        if event['event_time']:
            pass

        if isinstance(event['event_time'], mp_Mobspy_Parameter):
            for par in event['event_time']._parameter_set:
                parameters_in_events.add(par)

        events_for_sbml['e' + str(i)] = {'trigger': event['trigger'],
                                         'delay': str(event['event_time']),
                                         'assignments': assignments}

    return events_for_sbml, species_in_events
