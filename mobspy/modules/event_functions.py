import mobspy.modules.meta_class_utils as mcu
import mobspy.modules.unit_handler as uh
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.modules.species_string_generator as ssg
from pint import Quantity


def format_event_dictionary_for_sbml(species_for_sbml, event_list, characteristics_to_object,
                                     volume, dimension, meta_species_to_simulate):
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

    # Convert count from triggers
    for ev in event_list:
        if ev['trigger'] != 'true':
            for i, e in enumerate(ev['trigger'].operation):
                if isinstance(e, Quantity):
                    ev['trigger'].operation[i] = uh.convert_counts(e, volume, dimension)

    for ev in event_list:
        if not ev['event_counts']:
            continue
        event_dictionary = {}
        for ec in ev['event_counts']:
            dummy = ssg.construct_species_char_list(ec['species'], ec['characteristics'],
                                                    characteristics_to_object, symbol='_dot_')
            event_dictionary[dummy] = uh.convert_counts(ec['quantity'], volume, dimension)
        if type(ev['trigger']) == str:
            reformed_event_list.append({'event_time': ev['event_time'], 'event_counts': event_dictionary,
                                        'trigger': ev['trigger']})
        else:
            for e in ev['trigger'].operation:
                if type(e) == dict:
                    if e['object'] not in meta_species_to_simulate:
                        simlog.error(f'Meta species {e["object"]} was used in an event but is not in the model')
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
            else:
                simlog.error(f'Species {key} was not found in the model')

        events_for_sbml['e' + str(i)] = {'trigger': event['trigger'],
                                         'delay': str(event['event_time']),
                                         'assignments': assignments}

    return events_for_sbml

