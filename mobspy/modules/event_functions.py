import mobspy.modules.meta_class_utils as mcu
import mobspy.modules.unit_handler as uh
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.modules.species_string_generator as ssg
from pint import Quantity


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


def format_event_dictionary_for_sbml(species_for_sbml, event_list, characteristics_to_object,
                                     volume, dimension, meta_species_to_simulate):
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

