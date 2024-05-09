from mobspy.simulation_logging.log_scripts import error as simlog_error, warning as simlog_warning
from mobspy.modules.mobspy_parameters import Mobspy_Parameter as mp_Mobspy_Parameter
from mobspy.modules.species_string_generator import construct_all_combinations as ssg_construct_all_combinations, \
    construct_species_char_list as ssg_construct_species_char_list
from pint import Quantity
from mobspy.modules.unit_handler import convert_counts as uh_convert_counts, \
    extract_length_dimension as uh_extract_length_dimension, convert_volume as uh_convert_volume
from mobspy.modules.meta_class import EndFlagSpecies
from mobspy.modules.reaction_construction_nb import create_all_reactions as rc_create_all_reactions
from mobspy.modules.event_functions import format_event_dictionary_for_sbml as eh_format_event_dictionary_for_sbml
from mobspy.modules.assignments_implementation import Assign as asgi_Assign, Asg as asgi_Asg
from copy import deepcopy


class Compiler:
    """
        The compiler is responsible for constructing the model from the meta-species and meta-reactions given
        and checking if it is a valid MobsPy model
    """

    # Basic storage of defined variables for Compiler
    entity_counter = 0
    ref_characteristics_to_object = {}
    last_rate = None

    @classmethod
    def add_to_parameters_to_sbml(cls, parameters_used, parameters_for_sbml, parameters_to_add):
        for parameter in parameters_to_add:
            if parameter.name in parameters_used:
                parameters_used[parameter.name]['used_in'].add('$sbml')
            else:
                temp = {'name': parameter.name, 'values': parameter.value, 'used_in': {'$sbml'}}
                parameters_used[parameter.name] = temp

            try:
                parameters_for_sbml[parameter.name] = (parameter.value[0], f'dimensionless')
            except:
                parameters_for_sbml[parameter.name] = (parameter.value, f'dimensionless')

    @classmethod
    def override_get_item(cls, object_to_return, item):
        """
            Due to priority in Python the item is stored before the reaction
            So it is stored in the Compiler level and passed to the reaction object in the end
            Important: the rate is used before the compilation level, it's used when the reaction has been completely
            defined

            :param object_to_return: (Species or Reacting) returns the object that the __getitem__ was performed on
            :param item: (int, float, callable, Quantity) stored reaction rate
        """
        cls.last_rate = item
        return object_to_return

    @classmethod
    def add_phantom_reactions(cls, reactions_for_sbml, species_to_add_phantom_reaction):
        """
            This method is here because basico cannot compile species that have been assigned to events but are not
            present in reactions.
            Interestingly enough this error does not happen in python 3.10.
            However, to be compatible with earlier versions we add an extremely slow phantom reaction for each one
            of those species with a coefficient of 1e-100. The idea is that the reaction is so slow that it should
            not make any sizable difference in any model. It adds them directly to the reactions for sbml dictionary

            :param reactions_for_sbml: reactions dictionary in python sbml format
            :param species_to_add_phantom_reaction: species that will be assigned a phantom reaction
        """
        for i, spe in enumerate(species_to_add_phantom_reaction):
            reactions_for_sbml['phantom_reaction_' + str(i)] = {'re': [(1, spe)], 'pr': [],
                                                                'kin': str(spe) + ' * 1e-100'}

    @classmethod
    def compile(cls, meta_species_to_simulate, reactions_set, species_counts, orthogonal_vector_structure,
                volume=1, dimension=None,
                type_of_model='deterministic', verbose=True, event_dictionary=None,
                continuous_sim=False, ending_condition=None, skip_expression_check=False):
        """
            Here we construct the species for Copasi using their objects

            It handles the construction of species_for_sbml, reactions_for_sbml, parameters_for_sbml,
            mappings_for_sbml, model_str which will be used to write the sbml_string for basiCO
            It also performs checks to see if the model is valid

            :param meta_species_to_simulate: (List_Species or Species) Species to generate the model
            :param reactions_set: (set) Set of meta-reactions collected when the simulation objected was constructed
            :param species_counts: (dict) All counts assigned to species in the model
             before the simulation object was constructed
            :param orthogonal_vector_structure: (dict) ref_characteristics_to_objects in other modules. Dictionary
            with the characteristics as keys and the objects as values. Characteritics pointing to their coordinate
            in the vector space
            :param volume: (int, flot) Simulation volume
            :param type_of_model: (str) deterministic or stochastic
            :param verbose: (bool) print or not the model results by generating a model_str
            meta-reaction ( or Default we have the round-robin) and how MobsPy will handle meta-species in
            the products that are not referenced in the reactants - see order_operators.py
            :param event_dictionary: (dict) dictionary with all the events to be added to the model
            :param continuous_sim: (bool) simulation with conditional duration or not
            :param ending_condition: (MetaSpeciesLogicResolver) trigger for the ending contion

            :returns: species_for_sbml (dict) = dictionary where the species strings (in MobsPy format) are the keys and
            the values are their counts inside the model. reactions_for_sbml (dict) = Reaction in dictionary format
            with keys 're', 'pr' and 'rate'. Values are strings parameters_for_sbml (dict) = value with parameter
            name and tuple with quantity and unit. For standard models it only holds the volume mappings_for_sbml
            (dict) = Dictionary that maps a meta-species to the set of species strings that belong to it. model_str
            (str) = str of the variables defined above in a user-readable format. events_for_sbml (dict) = dictionary
            containing the event trigger and count assignments. assigned_species (list) = list of meta-species which
            had counts assigned to them during this model execution (used when composing simulations)
        """
        # Check to see if all species are named
        # Parameter compilation as well
        names_used = set()

        # Removing repeated elements to ensure only one meta-species in the simulation
        meta_species_to_simulate = meta_species_to_simulate.remove_repeated_elements()
        black_listed_names = {'Time', 'Rev', 'All'}
        for i, species in enumerate(meta_species_to_simulate):
            if '_dot_' in species.get_name():
                simlog_error(f'In species: {species.get_name()} \n _dot_ cannot be used in meta-species names')
            if species.get_name() in black_listed_names:
                simlog_error(f'The name {species.get_name()} is not allowed for meta-species please change it')
            if '$' in species.get_name():
                simlog_error(f'In species: {species.get_name()} \n'
                             f'An error has occurred and one of the species was either not named or named with the '
                             f'restricted $ symbol')
            if species.get_name() in names_used:
                simlog_error(f'Names must be unique for all species\n' +
                             f'The repeated name is {species.get_name()} in position {i}\n' +
                             f'Another possibility could be a repeated meta-species in the model')
            names_used.add(species.get_name())

        # Define parameter dictionary and get parameter stack
        parameters_used = {}
        parameter_exist = {}
        if mp_Mobspy_Parameter.parameter_stack != {}:
            parameter_exist = mp_Mobspy_Parameter.parameter_stack

        # Order the species references for later usage
        for species in meta_species_to_simulate:
            species.order_references()

        # Start by creating the Mappings for the SBML
        # Convert to user friendly format as well
        mappings_for_sbml = {}
        for spe_object in meta_species_to_simulate:
            mappings_for_sbml[spe_object.get_name()] = []

        # List of Species objects
        species_for_sbml = {}
        mappings_for_sbml = {}
        for spe_object in meta_species_to_simulate:
            species_string_list = ssg_construct_all_combinations(spe_object, 'std$',
                                                                 orthogonal_vector_structure)
            for x in species_string_list:
                x[0] = x[0].get_name()
            mappings_for_sbml[spe_object.get_name()] = ['.'.join(x) for x in species_string_list]
            for species_string in species_string_list:
                species_for_sbml['_dot_'.join(species_string)] = 0

        # Default dimension equal to three after update
        if isinstance(volume, Quantity):
            dimension = uh_extract_length_dimension(str(volume.dimensionality), dimension)
        elif dimension is None:
            dimension = 3

        for count in species_counts:
            if isinstance(count['quantity'], Quantity):
                if uh_extract_length_dimension(str(count['quantity'].dimensionality), dimension):
                    dimension = uh_extract_length_dimension(str(count['quantity'].dimensionality), dimension)

        # Check volume:
        volume = uh_convert_volume(volume, dimension)

        # Add the flag species used for verifying if the simulation is over
        if continuous_sim:
            species_for_sbml[EndFlagSpecies.get_name()] = 0

        # Check if there are any mols in the units
        has_mole = False
        for count in species_counts:
            if isinstance(count['quantity'], Quantity):
                if 'substance' in str(count['quantity'].dimensionality):
                    has_mole = True

        # Assignments with all do not take priority
        # This allows specific assignments to override All assignments
        assigned_species = []
        parameters_in_counts = set()
        for count in species_counts:

            if 'all$' not in count['characteristics']:
                continue

            temp_set = set(count['characteristics'])
            temp_set.remove('all$')
            species_strings = ssg_construct_all_combinations(count['object'], temp_set,
                                                             orthogonal_vector_structure,
                                                             symbol='_dot_')

            if isinstance(count['quantity'], mp_Mobspy_Parameter):
                parameters_in_counts.add(count['quantity'])
                if count['quantity'].name in parameters_used:
                    parameters_used[count['quantity'].name]['used_in'] \
                        = parameters_used[count['quantity'].name]['used_in'].union(set(species_strings))
                else:
                    temp = {'name': count['quantity'].name, 'values': count['quantity'].value,
                            'used_in': set(species_strings)}
                    parameters_used[count['quantity'].name] = temp

            temp_count = uh_convert_counts(count['quantity'], volume, dimension)
            for spe_str in species_strings:
                if type(temp_count) == float and not type_of_model == 'deterministic':
                    simlog_warning('The stochastic simulation rounds floats to integers')
                    species_for_sbml[spe_str] = int(temp_count)
                    assigned_species.append(spe_str)
                else:
                    species_for_sbml[spe_str] = temp_count
                    assigned_species.append(spe_str)

        # Set initial counts for SBML
        # Create the list HERE
        for count in species_counts:

            if 'all$' in count['characteristics']:
                continue

            species_string = ssg_construct_species_char_list(count['object'], count['characteristics'],
                                                             orthogonal_vector_structure,
                                                             symbol='_dot_')

            if isinstance(count['quantity'], mp_Mobspy_Parameter):
                parameters_in_counts.add(count['quantity'])
                if count['quantity'].name in parameters_used:
                    parameters_used[count['quantity'].name]['used_in'].add(species_string)
                else:
                    temp = {'name': count['quantity'].name, 'values': count['quantity'].value,
                            'used_in': {species_string}}
                    parameters_used[count['quantity'].name] = temp

            temp_count = uh_convert_counts(count['quantity'], volume, dimension)
            if type(temp_count) == float and not type_of_model == 'deterministic':
                simlog_warning('The stochastic simulation rounds floats to integers')
                species_for_sbml[species_string] = int(temp_count)
                assigned_species.append(species_string)
            else:
                species_for_sbml[species_string] = temp_count
                assigned_species.append(species_string)

        # BaseSpecies reactions for SBML with theirs respective parameters and rates
        # What do I have so far
        # Species_String_Dict and a set of reaction objects in Reactions_Set
        parameters_in_reaction = set()
        reactions_for_sbml, parameters_in_reaction = rc_create_all_reactions(reactions_set,
                                                                             meta_species_to_simulate,
                                                                             orthogonal_vector_structure,
                                                                             type_of_model, dimension,
                                                                             parameter_exist, parameters_in_reaction,
                                                                             skip_expression_check)

        parameters_for_sbml = {'volume': (volume, f'dimensionless')}
        cls.add_to_parameters_to_sbml(parameters_used, parameters_for_sbml, parameters_in_reaction)

        # O(n^2) reaction check for doubles
        for i, r1 in enumerate(reactions_for_sbml):
            for j, r2 in enumerate(reactions_for_sbml):
                if i == j:
                    continue

                if reactions_for_sbml[r1]['re'] == reactions_for_sbml[r2]['re'] \
                        and reactions_for_sbml[r1]['pr'] == reactions_for_sbml[r2]['pr'] \
                        and reactions_for_sbml[r1]['kin'] == reactions_for_sbml[r2]['kin']:
                    simlog_warning('The following reaction: \n' +
                                   f'{reactions_for_sbml[r1]} \n' +
                                   'Is doubled. Was that intentional? \n')

        # Event implementation here
        # Basico does not resolve species that have no reactions but were added to events
        # So we add "phantom" reactions to fix this
        parameters_in_events = set()
        events_for_sbml, species_in_events = eh_format_event_dictionary_for_sbml(species_for_sbml, event_dictionary,
                                                                                 orthogonal_vector_structure,
                                                                                 volume, dimension,
                                                                                 meta_species_to_simulate,
                                                                                 parameter_exist,
                                                                                 parameters_in_events)
        cls.add_to_parameters_to_sbml(parameters_used, parameters_for_sbml, parameters_in_events)

        # Check to see if parameters are names are repeated or used as meta-species
        for p in parameters_for_sbml:
            if p in names_used:
                simlog_error('Parameters names must be unique and they must not share a name with a species')
            names_used.add(p)

        # Store the parameter objects for possible unit conversion - or others
        parameter_object_dict = {}
        for key in parameters_used:
            parameter_object_dict[key] = mp_Mobspy_Parameter.parameter_stack[key]

        species_in_reactions = set()
        for key, reaction in reactions_for_sbml.items():
            for reactant in reaction['re']:
                species_in_reactions.add(reactant[1])
            for product in reaction['pr']:
                species_in_reactions.add(product[1])
        cls.add_phantom_reactions(reactions_for_sbml, species_in_events.difference(species_in_reactions))

        if continuous_sim:
            end_event = {'trigger': ending_condition.generate_string(orthogonal_vector_structure),
                         'delay': '0',
                         'assignments': [('End_Flag_MetaSpecies', '1')]}

            reactions_for_sbml['phantom_reaction_end'] = {'re': [(10, 'End_Flag_MetaSpecies')], 'pr': [],
                                                          'kin': 'End_Flag_MetaSpecies * 1e-100'}
            events_for_sbml['end_event'] = end_event

        set_to_double_parameter = \
            set().union(parameters_in_counts).union(parameters_in_reaction).union(parameters_in_events)
        for p1 in set_to_double_parameter:
            for p2 in set_to_double_parameter:
                if p1 == p2:
                    continue
                else:
                    if p1.name == p2.name:
                        simlog_error('There are two different Parameter Objects with the same name')

        non_processed_assignments = {}
        for spe in meta_species_to_simulate:
            for asgn, expression in spe._assignments.items():
                non_processed_assignments[asgn] = expression
        assignments_for_sbml = \
            asgi_Assign.compile_assignments_for_sbml(non_processed_assignments,
                                                     orthogonal_vector_structure,
                                                     meta_species_to_simulate)

        model_str = ''
        if verbose:
            model_str = '\n'
            model_str += 'Species' + '\n'
            species_alpha = list(sorted(species_for_sbml.keys()))
            for spe in species_alpha:
                model_str += spe.replace('_dot_', '.') + ',' + str(species_for_sbml[spe]) + '\n'

            model_str += '\n'
            model_str += 'Mappings' + '\n'
            mappings_alpha = list(sorted(mappings_for_sbml.keys()))
            for map in mappings_alpha:
                model_str += map + ' :' + '\n'
                for element in sorted(mappings_for_sbml[map]):
                    model_str += element + '\n'

            model_str += '\n'
            model_str += 'Parameters' + '\n'
            parameters_alpha = list(sorted(parameters_for_sbml.keys()))
            for par in parameters_alpha:
                model_str += par + ',' + str(parameters_for_sbml[par][0]) + '\n'

            model_str += '\n'
            model_str += 'Reactions' + '\n'
            remove_phantom_reactions = deepcopy(reactions_for_sbml)
            to_remove = []
            for reaction in remove_phantom_reactions:
                if 'phantom' in reaction:
                    to_remove.append(reaction)
            for r in to_remove:
                remove_phantom_reactions.pop(r, None)
            reaction_alpha = [str(x[1]).replace('_dot_', '.') for x in
                              list(sorted(remove_phantom_reactions.items(), key=lambda x: str(x[1])))]

            for i, reac in enumerate(reaction_alpha):
                model_str += 'reaction_' + str(i) + ',' + reac + '\n'

            if events_for_sbml != {}:
                model_str += '\n'
                model_str += 'Events' + '\n'
                list_to_sort = [str(events_for_sbml[key]) for key in events_for_sbml]
                list_to_sort = sorted(list_to_sort)
                for i in range(len(list_to_sort)):
                    model_str += ('event_' + str(i) + ',' + list_to_sort[i] + '\n').replace('_dot_', '.')

            if assignments_for_sbml != {}:
                model_str += '\n'
                model_str += 'Assignments' + '\n'
                list_to_sort = [str(assignments_for_sbml[key]) for key in assignments_for_sbml]
                list_to_sort = sorted(list_to_sort)
                for i in range(len(list_to_sort)):
                    model_str += ('assignment_' + str(i) + ',' + list_to_sort[i] + '\n').replace('_dot_', '.')

        return species_for_sbml, reactions_for_sbml, parameters_for_sbml, mappings_for_sbml, model_str, \
               events_for_sbml, assigned_species, parameters_used, parameter_object_dict, assignments_for_sbml, \
               has_mole
