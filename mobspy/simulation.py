"""
    Main MobsPy module. It stocks the Simulation class which is responsible for simulating a Model
"""
from copy import deepcopy
from contextlib import contextmanager
from mobspy.modules import meta_class
from mobspy.modules import event_functions
from mobspy.sbml_simulator import run
import mobspy.simulation_logging.log_scripts as simlog
from mobspy.modules.meta_class import *
from mobspy.parameter_scripts import parameter_reader as pr
from mobspy.parameters.default_reader import get_default_parameters
from mobspy.parameters.example_reader import get_example_parameters
from mobspy.plot_params.default_plot_reader import get_default_plot_parameters
import mobspy.sbml_simulator.builder as sbml_builder
import mobspy.sbml_simulator.run as sbml_run
import mobspy.plot_scripts.default_plots as dp
import mobspy.data_handler.process_result_data as dh
from mobspy.modules.set_counts_module import set_counts
from mobspy.data_handler.time_series_object import MobsPyTimeSeries
import json
import os
import inspect
import mobspy.modules.unit_handler as uh
from pint import UnitRegistry

# u is reserved for units
u = UnitRegistry()


class Simulation:

    # Event Implementation
    @classmethod
    def event_compilation_error(cls):
        simlog.error('The event condition did not compile.\n'
                     'Please make sure it follows the following format:\n'
                     'For simple conditions - if C1 \n'
                     'For and based condition - if (C1) & (C2)\n'
                     'For or based conditions - if (C1) & (C2)\n'
                     'Please include the parentheses')

    def event_context_finish(self):
        self._event_time = 0
        for meta_species in self._species_to_set:
            meta_species.reset_simulation_context()
        self._species_to_set = set()
        self._context_not_active = True

    def event_context_add(self, time, trigger):

        event_data = {'event_time': time, 'event_counts': list(self.current_event_count_data),
                      'trigger': trigger}

        self.current_event_count_data = []
        self.pre_number_of_context_comparisons = self.number_of_context_comparisons
        self.number_of_context_comparisons = 0

        if len(event_data['event_counts']) != 0:
            self.total_packed_events.append(event_data)

        self.event_context_finish()

    def event_context_initiator(self):

        # Set context in all meta-species
        if len(self._species_to_set) == 0:
            for i in range(len(inspect.stack())):
                local_names = inspect.stack()[i][0].f_locals
                global_names = inspect.stack()[i][0].f_globals
                for key, item in global_names.items():
                    if isinstance(item, Species):
                        self._species_to_set.add(item)
                for key, item in local_names.items():
                    if isinstance(item, Species):
                        self._species_to_set.add(item)

            for meta_species in self._species_to_set:
                meta_species.set_simulation_context(self)
        else:
            pass

    def _event_handler(self):
        if self._context_not_active:
            self._context_not_active = False
            self.__dict__['parameters']['_with_event'] = True
            self.event_context_initiator()
        else:
            simlog.error('MobsPy does not support multiple context calls')

    @contextmanager
    def event_condition(self, trigger, delay=0):
        try:
            code_line = inspect.stack()[2].code_context[0][:-1]
            if '==' in code_line:
                simlog.error('Equality comparison operator (==) not allowed for MobsPy events \n' +
                             'Please use (A <= n) & (A >= n) if necessary', stack_index=3)
            if type(trigger) == bool or type(trigger) == float or type(trigger) == int:
                simlog.error(f'MobsPy has received an invalid trigger type: {type(trigger)} \n' +
                             f'Please make sure you are not using the operator == for creating event conditions \n'
                             , stack_index=3)
            self._conditional_event = True
            self._event_handler()
            yield 0
        finally:
            delay = uh.convert_time(delay)
            self._conditional_event = False
            self.event_context_add(delay, trigger)

    @contextmanager
    def event_time(self, time):
        try:
            self._event_handler()
            yield 0
        finally:
            time = uh.convert_time(time)
            self.event_context_add(time, 'true')

    def __init__(self, model, names=None, parameters=None, plot_parameters=None):
        """
            Constructor of the simulation object

            Parameters:
            :param model: (List_Species object) Meta-Species object for modeling
            :param names: (dict) names of the meta-species in globals() format. If none it uses the variable names
            :param parameters: (dict) Simulation object parameters. If none takes default parameters
            :param plot_parameters: (dict) Parameters for plotting. If none takes default
        """

        # Event Variable Definitions
        self._species_to_set = set()
        self._event_time = 0
        self.previous_trigger = None
        self.current_event_count_data = []
        self.total_packed_events = []
        self.number_of_context_comparisons = 0
        self.pre_number_of_context_comparisons = 0
        self._list_of_models = []
        self._list_of_parameters = []
        self._context_not_active = True
        self._assigned_species_list = []
        self._conditional_event = False
        self._end_condition = None

        # HERE FABRICIO
        # Get all names
        #if names is None:
        #    local_names = inspect.stack()[1][0].f_locals
        #    global_names = inspect.stack()[1][0].f_globals
        #    names = {}
        #    for key, item in global_names.items():
        #        names[key] = item
        #    for key, item in local_names.items():
        #        names[key] = item

        self.model = model
        self.names = names

        if not isinstance(model, Species) and not isinstance(model, List_Species):
            simlog.error('Model must be formed only by Species objects or List_Species objects \n'
                         f'Model type {type(model)} and it is {model}')

        self.orthogonal_vector_structure = mcu.create_orthogonal_vector_structure(model)

        # Get all meta - reactions
        self._reactions_set = set()
        for spe_object in self.model:
            for reference in spe_object.get_references():
                self._reactions_set = self._reactions_set.union(reference.get_reactions())

        self._species_counts = []
        for spe_object in self.model:
            for count in spe_object.get_quantities():
                self._species_counts.append({'object': spe_object, 'characteristics': count['characteristics'],
                                             'quantity': count['quantity']})
            spe_object.reset_counts()

        if not parameters:
            self.parameters = get_default_parameters()

        if not plot_parameters:
            self.plot_parameters = get_default_plot_parameters()

        # Other needed things for simulating
        self.results = {}
        self.default_order = Default

        self._species_for_sbml = None
        self._reactions_for_sbml = None
        self._parameters_for_sbml = None
        self._mappings_for_sbml = None
        self._events_for_sbml = None
        self.model_string = ''

    def compile(self, verbose=True):
        """
            Compiler method that calls the Compiler class in the modules directory

            :param verbose: (bool) = print or not the results of the compilation
        """
        simlog.global_simlog_level = self.parameters['level']
        simlog.debug('Compiling model')

        pr.parameter_process(self.parameters)
        if self.parameters['method'] is not None:
            self.parameters['simulation_method'] = self.parameters['method']

        if self.parameters['simulation_method'].lower() == 'deterministic':
            self.plot_parameters['simulation_method'] = 'deterministic'
        elif self.parameters['simulation_method'].lower() == 'stochastic':
            self.plot_parameters['simulation_method'] = 'stochastic'

        # Pass end condition to dict parameters - It is stored outside of parameters to the parameters serializable
        # However, it is necessary for the compilation so it is passed as a parameter
        self.parameters['_end_condition'] = self._end_condition

        self._species_for_sbml, self._reactions_for_sbml, \
        self._parameters_for_sbml, self._mappings_for_sbml, \
        self.model_string, self._events_for_sbml, self._assigned_species_list = \
            Compiler.compile(self.model,
                             reactions_set=self._reactions_set,
                             species_counts=self._species_counts,
                             orthogonal_vector_structure=self.orthogonal_vector_structure,
                             volume=self.parameters['volume'],
                             type_of_model=self.parameters[
                                 "simulation_method"],
                             verbose=verbose,
                             event_dictionary=self.total_packed_events,
                             continuous_sim=self.parameters['_continuous_simulation'],
                             ending_condition=self.parameters['_end_condition'])

        # The volume is converted to the proper unit at the compiler level
        self.parameters['volume'] = self._parameters_for_sbml['volume'][0]
        self.mappings = deepcopy(self._mappings_for_sbml)

        # Set common parameters for plot and simulation
        self.plot_parameters['unit_x'] = self.parameters['unit_x']
        self.plot_parameters['unit_y'] = self.parameters['unit_y']
        self.plot_parameters['output_concentration'] = self.parameters['output_concentration']

        self.all_species_not_mapped = {}
        for key in self._species_for_sbml:
            self.all_species_not_mapped[key.replace('_dot_', '.')] = self._species_for_sbml[key]

        self._list_of_models += [{'species_for_sbml': self._species_for_sbml,
                                  'parameters_for_sbml': self._parameters_for_sbml,
                                  'reactions_for_sbml': self._reactions_for_sbml,
                                  'events_for_sbml': self._events_for_sbml,
                                  'species_not_mapped': self.all_species_not_mapped,
                                  'mappings': self.mappings,
                                  'assigned_species': self._assigned_species_list}]

        self._list_of_parameters = [self.parameters]

        if self.model_string != '':
            return self.model_string

    def run(self):
        """
            Runs the simulation by colling the models in the sbml_simulator directory.
            Compiles the model if it was not yet compiled
        """
        # Base case - If there are no events we compile the model here
        if self._species_for_sbml is None:
            self.compile(verbose=False)

        simlog.debug('Starting Simulator')

        unprocessed_data = sbml_run.simulate(self._list_of_parameters, self._list_of_models)

        processed_data = []
        for updl in unprocessed_data:
            processed_data.append(dh.convert_data_to_desired_unit(updl, self.parameters['unit_x'],
                                                                  self.parameters['unit_y'],
                                                                  self.output_concentration,
                                                                  self.parameters['volume']))

        data_dict = {'data': processed_data,
                     'params': self.parameters,
                     'models': self._list_of_models}

        self.results = MobsPyTimeSeries(data_dict)

        if self.parameters['save_data']:
            self.save_data()

        if self.parameters['plot_data']:
            methods_list = [x['simulation_method'] for x in self._list_of_parameters]

            if ('stochastic' in methods_list or 'directedMethod' in methods_list) \
                    and self.parameters['repetitions'] > 1:
                self.plot_stochastic()
            else:
                self.plot_deterministic()

    def save_data(self, file=None):
        """
            Saves the simulation result data to a file in json format

            :param file: (str) name of the file to save the data to. If none a default name is provided
        """
        self._save_data(file=file)

    def _save_data(self, file=None):
        """
            Save results manually into file. Useful for jupyter notebook users

            Parameters
                file (str) = name of the file to create and save JSON data
        """
        if file is None:
            try:
                with open(self.parameters["absolute_output_file"], 'w') as f:
                    json.dump(self.results.to_dict(), f, indent=4)
            except Exception as e:
                simlog.warning("Error saving data. Potential solve: file name parameter")
                simlog.warning(str(e))
        else:
            file += '.json'
            with open(file, 'w') as jf:
                json.dump(self.results.to_dict(), jf, indent=4)

    def _pack_data(self, time_series_data):
        """
            Packs data from multiple simulations or external data into one simulation object

            :param time_series_data: (data in MobsPy format) data to be packed in the simulation object
        """
        self.packed_data.append(time_series_data)

    # Dealing with parameters
    def set_from_json(self, file_name):
        """
            Set simulation parameters from json file

            :param file_name: (str) name of the json file
        """
        with open(file_name) as json_file:
            data = json.load(json_file)
            for key in data:
                self.__setattr__(key, data[key])

    def __setattr__(self, name, value):
        """
            __setattr__ override. For setting simulation parameters using the _dot_ operator

            :param name: (str) name of the parameter to set
            :param value: value of the parameter
        """
        white_list = ['default_order', 'volume', 'model', 'names', 'parameters', 'model_string',
                      'plot_parameters', 'results', '_species_for_sbml',
                      '_reactions_for_sbml', '_parameters_for_sbml', '_mappings_for_sbml', 'mappings',
                      'all_species_not_mapped', 'self._species_for_sbml', 'self._reactions_for_sbml',
                      'self._parameters_for_sbml', 'self._mappings_for_sbml', 'self.model_string',
                      'event_times', 'event_models', 'event_count_dics', '_events_for_sbml',
                      'total_packed_events', 'species_initial_counts', '_species_to_set',
                      '_event_time', 'previous_trigger', 'current_event_count_data',
                      'current_condition', 'current_event_trigger_data',
                      'number_of_context_comparisons', 'pre_number_of_context_comparisons', '_continuous_simulation',
                      'initial_duration', '_reactions_set', '_list_of_models', '_list_of_parameters',
                      '_context_not_active', '_species_counts', '_assigned_species_list', '_conditional_event',
                      '_end_condition', 'orthogonal_vector_structure']

        plotted_flag = False
        if name in white_list:
            self.__dict__[name] = value

        if 'plot_flag' in self.__dict__ and self.__dict__['plot_flag']:
            self.__dict__["plot_parameters"][name] = value
            self.__dict__["plot_flag"] = False
            plotted_flag = True

        if not plotted_flag:
            example_parameters = get_example_parameters()
            if name in example_parameters.keys():
                if name == 'duration':
                    if type(value) == bool:
                        simlog.error(f'MobsPy has received an invalid trigger type: {type(value)} \n' +
                                     f'Please make sure you are not using the operator == for ' +
                                     f'creating event conditions \n'
                                     , stack_index=2)

                if name == 'duration' and isinstance(value, MetaSpeciesLogicResolver):
                    self.__dict__['parameters']['_continuous_simulation'] = True
                    self.__dict__['_end_condition'] = value
                    if 'initial_conditional_duration' not in self.__dict__['parameters']:
                        self.__dict__['parameters']['initial_conditional_duration'] = 1
                else:
                    self.__dict__['parameters'][name] = value
            elif name in white_list:
                pass
            else:
                simlog.error(f'Parameter {name} is not supported', stack_index=2)

    def __getattribute__(self, item):
        if item == 'results' and self.__dict__['results'] == {}:
            simlog.error('The results were accessed before the execution of the simulation', stack_index=2)
        return super().__getattribute__(item)

    def __getattr__(self, item):
        """
            __getattr__ override. For the user to be able to set plot parameters as MySim.plot.parameter
        """
        if item == 'plot':
            self.__dict__['plot_flag'] = True
        else:
            self.__dict__['plot_flag'] = False
        return self

    def configure_parameters(self, config):
        """
            Configure simulation parameters from json file or dictionary

            :param file_name: (str) name of the json file
        """
        self.parameters = self.__config_parameters(config)

    def configure_plot_parameters(self, config):
        """
            Configure plot parameters from json file or dictionary

            :param file_name: (str) name of the json file
        """
        self.plot_parameters = self.__config_parameters(config)

    @staticmethod
    def __config_parameters(config):
        """
            Encapsulation for config_plot and config_parameters
        """
        if type(config) == str:
            if os.path.splitext(config)[1] != '.json':
                simlog.error('Wrong file extension', stack_index=3)
            parameters_to_config = pr.read_json(config)
        elif type(config) == dict:
            parameters_to_config = config
        else:
            simlog.error("Parameters must be python dictionary or json file", stack_index=3)
        return parameters_to_config

    # Plotting encapsulation
    def extract_plot_essentials(self, *species):
        """
            Extract essential information for plotting

            :param species: (meta-species objects) meta-species objects to plot
            :return: species_strings (str) = species strings to be plotted, self.results = data resulting from the
            simulation, self.plot_parameters (dict) = parameters for plotting
        """
        if not species:
            species_strings = set()
            for model in self._list_of_models:
                species_strings = species_strings.union(model['mappings'])
        else:
            species_strings = set()

        for spe in species:
            if isinstance(spe, Species):
                species_strings.add(str(spe))
            elif isinstance(spe, Reacting_Species):
                species_strings.add(str(spe))
            elif type(spe) == str:
                species_strings.add(spe)
            else:
                simlog.error('Only species objects or strings for plotting arguments', stack_index=4)

        return species_strings, self.results, self.plot_parameters

    def plot_stochastic(self, *species):
        """
            Calls stochastic plot. See default_plots module in the plot_scripts directory

            :param species: (str or meta-species objects) list of species to be plotted
        """
        plot_essentials = self.extract_plot_essentials(*species)
        dp.stochastic_plot(plot_essentials[0], plot_essentials[1], plot_essentials[2])

    def plot_deterministic(self, *species):
        """
            Calls deterministic plot. See default_plots module in the plot_scripts directory

            :param species: (str or meta-species objects) list of species to be plotted
        """
        plot_essentials = self.extract_plot_essentials(*species)
        dp.deterministic_plot(plot_essentials[0], plot_essentials[1], plot_essentials[2])

    def plot(self, *species):
        """
            Another way of calling plot_deterministic for simplicity

            :param species: (str or meta-species objects) list of species to be plotted
        """
        self.plot_deterministic(*species)

    def plot_raw(self, parameters_or_file):
        """
            Calls raw plot. See default_plots module in the plot_scripts directory

            :param parameters_or_file: json file name with plot parameter configuration or dictionary with plot
            parameter configuration
        """
        dp.raw_plot(self.results, parameters_or_file)

    def __add__(self, other):
        return SimulationComposition(self, other)

    def generate_sbml(self):
        """
            Generates sbmls strings from the current stored models in the simulation

            "return: to_return (list of str) list of sbml files from all the simulations stored
        """
        to_return = []
        for sbml_data in self._list_of_models:
            to_return.append(sbml_builder.build(sbml_data['species_for_sbml'], sbml_data['parameters_for_sbml'],
                                                sbml_data['reactions_for_sbml'], sbml_data['events_for_sbml']))
        return to_return


class SimulationComposition:

    def _compile_multi_simulation(self):
        for sim1 in self.list_of_simulations:
            for sim2 in self.list_of_simulations:
                if sim1 == sim2:
                    continue

                for spe1 in sim1.model:
                    for spe2 in sim2.model:

                        if spe1.get_name() == spe2.get_name():
                            if spe1.get_all_characteristics() != spe2.get_all_characteristics():
                                simlog.error(f'Species {spe1.get_name()} was modified through simulations. \n' +
                                             f'Although reactions can be removed, the characteristics inherited must'
                                             f'remains the same')

    def __init__(self, S1, S2):
        if isinstance(S1, Simulation) and isinstance(S2, Simulation):
            self.list_of_simulations = [S1] + [S2]
        elif isinstance(S1, SimulationComposition) and isinstance(S2, Simulation):
            self.list_of_simulations = S1.list_of_simulations + [S2]
        elif isinstance(S1, Simulation) and isinstance(S2, SimulationComposition):
            self.list_of_simulations = [S1] + S2.list_of_simulations
        elif isinstance(S1, SimulationComposition) and isinstance(S2, SimulationComposition):
            self.list_of_simulations = S1.list_of_simulations + S2.list_of_simulations
        else:
            simlog.error('Simulation compositions can only be performed with other simulations', stack_index=3)
        self.results = None
        self.base_sim = self.list_of_simulations[0]

    def __add__(self, other):
        return SimulationComposition(self, other)

    # FIX THIS
    def __setattr__(self, name, value):
        white_list = ['list_of_simulations', 'results', 'base_sim']
        broad_cast_parameters = ['level', 'method', 'volume']

        if name == 'duration':
            simlog.error('The durations are to be defined specifically to each simulation and not for the concatenated'
                         ' object. \n'
                         'Please set the durations for each simulation object independently', stack_index=2)

        if name in broad_cast_parameters:
            if name == 'volume':
                for sim in self.list_of_simulations:
                    if sim.__dict__['parameters'][name] != 1:
                        simlog.error('Volumes must be defined only individually for each simulation or once in the '
                                     'concatenated simulation', stack_index=2)
            else:
                for sim in self.list_of_simulations:
                    sim.__dict__['parameters'][name] = value

        if name in white_list:
            self.__dict__[name] = value
        else:
            self.base_sim.__setattr__(name, value)

    def compile(self, verbose=True):
        str = ''
        for sim in self.list_of_simulations:
            str += sim.compile(verbose)
        self._compile_multi_simulation()
        if str != '':
            return str

    def run(self):
        for sim in self.list_of_simulations:
            if sim._species_for_sbml is None:
                sim.compile(verbose=False)
        self._compile_multi_simulation()

        for sim in self.list_of_simulations:
            if sim == self.base_sim:
                continue

            self.base_sim._list_of_models += sim._list_of_models
            self.base_sim._list_of_parameters += sim._list_of_parameters
        self.base_sim.run()
        self.results = self.base_sim.results

    def plot_deterministic(self, *species):
        self.base_sim.plot_deterministic(*species)

    def plot_stochastic(self, *species):
        self.base_sim.plot_stochastic(*species)

    def plot(self, *species):
        self.base_sim.plot(*species)

    def plot_raw(self, parameters_or_file):
        self.base_sim.plot_raw(parameters_or_file)

    def generate_sbml(self):
        return self.base_sim.generate_sbml()


if __name__ == '__main__':
    pass
