"""
    Main MobsPy module. It stocks the Simulation class which is responsible for simulating a Model
"""

# For user imports
from mobspy.modules.meta_class import BaseSpecies, ListSpecies, New, Zero
from mobspy.modules.set_counts_module import set_counts
from mobspy.modules.mobspy_parameters import ModelParameters
from mobspy.modules.assignments_implementation import Assign
from mobspy.modules.order_operators import All, Default, Rev
from mobspy.modules.class_of_meta_specie_named_any import Any
from mobspy.parameter_estimation_data_loader.parameter_estimation_scripts import basiCO_parameter_estimation

# Module imports
from mobspy.modules.meta_class_utils \
    import create_orthogonal_vector_structure as mcu_create_orthogonal_vector_structure
from mobspy.modules.mobspy_expressions import u
from mobspy.modules.logic_operator_objects import MetaSpeciesLogicResolver as lop_MetaSpeciesLogicResolver
from mobspy.modules.compiler import Compiler
from copy import deepcopy
from contextlib import contextmanager
from mobspy.parameter_estimation_data_loader.data_loader \
    import Experimental_Data_Holder as pdl_Experimental_Data_Holder
import mobspy.simulation_logging.log_scripts as simlog
from mobspy.modules.meta_class import Species, Reacting_Species, List_Species
from inspect import stack as inspect_stack
from mobspy.modules.unit_handler import convert_time as uh_convert_time, \
    extract_length_dimension as uh_extract_length_dimension
from mobspy.parameters.default_reader import get_default_parameters
from mobspy.parameters.example_reader import get_example_parameters
from mobspy.plot_params.default_plot_reader import get_default_plot_parameters
from mobspy.plot_params.example_plot_reader import get_example_plot_parameters
from pint import Quantity
from mobspy.parameter_scripts.parameter_reader import parameter_process as pr_parameter_process, \
    manually_process_each_parameter as pr_manually_process_each_parameter, \
    convert_time_parameters_after_compilation as pr_convert_time_parameters_after_compilation, \
    convert_volume_after_compilation as pr_convert_volume_after_compilation,\
    read_json as pr_read_json
from mobspy.parameter_scripts.parametric_sweeps import generate_all_sbml_models as ps_generate_all_sbml_models,\
    unite_parameter_dictionaries as ps_unite_parameter_dictionaries
from joblib import Parallel as joblib_Parallel, delayed as joblib_delayed
from mobspy.sbml_simulator.run import simulate as sbml_simulate
from mobspy.data_handler.process_result_data import extract_time_and_volume_list as dh_extract_time_and_volume_list, \
    convert_data_to_desired_unit as dh_convert_data_to_desired_unit
from mobspy.data_handler.time_series_object import MobsPyTimeSeries, MobsPyList_of_TS
from json import dump as json_dump, load as json_load
from mobspy.plot_scripts.default_plots import deterministic_plot as dp_deterministic_plot, \
    stochastic_plot as dp_stochastic_plot, raw_plot as dp_raw_plot, parametric_plot as dp_parametric_plot
from mobspy.sbml_simulator.builder import build as sbml_build
from os.path import splitext as os_path_splitext


class Simulation(pdl_Experimental_Data_Holder):

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
        """
            Removes the context in all meta-species and resets some varaiables. Called each time an event context is finished. 
        """
        self._event_time = 0
        Species.reset_simulation_context()
        self._context_not_active = True

    def event_context_add(self, time, trigger):
        """
            Adds an event to the event context

            :param trigger: () condition that triggers the event when fulfilled
            :param time: (int, float, Quantity) time to wait before triggering the event once the trigger condition has been fulfilled
        """

        event_data = {'event_time': time, 'event_counts': list(self.current_event_count_data),
                      'trigger': trigger}

        self.current_event_count_data = []
        self.pre_number_of_context_comparisons = self.number_of_context_comparisons
        self.number_of_context_comparisons = 0

        if len(event_data['event_counts']) != 0:
            self.total_packed_events.append(event_data)

        self.event_context_finish()

    def event_context_initiator(self):
        """
            Sets the context in all meta-species. Called each time an event context is initiated.
        """
        Species.set_simulation_context(self)

    def _event_handler(self):
        """
            Handles the event context by activating the current context and checking it is the only one active. It is called in every event context manager.

            :raise simlog.error: if the event context is called although another event is already active
        """
        if self._context_not_active:
            self._context_not_active = False
            self.__dict__['parameters']['_with_event'] = True
            self.event_context_initiator()
        else:
            simlog.error('MobsPy does not support multiple context calls')

    @contextmanager
    def event_condition(self, trigger, delay=0):
        """
            Context manager for condition events. Called in "with Simulation.event_condition(trigger) :" format. 

            :param trigger: () condition that triggers the event when fulfilled
            :param delay: (int, float, Quantity) time to wait before triggering the event once the trigger condition has been fulfilled
        """
        try:
            code_line = inspect_stack()[2].code_context[0][:-1]
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
            delay = uh_convert_time(delay)
            self._conditional_event = False
            self.event_context_add(delay, trigger)

    @contextmanager
    def event_time(self, time):
        """
            Context manager for time events. Called in "with Simulation.event_time(time) :" format. 

            :param time: (int, float, Quantity) any time used
        """
        try:
            self._event_handler()
            yield 0
        finally:
            time = uh_convert_time(time)
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
        super().__init__()

        # Event Variable Definitions
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
        self.model_parameters = {}
        self.sbml_data_list = []
        self._parameter_list_of_dic = []
        self._is_compiled = False
        self.dimension = None
        self.model_parameter_objects_dict = None

        # Must copy to avoid reference assignment
        self.model = List_Species(model)
        self.names = names

        if not isinstance(model, Species) and not isinstance(model, List_Species):
            simlog.error('Model must be formed only by Species objects or List_Species objects \n'
                         f'Model type {type(model)} and it is {model}')

        self.orthogonal_vector_structure = mcu_create_orthogonal_vector_structure(model)

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

        if not parameters:
            self.parameters = get_default_parameters()

        if not plot_parameters:
            self.plot_parameters = get_default_plot_parameters()

        # Other needed things for simulating
        self.results = {}
        self.fres = {}
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
        if self.dimension is None:
            if isinstance(self.volume, Quantity):
                self.dimension = uh_extract_length_dimension(str(self.volume.dimensionality), self.dimension)
            else:
                self.dimension = 3

        simlog.global_simlog_level = self.parameters['level']
        simlog.debug('Compiling model')

        pr_parameter_process(self.parameters)
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
        self.model_string, self._events_for_sbml, self._assigned_species_list, \
        self.model_parameters, self.model_parameter_objects_dict, \
        self._assignments_for_sbml, self._has_mole = \
            Compiler.compile(self.model,
                             reactions_set=self._reactions_set,
                             species_counts=self._species_counts,
                             orthogonal_vector_structure=self.orthogonal_vector_structure,
                             volume=self.parameters['volume'],
                             dimension=self.dimension,
                             type_of_model=self.parameters["rate_type"],
                             verbose=verbose,
                             event_dictionary=self.total_packed_events,
                             continuous_sim=self.parameters['_continuous_simulation'],
                             ending_condition=self.parameters['_end_condition'],
                             skip_expression_check=self.parameters['skip_expression_check'])

        # The volume is converted to the proper unit at the compiler level
        self.parameters['volume'] = self._parameters_for_sbml['volume'][0]
        self.mappings = deepcopy(self._mappings_for_sbml)

        self.all_species_not_mapped = {}
        for key in self._species_for_sbml:
            self.all_species_not_mapped[key.replace('_dot_', '.')] = self._species_for_sbml[key]

        self._list_of_models += [{'species_for_sbml': self._species_for_sbml,
                                  'parameters_for_sbml': self._parameters_for_sbml,
                                  'reactions_for_sbml': self._reactions_for_sbml,
                                  'events_for_sbml': self._events_for_sbml,
                                  'assignments_for_sbml': self._assignments_for_sbml,
                                  'species_not_mapped': self.all_species_not_mapped,
                                  'mappings': self.mappings,
                                  'assigned_species': self._assigned_species_list}]

        self._list_of_parameters = [self.parameters]

        self._is_compiled = True

        if self.model_string != '':
            return self.model_string

    def _assemble_multi_simulation_structure(self):

        if not self.sbml_data_list:
            data_for_sbml_construction, parameter_list_of_dic = ps_generate_all_sbml_models(self.model_parameters,
                                                                                            self._list_of_models)
            self.sbml_data_list = data_for_sbml_construction
            self._parameter_list_of_dic = parameter_list_of_dic

    def run(self, duration=None, volume=None, dimension=None,
            repetitions=None, level=None, simulation_method=None,
            rate_type=None, plot_type=None,
            start_time=None, r_tol=None, a_tol=None, seeds=None, step_size=None,
            jobs=None, unit_x=None, unit_y=None, output_concentration=None, output_event=None,
            output_file=None, save_data=None, plot_data=None) -> None:
        """
            runs the simulation, the arguments for this function can be set using the get item method or this way
            Parameters not set to run will override the previous parameters, causing the model to recompile
            Resulting data is found with the query Simulation_object.results

            :param duration: (float) duration of a simulation
            :param volume: (float) volume of the simulation - if none given 1 - liter is used
            :param repetitions: (int) number of times to reapeat a simulation
            :param level: (int) 0 - only error messages, 3 - errors, warnings, compilation info
            :param simulation_method: (str) stochastic, deterministic, direct_method - simulation method
            :param rate_type: (str) stochastic or deterministic - mass action kinetics or probability of encounter
            :param plot_type: (str) stochastic-style MobsPy plot or deterministic-style MobsPy plot
            :param start_time: (float) the simulation will only display data after the start time
            :param r_tol: (float) relative tolerance - basiCO simulation parameter
            :param a_tol: (float) absolute tolerance  - basiCO simulation parameter
            :param seeds: (list) list of seeds for stochastic simulation
            :param step_size: (float) time step-size for simulation
            :param jobs: (int) number of jobs to execute simulation
            :param unit_x: (unit) unit of the time x-axis
            :param unit_y: (unit) unit of the y-axis
            :param output_concentration: (bool) outputs the concentration instead of counts - to be changed
            :param output_event: (bool) exactly when an event happens, it adds the data point to the results
            :param output_file: (str) name of the file
            :param save_data: (bool) save data or not
            :param plot_data: (bool) plot data or not
        """
        # Level needs to be set before compilation
        if level is not None:
            self.level = level

        # Base case - If there are no events we compile the model here
        if self._species_for_sbml is None:
            self.compile(verbose=False)

        # This is only here so the ide gives the users tips about the function argument.
        # I wish there was a way to loop over all argument without args and kargs
        pr_manually_process_each_parameter(self, duration, volume, dimension,
                                           repetitions, level, simulation_method,
                                           start_time, r_tol, a_tol, seeds, step_size,
                                           jobs, unit_x, unit_y, output_concentration, output_event,
                                           output_file, save_data, plot_data, rate_type, plot_type)

        self._assemble_multi_simulation_structure()

        simlog.debug('Starting Simulator')
        jobs = self.set_job_number(self.parameters)
        simulation_function = lambda x: sbml_simulate(jobs, self._list_of_parameters, x)
        results = joblib_Parallel(n_jobs=jobs, prefer="threads")(joblib_delayed(simulation_function)(sbml)
                                                                 for sbml in self.sbml_data_list)

        simlog.debug("Simulation is Over")

        if self.parameters['unit_y'] is None and self.parameters['output_concentration'] and self._has_mole:
            self.parameters['unit_y'] = 1 * u.unit_registry_object.molar
        elif self.parameters['unit_y'] is None and not self.parameters['output_concentration'] and self._has_mole:
            self.parameters['unit_y'] = 1 * u.unit_registry_object.mol

        # Volume list and time list are to convert into concentrations
        # This section also checks if the output_concentration parameter is valid
        volume_list, time_list, flag_concentration = dh_extract_time_and_volume_list(self._list_of_parameters)
        tcb = self.parameters['unit_y'] is not None and '[length]' not in self.parameters['unit_y'].dimensionality
        if not flag_concentration or tcb:
            self.parameters['output_concentration'] = False

        def convert_one_ts_to_desired_unit(unconverted_data):
            # Convert all the data from a single ts to desired unit
            return dh_convert_data_to_desired_unit(unconverted_data, time_list, volume_list,
                                                   self.parameters['unit_x'],
                                                   self.parameters['unit_y'],
                                                   self.parameters['output_concentration'])

        def convert_all_ts_to_correct_format(single_ts, parameters, unit_convert=False):
            # Convert multiple ts_data into correct format
            if unit_convert:
                data_dict = {'data': convert_one_ts_to_desired_unit(single_ts),
                             'params': self.parameters,
                             'models': self._list_of_models}
            else:
                data_dict = {'data': single_ts,
                             'params': self.parameters,
                             'models': self._list_of_models}
            return MobsPyTimeSeries(data_dict, parameters)

        flatt_ts = []
        if self._parameter_list_of_dic:
            for r, params in zip(results, self._parameter_list_of_dic):
                for ts in r:
                    flatt_ts.append((ts, params))
        else:
            for r in results:
                for ts in r:
                    flatt_ts.append((ts, {}))

        ta = self.parameters['unit_x'] is not None
        tb = self.parameters['unit_y'] is not None
        tc = self.parameters['output_concentration'] if flag_concentration else False

        if ta or tb or tc:
            all_processed_data = joblib_Parallel(n_jobs=jobs, prefer="threads") \
                (joblib_delayed(convert_all_ts_to_correct_format)(ts, params, True) for ts, params in flatt_ts)
        else:
            all_processed_data = joblib_Parallel(n_jobs=jobs, prefer="threads") \
                (joblib_delayed(convert_all_ts_to_correct_format)(ts, params, False) for ts, params in flatt_ts)

        self.results = MobsPyList_of_TS(all_processed_data, self.model_parameter_objects_dict)
        self.fres = MobsPyList_of_TS([all_processed_data[0]], None, True)

        if self.parameters['save_data']:
            self.save_data()

        # Set common simulation and plot parameters
        self.plot_parameters['unit_x'] = self.parameters['unit_x']
        self.plot_parameters['unit_y'] = self.parameters['unit_y']
        self.plot_parameters['output_concentration'] = self.parameters['output_concentration']

        if self.parameters['plot_data']:

            methods_list = [x['plot_type'] for x in self._list_of_parameters]

            if len(self._parameter_list_of_dic) > 1:
                self.plot_parametric()
                return 0

            if 'stochastic' in methods_list:
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
                    json_dump(self.results.to_dict(), f, indent=4)
            except Exception as e:
                simlog.warning("Error saving data. Potential solve: file name parameter")
                simlog.warning(str(e))
        else:
            file += '.json'
            with open(file, 'w') as jf:
                json_dump(self.results.to_dict(), jf, indent=4)

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
            data = json_load(json_file)
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
                      'total_packed_events', 'species_initial_counts',
                      '_event_time', 'previous_trigger', 'current_event_count_data',
                      'current_condition', 'current_event_trigger_data',
                      'number_of_context_comparisons', 'pre_number_of_context_comparisons', '_continuous_simulation',
                      'initial_duration', '_reactions_set', '_list_of_models', '_list_of_parameters',
                      '_context_not_active', '_species_counts', '_assigned_species_list', '_conditional_event',
                      '_end_condition', 'orthogonal_vector_structure', 'model_parameters', 'fres',
                      'sbml_data_list', '_parameter_list_of_dic', '_is_compiled', 'dimension', 'experimental_data',
                      'model_parameter_objects_dict', '_assignments_for_sbml', '_has_mole']

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
                # If the model is already compiled, the change in parameters should be faster
                if self._is_compiled and name != 'unit_x' and name != 'unit_y':
                    value = pr_convert_time_parameters_after_compilation(value)
                if self._is_compiled and name == 'volume':
                    value = pr_convert_volume_after_compilation(self.dimension, self._parameters_for_sbml, value)

                if name == 'duration':
                    if type(value) == bool:
                        simlog.error(f'MobsPy has received an invalid trigger type: {type(value)} \n' +
                                     f'Please make sure you are not using the operator == for ' +
                                     f'creating event conditions \n'
                                     , stack_index=2)

                if name == 'duration' and isinstance(value, lop_MetaSpeciesLogicResolver):
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
        ta = item == 'results' and self.__dict__['results'] == {}
        tb = item == 'fres' and self.__dict__['fres'] == {}
        if ta or tb:
            simlog.error('The results were accessed before the execution of the simulation', stack_index=2)

        if item == 'plot_config':
            return self.__getattr__(item)

        return super().__getattribute__(item)

    def __getattr__(self, item):
        """
            __getattr__ override. For the user to be able to set plot parameters as MySim.plot.parameter
        """
        if item == 'plot_config':
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
            if os_path_splitext(config)[1] != '.json':
                simlog.error('Wrong file extension', stack_index=3)
            parameters_to_config = pr_read_json(config)
        elif type(config) == dict:
            parameters_to_config = config
        else:
            simlog.error("Parameters must be python dictionary or json file", stack_index=3)
        return parameters_to_config

    def add_plot_params(self, *args, **kwargs):

        for a in args:
            if type(a) == dict:
                for par in a:
                    self.base_sim.plot_parameters[par] = a[par]

        for key in kwargs:
            self.plot_parameters[key] = deepcopy(kwargs[key])

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
        return dp_stochastic_plot(plot_essentials[0], plot_essentials[1], plot_essentials[2])

    def plot_deterministic(self, *species):
        """
            Calls deterministic plot. See default_plots module in the plot_scripts directory

            :param species: (str or meta-species objects) list of species to be plotted
        """
        plot_essentials = self.extract_plot_essentials(*species)
        return dp_deterministic_plot(plot_essentials[0], plot_essentials[1], plot_essentials[2])

    def plot_parametric(self, *species):
        plot_essentials = self.extract_plot_essentials(*species)
        return dp_parametric_plot(plot_essentials[0], plot_essentials[1], plot_essentials[2])

    def plot(self, *species):
        """
            Another way of calling plot_deterministic for simplicity

            :param species: (str or meta-species objects) list of species to be plotted
        """
        return self.plot_deterministic(*species)

    def plot_raw(self, parameters_or_file, return_fig=False):
        """
            Calls raw plot. See default_plots module in the plot_scripts directory

            :param parameters_or_file: json file name with plot parameter configuration or dictionary with plot
            parameter configuration
        """
        return dp_raw_plot(self.results, parameters_or_file, return_fig=return_fig)

    def __add__(self, other):
        """
            The add operator is used to concatenate simulations.
        """
        return SimulationComposition(self, other)

    def generate_sbml(self):
        """
            Generates sbmls strings from the current stored models in the simulation

            "return: to_return (list of str) list of sbml files from all the simulations stored
        """
        to_return = []
        if self._species_for_sbml is None:
            self.compile(verbose=False)
        self._assemble_multi_simulation_structure()

        for parameter_sweep in self.sbml_data_list:
            for sbml_data in parameter_sweep:
                to_return.append(sbml_build(sbml_data['species_for_sbml'], sbml_data['parameters_for_sbml'],
                                            sbml_data['reactions_for_sbml'], sbml_data['events_for_sbml'],
                                            sbml_data['assignments_for_sbml']))
        return to_return

    @classmethod
    def is_simulation(cls):
        return True

    @classmethod
    def set_job_number(cls, params):
        # Run in parallel or sequentially
        # If nothing is specified just run it in parallel
        try:
            if params["jobs"] == 1:
                jobs = params["jobs"]
            else:
                jobs = params["jobs"]
        except KeyError:
            jobs = -1
        return jobs


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
                                             f' remain the same')

    def __len__(self):
        return len(self.list_of_simulations)

    def __iter__(self):
        for sim in self.list_of_simulations:
            yield sim

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
        self.fres = None
        self.base_sim = self.list_of_simulations[0]

    def __add__(self, other):
        return SimulationComposition(self, other)

    def __setattr__(self, name, value):
        white_list = ['list_of_simulations', 'results', 'base_sim', 'fres']
        multi_cast_parameters = ['duration']
        broad_cast_parameters = ['level', 'rate_type', 'plot_type', 'repetitions']
        double_cast_parameters = ['simulation_method', 'volume', 'method']

        if name in double_cast_parameters:

            # Broadcast if single value
            if type(value) == str or type(value) == int or type(value) == float or isinstance(value, Quantity):
                for sim in self:
                    sim.__dict__['parameters'][name] = value
            else:
                # Multicast if list
                try:
                    if not len(self) == len(value):
                        raise SystemExit
                except:
                    simlog.error(f'The parameter {name} was assigned non-accepted type.', stack_index=2)

                for par, sim in zip(value, self):
                    # DON'T ADD DIRECTLY to the object's dict, changes in volume, duration after compilation are
                    # checked in the setattr method in the individual simulations
                    if name == 'volume':
                        sim.volume = par
                    elif name == 'duration':
                        sim.duration = par
                    else:
                        sim.__dict__['parameters'][name] = par

        elif name in multi_cast_parameters:

            try:
                if not len(self) == len(value):
                    raise SystemExit
            except:
                simlog.error('From 2.4.4 duration must be assigned to each simulation individually or a list '
                             'with all durations must be assigned to the concatenated simulation', stack_index=2)

            for par, sim in zip(value, self):
                # DON'T ADD DIRECTLY to the object's dict, changes in volume, duration after compilation are
                # checked in the setattr method in the individual simulations
                if name == 'volume':
                    sim.volume = par
                elif name == 'duration':
                    sim.duration = par
                else:
                    sim.__dict__['parameters'][name] = par
        elif name in broad_cast_parameters:
            for sim in self:
                sim.__dict__['parameters'][name] = value
        else:
            if name in white_list:
                self.__dict__[name] = value
            else:
                self.base_sim.__setattr__(name, value)

    def __getattr__(self, item):
        if item == 'plot_config':
            self.base_sim.__dict__['plot_flag'] = True
            return self.base_sim

    def compile(self, verbose=True):
        str = ''
        for sim in self.list_of_simulations:
            str += sim.compile(verbose)

        self._compile_multi_simulation()
        self.base_sim._assemble_multi_simulation_structure()
        if str != '':
            return str

    def _check_all_sims_compilation(self):

        for sim in self.list_of_simulations:
            if sim._species_for_sbml is None:
                sim.compile(verbose=False)

    # This run is for the multiple simulations
    def run(self, duration=None, volume=None, dimension=None,
            repetitions=None, level=None, simulation_method=None,
            start_time=None, r_tol=None, a_tol=None, seeds=None, step_size=None,
            jobs=None, unit_x=None, unit_y=None, output_concentration=None, output_event=None,
            output_file=None, save_data=None, plot_data=None, rate_type=None, plot_type=None):
        """
            runs a concatenated simulation with multiple simulation objects summed
            Some inputs are different here, duration and volume must receive any iterable with the volume for
            each simulation

            :param duration: (iterable) duration of a simulation
            :param volume: (iterable) volume of the simulation - if none given 1 - liter is used
            :param repetitions: (int) number of times to reapeat a simulation
            :param level: (int) 0 - only error messages, 3 - errors, warnings, compilation info
            :param simulation_method: (iterable) stochastic, deterministic, direct_method - simulation method
            :param start_time: (float) the simulation will only display data after the start time
            :param r_tol: (float) relative tolerance - basiCO simulation parameter
            :param a_tol: (float) absolute tolerance  - basiCO simulation parameter
            :param seeds: (list) list of seeds for stochastic simulation
            :param step_size: (float) time step-size for simulation
            :param jobs: (int) number of jobs to execute simulation
            :param unit_x: (unit) unit of the time x-axis
            :param unit_y: (unit) unit of the y-axis
            :param output_concentration: (bool) outputs the concentration instead of counts - to be changed
            :param output_event: (bool) exactly when an event happens, it adds the data point to the results
            :param output_file: (str) name of the file
            :param save_data: (bool) save data or not
            :param plot_data: (bool) plot data or not
            :param rate_type: (str) stochastic or deterministic rate expression - mass action kinetic or prob. based
            :param plot_type: (str) stochastic or deterministic - style of MobsPy plot
        """
        if level is not None:
            self.level = level

        self._check_all_sims_compilation()
        self._compile_multi_simulation()

        pr_manually_process_each_parameter(self, duration, volume, dimension,
                                           repetitions, level, simulation_method,
                                           start_time, r_tol, a_tol, seeds, step_size,
                                           jobs, unit_x, unit_y, output_concentration, output_event,
                                           output_file, save_data, plot_data, rate_type, plot_type)

        multi_parameter_dictionary = {}

        for sim in self.list_of_simulations:
            multi_parameter_dictionary = ps_unite_parameter_dictionaries(multi_parameter_dictionary,
                                                                         sim.model_parameters)

        self.base_sim.model_parameters = multi_parameter_dictionary

        for sim in self.list_of_simulations:

            if sim == self.base_sim:
                continue

            self.base_sim._list_of_models += sim._list_of_models
            self.base_sim._list_of_parameters += sim._list_of_parameters

        self.base_sim.run()
        self.results = self.base_sim.results
        self.fres = self.base_sim.fres

    def plot_deterministic(self, *species):
        self.base_sim.plot_deterministic(*species)

    def plot_stochastic(self, *species):
        self.base_sim.plot_stochastic(*species)

    def plot(self, *species):
        self.base_sim.plot(*species)

    def plot_raw(self, parameters_or_file):
        self.base_sim.plot_raw(parameters_or_file)

    def add_plot_params(self, *args, **kwargs):

        for a in args:
            if type(a) == dict:
                for par in a:
                    self.base_sim.plot_parameters[par] = a[par]

        for key in kwargs:
            self.base_sim.plot_parameters[key] = deepcopy(kwargs[key])

    def generate_sbml(self):

        self._check_all_sims_compilation()
        self._compile_multi_simulation()

        for sim in self.list_of_simulations:

            if sim == self.base_sim:
                continue

            self.base_sim._list_of_models += sim._list_of_models
            self.base_sim._list_of_parameters += sim._list_of_parameters

        return self.base_sim.generate_sbml()

    @classmethod
    def is_simulation(cls):
        return True


if __name__ == '__main__':
    pass
