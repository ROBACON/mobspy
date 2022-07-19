"""
    Main MobsPy module. It stocks the Simulation class which is responsible for simulating a Model
"""
from mobspy.modules import meta_class
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
import json
import os
import inspect
from pint import UnitRegistry

# u is reserved for units
u = UnitRegistry()


class Simulation:

    def __init__(self, model, names=None, parameters=None, plot_parameters=None):
        """
            Constructor of the simulation object

            Parameters:
            model (ParallelSpecies object) = Meta-Species object for modeling
            names (dict) = names of the meta-species in globals() format
            parameters (dict) = Simulation object parameters
            plot_parameters (dict) = Parameters for plotting
        """
        if names is None:
            local_names = inspect.stack()[1][0].f_locals
            global_names = inspect.stack()[1][0].f_globals
            names = {}
            for key, item in global_names.items():
                names[key] = item
            for key, item in local_names.items():
                names[key] = item

        self.model = model
        self.names = names

        if not isinstance(model, Species) and not isinstance(model, ParallelSpecies):
            simlog.error('Model must be formed by Species objects')

        if not parameters:
            self.parameters = get_default_parameters()

        if not plot_parameters:
            self.plot_parameters = get_default_plot_parameters()

        # Other needed things for simulating
        self.sbml_string = None
        self.results = None
        self.packed_data = []
        self.default_order = Default

        self._species_for_sbml = None
        self._reactions_for_sbml = None
        self._parameters_for_sbml = None
        self._mappings_for_sbml = None
        self.model_string = ''

    def compile(self, verbose=True):
        """
            Compiler method that calls the Compiler class in the modules directory

            Parameters:
                verbose (bool) = print or not the results of the compilation
        """
        simlog.global_simlog_level = self.parameters['level']
        simlog.debug('Compiling model')

        pr.parameter_process(self.parameters)
        if self.parameters['simulation_method'].lower() == 'deterministic':
            self.parameters['repetitions'] = 1
            self.plot_parameters['simulation_method'] = 'deterministic'
        elif self.parameters['simulation_method'].lower() == 'stochastic':
            self.plot_parameters['simulation_method'] = 'stochastic'

        self._species_for_sbml, self._reactions_for_sbml, \
        self._parameters_for_sbml, self._mappings_for_sbml, \
        self.model_string = Compiler.compile(self.model, names=self.names,
                                             volume=self.parameters['volume'],
                                             type_of_model=self.parameters[
                                                 "simulation_method"],
                                             verbose=verbose,
                                             default_order=self.default_order)

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

        self.sbml_string = sbml_builder.build(self._species_for_sbml,
                                              self._parameters_for_sbml,
                                              self._reactions_for_sbml)

        if self.model_string != '':
            return self.model_string

        return self.sbml_string

    def run(self):
        """
            Runs the simulation by colling the models in the sbml_simulator directory.
            Compiles the model if it was not yet compiled
        """
        if self._species_for_sbml is None:
            self.compile(verbose=False)

        simlog.debug('Starting Simulator')
        self.sbml_string = sbml_builder.build(self._species_for_sbml,
                                              self._parameters_for_sbml,
                                              self._reactions_for_sbml)

        self.results = sbml_run.simulate(self.sbml_string, self.parameters, self.mappings, self.all_species_not_mapped)
        self.results['data'] = dh.convert_data_to_desired_unit(self.results['data'],
                                                               self.parameters['unit_x'], self.parameters['unit_y'],
                                                               self.output_concentration, self.parameters['volume'])

        self._pack_data(self.results['data'])

        if self.parameters['save_data']:
            simlog.debug("Saving data (reason: parameter <save_data>)")
            pickled = {'data': self.results,
                       'mappings': self._mappings_for_sbml,
                       'params': self.parameters}
            # Save pickle data
            # TODO Add error here
            if not os.path.isdir(self.parameters['output_absolute_directory']):
                simlog.debug("Creating output directory: %s..." % (self.parameters['output_absolute_directory']))
                os.makedirs(self.parameters['output_absolute_directory'], exist_ok=True)

            with open(self.parameters['output_absolute_file'], 'w') as jf:
                json.dump(pickled, jf, indent=4)
        else:
            simlog.warning("NOT saving data (reason: parameter <save_data>)")

        if self.parameters['plot_data']:
            if self.plot_parameters['simulation_method'] == 'stochastic':
                self.plot_stochastic()
            elif self.plot_parameters['simulation_method'] == 'deterministic':
                self.plot_deterministic()

    def save_results(self, file):
        """
            Save results manually into file. Useful for jupyter notebook users

            Parameters
                file (str) = name of the file to create and save JSON data
        """
        simlog.warning('Only THIS model data will be saved not externally added data')
        with open(file, 'w') as jf:
            json.dump(self.results, jf, indent=4)

    def _pack_data(self, time_series_data):
        """
            Packs data from multiple simulations or external data into one simulation object

            Parameters:
                time_series_data (data in MobsPy format) = data to be packed in the simulation object
        """
        self.packed_data.append(time_series_data)

    def add_simulation_data(self, results):
        """
            Encapsulation. See pack data

            Parameters:
                results (dict) = result data MobsPy.results
        """
        self._pack_data(results['data'])
        self.mappings = results['mappings']

    def add_external_data(self, data):
        """
            Encapsulation. See pack data

            Parameters:
                results (dict) = result data MobsPy.results["data"]
        """
        self._pack_data(data)

    # Dealing with parameters
    def set_from_json(self, file_name):
        """
            Set simulation parameters from json file

            Parameters:
                file_name (str) = name of the json file
        """
        with open(file_name) as json_file:
            data = json.load(json_file)
            for key in data:
                self.__setattr__(key, data[key])

    def __setattr__(self, name, value):
        """
            __setattr__ override. For seeting simulation parameters using the _dot_ operator

            Parameters:
                name (str) = name of the parameter to set
                value = value of the parameter
        """
        white_list = ['default_order', 'volume', 'model', 'names', 'parameters', 'model_string',
                      'plot_parameters', 'sbml_string', 'results', 'packed_data', '_species_for_sbml',
                      '_reactions_for_sbml', '_parameters_for_sbml', '_mappings_for_sbml', 'mappings',
                      'all_species_not_mapped', 'self._species_for_sbml', 'self._reactions_for_sbml',
                      'self._parameters_for_sbml', 'self._mappings_for_sbml', 'self.model_string']
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
                self.__dict__['parameters'][name] = value
            elif name in white_list:
                pass
            else:
                simlog.error(f'Parameter {name} is not supported')

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
            Configure simulation parameters from json file. Different from set as it overrides previous parameters

            Parameters:
                file_name (str) = name of the json file
        """
        self.parameters = self.__config_parameters(config)

    def configure_plot_parameters(self, config):
        """
            Configure plot parameters from json file. Different from set as it overrides previous parameters

            Parameters:
                file_name (str) = name of the json file
        """
        self.plot_parameters = self.__config_parameters(config)

    @staticmethod
    def __config_parameters(config):
        """
            Encapsulation for config_plot and config_parameters
        """
        if type(config) == str:
            if os.path.splitext(config)[1] != '.json':
                simlog.error('Wrong file extension')
            parameters_to_config = pr.read_json(config)
        elif type(config) == dict:
            parameters_to_config = config
        else:
            simlog.error("Parameters must be python dictionary or json file")
        return parameters_to_config

    # Plotting encapsulation
    def extract_plot_essentials(self, *species):
        """
            Extract essential information for the hierarchical plotting tool

            Parameters:
                *species (meta-species objects) = meta-species objects to plot
        """
        if not species:
            species_strings = list(self.mappings.keys())
        else:
            species_strings = []

        for spe in species:
            if isinstance(spe, Species):
                species_strings.append(str(spe))
            elif isinstance(spe, Reacting_Species):
                species_strings.append(str(spe))
            elif type(spe) == str:
                species_strings.append(spe)
            else:
                simlog.error('Only species objects or strings for plotting arguments')
        return species_strings, self.packed_data, self.plot_parameters

    def plot_stochastic(self, *species):
        """
            Calls stochastic plot. See default_plots module in the plot_scripts directory
        """
        plot_essentials = self.extract_plot_essentials(*species)
        dp.stochastic_plot(plot_essentials[0], plot_essentials[1], plot_essentials[2])

    def plot_deterministic(self, *species):
        """
            Calls deterministic plot. See default_plots module in the plot_scripts directory
        """
        plot_essentials = self.extract_plot_essentials(*species)
        dp.deterministic_plot(plot_essentials[0], plot_essentials[1], plot_essentials[2])

    def plot_raw(self, parameters_or_file):
        """
            Calls raw plot. See default_plots module in the plot_scripts directory
        """
        dp.raw_plot(self.packed_data, parameters_or_file)


if __name__ == '__main__':
    pass
