import modules.meta_class
import sbml_simulator.run
import simulation_logging.log_scripts as simlog
from modules.meta_class import *
from parameter_scripts import parameter_reader as pr
from parameters.default_reader import get_default_parameters
from parameters.example_reader import get_example_parameters
from plot_params.default_plot_reader import get_default_plot_parameters
import sbml_simulator.builder as sbml_builder
import plot_scripts.default_plots as dp
import json
import os
import inspect


class Simulation:

    def __init__(self, model, names=None, parameters=None, plot_parameters=None):
        """
            This part of the code just creates a model
        :param model: Species or Parallel species object instance for modeling
        :param names: names of the variables. Set to globals() for easy naming
        :param parameters: Parameters for a model
        :param plot_parameters: Plot parameters for plotting
        """
        if names is None:
            names = inspect.stack()[1][0].f_globals

        self.model = model
        self.names = names
        self.plot_call_flag = False

        if not isinstance(model, Species) and not isinstance(model, ParallelSpecies):
            simlog.error('Model must be formed by Species objects')

        if not parameters:
            self.parameters = get_default_parameters()

        if not plot_parameters:
            self.plot_parameters = get_default_plot_parameters()

        simlog.debug('Compiling model')
        self.compile(verbose=False)
        self.mappings = deepcopy(self._mappings_for_sbml)

        # Other needed things for simulating
        self.sbml_string = None
        self.results = None
        self.packed_data = []

    def compile(self, verbose=True):
        self._species_for_sbml, self._reactions_for_sbml, \
        self._parameters_for_sbml, self._mappings_for_sbml = Compiler.compile(self.model, names=self.names,
                                                                              volume_ml=self.parameters['volume_ml'],
                                                                              type_of_model=self.parameters[
                                                                                  "simulation_method"],
                                                                              verbose=verbose)

    def run(self):
        """
            Just calls the simulator part of the codes for running
        :return: nothing, data is saved automaticaly or in self.results
        """
        # We process the parameters here in case there were updates
        pr.parameter_process(self.parameters, self._mappings_for_sbml, self._parameters_for_sbml)
        if self.parameters['simulation_method'].lower() == 'deterministic':
            self.parameters['repetitions'] = 1
            self.plot_parameters['simulation_method'] = 'deterministic'
        elif self.parameters['simulation_method'].lower() == 'stochastic':
            self.plot_parameters['simulation_method'] = 'stochastic'

        simlog.debug('Starting Simulator')
        self.sbml_string = sbml_builder.build(self._species_for_sbml,
                                              self._parameters_for_sbml,
                                              self._reactions_for_sbml)

        self.results = sbml_simulator.run.simulate(self.sbml_string, self.parameters, self._mappings_for_sbml)
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

    def save_data(self):
        simlog.warning('Only THIS model data will be saved not externally added data')
        with open(self.parameters['output_absolute_file'], 'w') as jf:
            json.dump(self.results, jf, indent=4)

    def _pack_data(self, time_series_data):
        self.packed_data.append(time_series_data)

    def add_simulation_data(self, results):
        self._pack_data(results['data'])
        self.mappings = results['mappings']

    def add_external_data(self, data):
        self._pack_data(data)

    # Dealing with parameters
    def __setattr__(self, name, value):

        if 'plot_flag' in self.__dict__ and self.__dict__['plot_flag']:
            self.__dict__["plot_parameters"][name] = value
            self.__dict__["plot_flag"] = False

        example_parameters = get_example_parameters()
        if name in example_parameters.keys():
            self.__dict__['parameters'][name] = value
        else:
            self.__dict__[name] = value

    def __getattr__(self, item):
        if item == 'plot':
            self.__dict__['plot_flag'] = True
        else:
            self.__dict__['plot_flag'] = False
        return self

    def configure_parameters(self, config):
        self.parameters = self.__config_parameters(config)

    def configure_plot_parameters(self, config):
        self.plot_parameters = self.__config_parameters(config)

    @staticmethod
    def __config_parameters(config):
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
        if not species:
            species_strings = list(self.mappings.keys())

        for spe in species:
            species_strings = []
            if isinstance(spe, Species):
                species_strings.append(spe.get_name())
            elif type(spe) == str:
                species_strings.append(spe)
            else:
                simlog.error('Only species objects or strings for plotting arguments')
        return species_strings, self.packed_data, self.plot_parameters

    def plot_stochastic(self, *species):
        plot_essentials = self.extract_plot_essentials(*species)
        dp.stochastic_plot(plot_essentials[0], plot_essentials[1], plot_essentials[2])

    def plot_deterministic(self, *species):
        plot_essentials = self.extract_plot_essentials(*species)
        dp.deterministic_plot(plot_essentials[0], plot_essentials[1], plot_essentials[2])

    def plot_raw(self, parameters_or_file):
        dp.raw_plot(self.packed_data, parameters_or_file)


if __name__ == '__main__':
    pass
