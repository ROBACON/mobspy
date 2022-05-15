import mobspy.plot_scripts.statistics_calculations as sc
from copy import deepcopy
import mobspy.plot_scripts.hierarchical_plot as hp
import json
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.plot_scripts.process_plot_data as ppd


def read_plot_json(plot_json_filename):
    """
        This function converts a plot_json file into a dictionary

        Parameters:
            plot_json_filename (str) = JSON file name

        Returns:
            json_data (dict) = Converted JSON as dictionary
    """
    with open(plot_json_filename, 'r') as file:
        try:
            json_data = json.load(file)
        except Exception:
            simlog.error(f'Error while decoding json file "{plot_json_filename}".')

    return json_data


def set_plot_units(new_plot_params):
    """
        Sets the plot labels to the unit names by adding to the xlabel and ylabel

        Parameters:
            new_plot_params (dict) = plot parameters after some changes
    """
    new_plot_params['xlabel'] = 'Time'
    if new_plot_params['unit_x'] is not None:
        new_plot_params['xlabel'] += f'({new_plot_params["unit_x"]})'
    else:
        new_plot_params['xlabel'] += '(second)'

    if new_plot_params['output_concentration']:
        new_plot_params['ylabel'] = 'Conc.'
    else:
        new_plot_params['ylabel'] = 'Counts'

    if new_plot_params['unit_y'] is not None:
        new_plot_params['ylabel'] += f'({new_plot_params["unit_y"]})'
    else:
        if new_plot_params['output_concentration']:
            new_plot_params['ylabel'] += '(molar)'


def stochastic_plot(species, data, plot_params):
    """
        Design default stochastic plot using MobsPy plotting hierarchy. It them passes the parameters for plotting
        in the hierarchical plot module

        Parameters:
            species (list of str) = list of species names in MobsPy str format (queries are performed with the query plot
            data function)
            data (dict) = data in MobsPy format Simulation.results['data']
            plot_params (dict) = dictionary with the plot parameters supplied by the user before the modifications
            from this function
    """
    # Data Handling
    species, data = ppd.query_plot_data(species, data)
    ppd.check_plot_parameters(species, plot_params)

    data_to_plot = deepcopy(data)
    new_plot_params = deepcopy(plot_params)
    set_plot_units(new_plot_params)

    # Plot config
    new_plot_params['frameon'] = False
    new_plot_params['figures'] = []
    new_plot_params['pad'] = 1.5
    color_cycler = hp.Color_cycle()
    for spe in species:
        # We define new 'mappings' with the resulting runs for the statics for the plot structure
        try:
            plots_for_spe_i = []
            plots_for_spe_i_sta = []
            for i, time_series in enumerate(data):
                if len(time_series) == 1:
                    i = ''

                processed_runs = sc.average_plus_standard_deviation(time_series[spe]['runs'])
                key_average = spe + '$' + 'average'
                key_dev = spe + '$' + 'deviation'
                data_to_plot[i][key_average] = {'runs': [processed_runs[0]]}
                data_to_plot[i][key_dev] = {'runs': [processed_runs[1], processed_runs[2]]}

                # We define the standard plot for the average and deviation
                color = color_cycler(1)
                plots_for_spe_i.append({'species_to_plot': [spe], 'time_series': i + 1,
                                        spe: {'color': color, 'label': spe}})

                plots_for_spe_i_sta.append({'species_to_plot': [key_average], 'time_series': i + 1,
                                            key_average: {'color': color, 'linestyle': '-', 'label': 'mean'}})
                plots_for_spe_i_sta.append({'species_to_plot': [key_dev], 'time_series': i + 1, 'fill_between': True,
                                            key_dev: {'color': (0.8, 0.8, 0.8), 'linestyle': ':', 'label': 'std. dev'}})
        except KeyError:
            simlog.error(f'{spe} species not found in data')
        new_plot_params['figures'].append({'ylabel': spe + ' ' + new_plot_params['ylabel'], 'plots': plots_for_spe_i })
        new_plot_params['figures'].append({'ylabel': spe + ' ' + new_plot_params['ylabel'], 'plots': plots_for_spe_i_sta})

    hp.plot_data(data_to_plot, new_plot_params)


def deterministic_plot(species, data, plot_params):
    """
        Design default deterministic plot using MobsPy plotting hierarchy. It them passes the parameters for plotting
        in the hierarchical plot module

        Parameters:
            species (list of str) = list of species names in MobsPy str format (queries are performed with the query plot
            data function)
            data (dict) = data in MobsPy format Simulation.results['data']
            plot_params (dict) = dictionary with the plot parameters supplied by the user before the modifications
            from this function
    """
    # Data Handling
    species, data = ppd.query_plot_data(species, data)
    ppd.check_plot_parameters(species, plot_params)

    new_plot_params = deepcopy(plot_params)
    set_plot_units(new_plot_params)

    #  Plot Config
    new_plot_params['frameon'] = False
    new_plot_params['species_to_plot'] = species
    color_cycler = hp.Color_cycle()
    for spe in species:
        new_plot_params[spe] = {'label': spe, 'color': color_cycler(1)}
    hp.plot_data(data, new_plot_params)


def raw_plot(data, parameters_or_file):
    """
        Plots data from a json or parameter dictionary configured according to the hierarchical plot structure
        Does not accept parameters from a Simulation object, it must be given the parameters in it's entirety

        Parameters:
            data (dict) = data in MobsPy format Simulation.results['data']
            parameters_or_file (dict, str) = Dictionary originated from a JSON or JSON file name
    """
    if type(parameters_or_file) == str and parameters_or_file[-5:] == '.json':
        plot_params = read_plot_json(parameters_or_file)
    elif type(parameters_or_file) == dict:
        plot_params = parameters_or_file
    else:
        simlog.error('Raw plot only takes json files or parameters for configuration')

    species = list(data[0].keys())
    ppd.check_plot_parameters(species, plot_params)
    hp.plot_data(data, plot_params)


if __name__ == '__main__':
    raw_plot({}, 'default_parameters.json')
