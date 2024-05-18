import mobspy.plot_scripts.statistics_calculations as sc
from copy import deepcopy
import mobspy.plot_scripts.hierarchical_plot as hp
import json
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.plot_scripts.process_plot_data as ppd
from pint import Quantity


def read_plot_json(plot_json_filename):
    """
        This function converts a plot_json file into a dictionary

        :param plot_json_filename: (str) JSON file name

        :return: json_data (dict) converted JSON as dictionary
    """
    with open(plot_json_filename, 'r') as file:
        try:
            json_data = json.load(file)
        except Exception as e:
            simlog.error(f'The following error happened while decoding json file "{plot_json_filename}":\n' +
                         str(e))

    return json_data


def set_plot_units(new_plot_params):
    """
        Sets the plot labels to the unit names by adding to the xlabel and ylabel

        :param new_plot_params: (dict) plot parameters after some changes
    """
    if 'xlabel' not in new_plot_params:
        new_plot_params['xlabel'] = 'Time'

        if new_plot_params['unit_x'] is not None and not ("ignore_unit_label_x" in new_plot_params and
                                                          new_plot_params["ignore_unit_label_x"]):
            if not isinstance(new_plot_params['unit_x'], Quantity):
                new_plot_params['xlabel'] += f' ({new_plot_params["unit_x"]}s)'
            else:
                new_plot_params['xlabel'] += f' ({new_plot_params["unit_x"].units}s)'

    if 'ylabel' not in new_plot_params:
        if new_plot_params['output_concentration']:
            new_plot_params['ylabel'] = 'Conc.'
        else:
            new_plot_params['ylabel'] = 'Counts'

        if new_plot_params['unit_y'] is not None and not ("ignore_unit_label_x" in new_plot_params and
                                                          new_plot_params["ignore_unit_label_y"]):
            if not isinstance(new_plot_params['unit_y'], Quantity):
                new_plot_params['ylabel'] += f' ({new_plot_params["unit_y"]})'
            else:
                new_plot_params['ylabel'] += f' ({new_plot_params["unit_y"].units})'


def stochastic_plot(species, data, plot_params):
    """
        Design default stochastic plot using MobsPy plotting hierarchy. It them passes the parameters for plotting
        in the hierarchical plot module

        :param species: (list of str) list of species names in MobsPy str format (queries are performed with the
        query plot data function)
        :param data: (dict) data in MobsPy format Simulation.results['data']
        :param plot_params: (dict) dictionary with the plot parameters supplied by the user before the modifications
        from this function
    """

    # Data Handling - Copy data object to not interfere with simulation data
    species, data = ppd.query_plot_data(species, data)
    ppd.check_plot_parameters(species, plot_params)

    data_to_plot = data

    try:
        new_plot_params = deepcopy(plot_params)
    except:
        new_plot_params = plot_params
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

            processed_runs = sc.average_plus_standard_deviation(spe, data_to_plot)

            key_average = spe + '$' + 'average'
            key_dev = spe + '$' + 'deviation'

            data_to_plot[0][key_average] = processed_runs[0]
            data_to_plot[0][key_dev] = [processed_runs[1], processed_runs[2]]

            # We define the standard plot for the average and deviation
            plots_for_spe_i.append({'species_to_plot': [spe]})

            plots_for_spe_i_sta.append({'species_to_plot': [key_average], 'time_series': [0]})
            plots_for_spe_i_sta.append({'species_to_plot': [key_dev], 'fill_between': True, 'time_series': [0]})

        except ValueError:
            simlog.error(f'{spe} species not found in data')
        new_plot_params['figures'].append({'ylabel': spe + ' ' + new_plot_params['ylabel'], 'plots': plots_for_spe_i})
        new_plot_params['figures'].append(
            {'ylabel': spe + ' ' + new_plot_params['ylabel'], 'plots': plots_for_spe_i_sta})

    # Setting species parameters
    for spe in species:
        key_average = spe + '$' + 'average'
        key_dev = spe + '$' + 'deviation'

        color = color_cycler(1)

        new_plot_params[key_dev] = {'color': (0.8, 0.8, 0.8), 'linestyle': ':', 'label': 'std. dev'}
        if spe not in new_plot_params:
            new_plot_params[spe] = {'color': color, 'label': spe}
            new_plot_params[key_average] = {'color': color, 'linestyle': '-', 'label': 'mean'}
        else:
            new_plot_params[key_average] = {}

            for par in new_plot_params[spe]:
                new_plot_params[key_average][par] = new_plot_params[spe][par]
                new_plot_params[key_dev][par] = new_plot_params[spe][par]

            if 'label' not in new_plot_params[spe]:
                new_plot_params[spe]['label'] = spe

            if 'color' not in new_plot_params[spe]:
                new_plot_params[spe]['color'] = color
                new_plot_params[key_average]['color'] = color
                new_plot_params[key_average]['linestyle'] = '-'
                new_plot_params[key_average]['label'] = 'mean'
            else:
                new_plot_params[key_average]['color'] = new_plot_params[spe]['color']
                new_plot_params[key_average]['linestyle'] = '-'
                new_plot_params[key_average]['label'] = 'mean'

    return hp.plot_data(data_to_plot, new_plot_params)


def deterministic_plot(species, data, plot_params):
    """
        Design default deterministic plot using MobsPy plotting hierarchy. It them passes the parameters for plotting
        in the hierarchical plot module

        :param species: (list of str) list of species names in MobsPy str format (queries are performed
        with the query plot data function)
        :param data: (dict) data in MobsPy format Simulation.results['data']
        :param plot_params: (dict) dictionary with the plot parameters supplied by the user before the modifications
        from this function
    """
    # Data Handling
    species, data = ppd.query_plot_data(species, data)
    ppd.check_plot_parameters(species, plot_params)

    try:
        new_plot_params = deepcopy(plot_params)
    except:
        new_plot_params = plot_params
    set_plot_units(new_plot_params)

    #  Plot Config
    new_plot_params['frameon'] = False
    new_plot_params['species_to_plot'] = species
    color_cycler = hp.Color_cycle()
    for spe in species:
        if spe not in new_plot_params:
            new_plot_params[spe] = {'label': spe, 'color': color_cycler(1)}
        else:
            if 'label' not in new_plot_params[spe]:
                new_plot_params[spe]['label'] = spe
            if 'color' not in new_plot_params[spe]:
                new_plot_params[spe]['color'] = color_cycler(1)

    return hp.plot_data(data, new_plot_params)


def parametric_plot(species, data, plot_params):

    max_labels = 15
    current_labels = 0

    try:
        new_plot_params = deepcopy(plot_params)
    except:
        new_plot_params = plot_params
    set_plot_units(new_plot_params)

    new_plot_params['frameon'] = False
    new_plot_params['figures'] = []
    new_plot_params['pad'] = 1.5

    # color_cycler = hp.Color_cycle()
    previous_parameter = data.ts_model_parameters[0]

    # Update plot to add new curve
    def update_plot(spe, temp_ts, i, p, previous_parameter, current_labels):
        label = str(spe) + ' ' + str(previous_parameter) if current_labels < max_labels else None
        temp_dict = {'species_to_plot': spe, 'time_series': temp_ts, spe: {'label': label}}
        return temp_dict, [i], p

    # Extracting time-series per parameter in sweep
    for spe in species:
        temp_ts = []
        plots = []
        for i, p in enumerate(data.ts_model_parameters):
            if str(p) == str(previous_parameter):
                temp_ts.append(i)
            else:
                current_labels += 1
                temp_dict, temp_ts, previous_parameter = update_plot(spe, temp_ts, i, p,
                                                                     previous_parameter, current_labels)
                plots.append(temp_dict)

            # For the final value
            if i == len(data) - 1:
                if str(p) != str(previous_parameter):
                    current_labels += 1
                    temp_dict, temp_ts, previous_parameter = \
                        update_plot(spe, [len(data) - 1], len(data) - 1, p, previous_parameter,
                                    current_labels)
                else:
                    current_labels += 1
                    temp_dict, temp_ts, previous_parameter = \
                        update_plot(spe, temp_ts, len(data) - 1, p, previous_parameter,
                                    current_labels)

                current_labels = 0
                plots.append(temp_dict)

        new_plot_params['figures'].append({'plots': plots})

    # Adjusting for label size
    if len(data.ts_model_parameters) < 5:
        prop = {'size': 10}
    elif len(data.ts_model_parameters) < 10:
        prop = {'size': 8}
    else:
        prop = {'size': 6}

    new_plot_params['prop'] = prop
    return hp.plot_data(data, new_plot_params)


def raw_plot(data, parameters_or_file, return_fig=False):
    """
        Plots data from a json or parameter dictionary configured according to the hierarchical plot structure
        Does not accept parameters from a Simulation object, it must be given the parameters in it's entirety

        :param data: (dict) data in MobsPy format Simulation.results['data']
        :param parameters_or_file: (dict, str) Dictionary originated from a JSON or JSON file name
        :param return_fig: (bool) return figure instead of plotting
    """
    if type(parameters_or_file) == str and parameters_or_file[-5:] == '.json':
        plot_params = read_plot_json(parameters_or_file)
    elif type(parameters_or_file) == dict:
        plot_params = parameters_or_file
    else:
        simlog.error('Raw plot only takes json files or parameters for configuration')

    species = list(data.ts_data[0].keys())
    ppd.check_plot_parameters(species, plot_params)

    # Here we add some base data just in case the user has not supplied it
    return hp.plot_data(data, plot_params, return_fig_object=return_fig)


if __name__ == '__main__':
    raw_plot({}, 'default_parameters.json')
