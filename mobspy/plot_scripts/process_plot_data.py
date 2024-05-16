from copy import deepcopy
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.plot_params.example_plot_reader as epr


def query_plot_data(species, data):
    """
        Performs a query of the plot data, when one wishes to plot species with characteristics
        It creates a new data structure with the results from the query added to it for the plotting structure

        :param species: (str) Species name in string format
        :param data: (dict) Data in MobsPy dictionary format
    """
    new_data = deepcopy(data)

    species_to_plot = set()
    for time_series in data.ts_data:
        for key in time_series.keys():
            if key in species:
                species_to_plot.add(key)

    for spe in species:
        if spe in species_to_plot:
            continue

        for i, _ in enumerate(data.ts_data):
            new_data[i][spe] = data[spe][i]

        species_to_plot.add(spe)

    species_to_plot = sorted([spe for spe in species])
    return species_to_plot, new_data


def check_plot_parameters(species, plot_params):
    """
        Performs a check of the plot_parameters given. To see if the parameters are correctly named

        :param species: (str) Species in str format
        :param plot_params: (dict) Plot parameter dictionary
    """
    dictionary = epr.get_example_plot_parameters()

    if 'Time' in plot_params:
        simlog.error('Time must not be a plot parameter name')

    for spe in species:
        if spe in dictionary.keys():
            simlog.error(f'Plotting is impossible, species {spe} is a parameter name')

    # Check if parameters are valid
    validated_keys = set()
    for key in plot_params:
        if key not in dictionary and key not in species:
            continue
        else:
            validated_keys.add(key)

    # Check if query is present
    for key in plot_params:
        if key in validated_keys:
            continue
        else:
            spe_name = key.split('.')[0]
            if spe_name not in species:
                simlog.warning(f'Parameter {key} not supported')
            validated_keys.add(key)


def time_filter_operation(low, high, time_data, data):

    new_time_data = []
    new_data = []

    for t, d in zip(time_data, data):
        if t < low:
            continue
        elif low < t < high:
            new_time_data.append(t)
            new_data.append(d)
        elif t > high:
            break

    return new_time_data, new_data


def y_filter_operation(low_y, high_y, time_data, data):

    new_time_data = []
    new_data = []

    for t, d in zip(time_data, data):
        if d < low_y:
            continue
        elif low_y <= d <= high_y:
            new_time_data.append(t)
            new_data.append(d)
        elif d > high_y:
            continue

    return new_time_data, new_data





