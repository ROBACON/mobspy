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

        temp_list = [0 for _ in range(len(new_data.ts_data))]
        for i, ts in enumerate(new_data.ts_data):
            temp_list[i] = new_data[spe, i]
            species_to_plot.add(spe)

        for res, ts in zip(temp_list, new_data.ts_data):
            ts[spe] = res

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

    for key in plot_params:
        if key not in dictionary and key not in species:
            simlog.warning(f'Parameter {key} not supported')



