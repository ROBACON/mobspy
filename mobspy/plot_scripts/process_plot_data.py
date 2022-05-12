from copy import deepcopy
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.plot_params.example_plot_reader as epr


def query_plot_data(species, data):
    new_data = deepcopy(data)

    species_to_plot = []
    for spe in species:
        for time_series in data:
            for key in time_series.keys():
                if spe == key:
                    species_to_plot.append(key)
                    break

    for spe in species:
        if spe in species_to_plot:
            continue

        species_string_split = spe.split('.')
        for i, time_series in enumerate(data):
            new_data[i][spe] = {'runs': []}
            # Extract necessary data
            for key in time_series:
                if key == 'Time':
                    T = range(len(time_series['Time']))
                else:
                    runs = range(len(time_series[key]['runs']))

            species_to_sum = []
            for key in time_series.keys():
                if all([char in key.split('.') for char in species_string_split]):
                    species_to_sum.append(key)

            # Can be optimized with list usage
            for run in runs:
                this_run = []
                for t in T:
                    mapping_sum = 0
                    for spe_sum in species_to_sum:
                        mapping_sum = mapping_sum + time_series[spe_sum]['runs'][run][t]
                    this_run.append(mapping_sum)
                new_data[i][spe]['runs'].append(this_run)
        species_to_plot.append(spe)

    return species_to_plot, new_data


def check_plot_parameters(species, plot_params):

    dictionary = epr.get_example_plot_parameters()

    if 'Time' in plot_params:
        simlog.error('Time must not be a plot parameter name')

    for spe in species:
        if spe in dictionary.keys():
            simlog.error(f'Plotting is impossible, species {spe} is a parameter name')

    for key in plot_params:
        if key not in dictionary and key not in species:
            simlog.warning(f'Parameter {key} not supported')



