import json
import os.path
import math
import mobspy.plot_scripts.statistics_calculations as spd
import mobspy.simulation_logging.log_scripts as simlog

import numpy as np
import pickle as pkl
from pathlib import Path
import inspect
import matplotlib.pyplot as plt


####################### PRACTICAL FUNCTIONS
class Color_cycle:
    """
    :return: a simple color cycle fuction for different colors for different species
    """

    def __init__(self):
        self.index = 0
        self.color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    def __call__(self, n):
        self.index = (self.index + n) % len(self.color_list)
        return self.color_list[self.index]


def find_species_time_series(spe, data):
    for time_series in data:
        if spe in time_series:
            yield time_series


# This function only works for the standart time we are currently using
def convert_to_time_unit(time_data, unit):
    """
    :param time_data: the time data for the used axis in minutes
    :param unit: the unit desired to be converted
    :return: list of time values on the desired unit
    """

    unit = unit.lower()
    time_unit = unit
    if unit == 's':
        factor = 60
    elif unit == 'min':
        factor = 1
    elif unit == 'h':
        factor = 1 / 60
    elif unit == 'd':
        factor = 1 / 60 * 1 / 24
    else:
        simlog.error("Time unit not recognized")

    # use np for speed
    time_data = list(np.array(time_data) * factor)
    return time_data, time_unit


def get_total_figure_number(axis_matrix):
    # Get the number of figures, using the axis_matrix
    try:
        total_figure_number = axis_matrix.shape[0] * axis_matrix.shape[1]
    except IndexError:
        total_figure_number = axis_matrix.shape[0]
    return total_figure_number


# Hash for converting linear figure number into index
def figure_hash(current_figure, axis_matrix):
    """
        This function allows one to acess the figure grid with a linear input
        For instance one can access a 2x2 grid using 0, 1, 2, 3
        0 becomes 0,0
        1 becomes 1,0
        2 becomes 0,1
        3 becomes 1,1

    :param current_figure: linear number of the figure
    :param axis_matrix: a list with all the created axis on the multiple figure subplot
    :return: the correct axis based on the number provided
    """

    # Get the number of lines
    max_lines = len(axis_matrix)
    total_figure_number = get_total_figure_number(axis_matrix)

    # Correction for index purposes
    current_figure = current_figure

    col = math.floor(current_figure / max_lines)

    if total_figure_number <= max_lines:
        return axis_matrix[int(current_figure % max_lines)]
    else:
        return axis_matrix[int(current_figure % max_lines), int(col)]


# Hash to convert total figure number into grid
def figure_hash_creation(total_figure_number, max_lines=None):
    """

        This is hash used to create the figure grid automatically according to the number of figures
        and the maximum number of lines

        For instance 4 figures with 2 as max_lines creates a 2x2 figure grid automatically

    :param total_figure_number: Number of total figures to create
    :param max_lines: Maximum number of lines in the grid
    :return: Figures grid will all the respective axis
    """

    # Default value if nothing is set up
    if max_lines is None:
        max_lines = 2

    if total_figure_number <= max_lines:
        fig, axs = plt.subplots(total_figure_number)
    else:
        column_number = int((total_figure_number + total_figure_number % max_lines) / max_lines)
        fig, axs = plt.subplots(max_lines, column_number)

    if total_figure_number == 1:
        axs = np.array([axs])

    return fig, axs


def find_parameter(params, key, index=None):
    """

        This is the hearth of the plotting structure, this functions allows one to simply set multiple characteristics
        for only one figure

        The priority for parameter search is plots => figures => global, with plot overrinding others and so on

        If a parameter is defined globaly it will be applied to all figures, if it defined inside a figure element
        it will only apply to that figure, if it is defined in a plot element it will only apply to that curve

        Check the readme or the tutorials for more details on the plotting structure. It is simple and versatile

    :param params: Plot parameters from python dictionary (after json conversion)
    :param key: Key necessary to acess the parameters
    :param index: None for global search, one index for figure search, and two for figure curve search
    :return: the parameter if found, and None if nothing is found
    """

    # No index is given, look global
    if index is None:
        try:
            return params[key]
        except (KeyError, IndexError):
            return None
    # Search a parameter
    # If local return local, otherwise return global
    # If not found return nothing
    elif type(index) == int:

        try:
            return params['figures'][index][key]
        except (KeyError, IndexError):

            try:
                return params[key]
            except (KeyError, IndexError):
                return None
    # If two indexes are given
    # Look inside the plot
    # Than figure
    # Than global
    elif type(index) == tuple:

        try:
            return params['figures'][index[0]]['plots'][index[1]][key]

        except (KeyError, IndexError):

            try:
                return params['figures'][index[0]][key]

            except (KeyError, IndexError):

                try:
                    return params[key]
                except (KeyError, IndexError):
                    return None


####################### PLOTING FUNCTIONS

def plot_curves(data, axs, figure_index, plot_params):
    """

        This function plots the curves in the assigned axis
        Also it looks for the final parameters inside the plot json, at the plot level

    :param axs: axs to plot the data in
    :param data: data given to be plot from pickle file
    :param figure_index: index of the current figure analysed
    :param plot_params: parameters for plotting
    :param species: species mappings set by the user
    """

    # Get the plot number from the list of plots
    try:
        plot_number = len(find_parameter(plot_params, 'plots', figure_index))
    except TypeError:
        # No figures plot only the default
        plot_number = 1
    if plot_number == 0:
        plot_number = 1

    # For all plots in the figure
    # Set all parameters and plot
    legend_flag = False
    for plot_index in range(plot_number):

        # Get the species from plot parameters
        if find_parameter(plot_params, key='species_to_plot', index=(figure_index, plot_index)) is not None:
            species = find_parameter(plot_params, key='species_to_plot', index=(figure_index, plot_index))
        else:
            simlog.error('No species found for plotting in the Plotting Parameters')

        # Get the time series to plot
        if find_parameter(plot_params, key='time_series', index=(figure_index, plot_index)) is not None:
            time_series = find_parameter(plot_params, key='time_series', index=(figure_index, plot_index))
            if type(time_series) == int:
                time_series = [time_series]
        else:
            time_series = list(range(len(data)))

        for ts in time_series:
            ts = data[ts - 1]

            for spe in species:

                color_index = 0

                # Get the parameters assigned to the species, if not assign empty for None returns
                if find_parameter(plot_params, key=spe, index=(figure_index, plot_index)) is not None:
                    species_characteristics = find_parameter(plot_params, key=spe, index=(figure_index, plot_index))
                else:
                    species_characteristics = {}

                # Now we search the species characteristics for parameters
                if find_parameter(species_characteristics, key='color') is not None:
                    curve_color = find_parameter(species_characteristics, key='color')
                else:
                    color_index = (color_index + 1) % len(Color_cycle().color_list)
                    curve_color = Color_cycle()(color_index)

                if find_parameter(species_characteristics, key='linestyle') is not None:
                    linestyle = find_parameter(species_characteristics, key='linestyle')
                else:
                    linestyle = '-'

                if find_parameter(species_characteristics, key='linewidth') is not None:
                    linewidth = find_parameter(species_characteristics, key='linewidth')
                else:
                    linewidth = None

                if find_parameter(species_characteristics, key='label') is not None:
                    label = find_parameter(species_characteristics, key='label')
                    legend_flag = True
                else:
                    label = None

                if find_parameter(plot_params, key='runs', index=(figure_index, plot_index)) is not None:
                    runs = find_parameter(plot_params, key='runs', index=(figure_index, plot_index))
                else:
                    runs = range(len(ts[spe]['runs']))

                if find_parameter(plot_params, key='fill_between', index=(figure_index, plot_index)) is not None and \
                        find_parameter(plot_params, key='fill_between', index=(figure_index, plot_index)):
                    try:
                        axs.fill_between(ts['Time'], ts[spe]['runs'][0], ts[spe]['runs'][1], color=curve_color,
                                         label=label)
                    except IndexError:
                        simlog.error('Fill_between must only have two or less runs referring to it')
                else:
                    for run in runs:
                        axs.plot(ts['Time'], ts[spe]['runs'][run], color=curve_color,
                            linestyle=linestyle, linewidth=linewidth, label=label)
                        label = None

    if legend_flag:
        if find_parameter(plot_params, key='frameon', index=figure_index) is not None:
            frameon = find_parameter(plot_params, key='frameon', index=figure_index)
        else:
            frameon = True
        axs.legend(frameon=frameon)


def set_figure_characteristics(axis_matrix, plot_params):
    """

        Sets the figure parameters

    :param axis_matrix: list of all axis in the grid
    :param params: plot parameters received
    :return: nothing axis work by register, so their methods change them outside the function scope
    """
    total_figure_number = get_total_figure_number(axis_matrix)
    # Loop through all axis
    for i in range(total_figure_number):

        if find_parameter(plot_params, 'xlim', i) is not None:
            figure_hash(i, axis_matrix).set_xlim(find_parameter(plot_params, 'xlim', i))

        if find_parameter(plot_params, 'ylim', i) is not None:
            figure_hash(i, axis_matrix).set_ylim(find_parameter(plot_params, 'ylim'))

        if find_parameter(plot_params, 'logscale', i) is not None and 'X' in find_parameter(plot_params, 'logscale', i):
            figure_hash(i, axis_matrix).set_xscale('log')

        if find_parameter(plot_params, 'logscale', i) is not None and 'Y' in find_parameter(plot_params, 'logscale', i):
            figure_hash(i, axis_matrix).set_yscale('log')

        if find_parameter(plot_params, 'title', i) is not None:
            figure_hash(i, axis_matrix).set_title(find_parameter(plot_params, 'title', i))

        if find_parameter(plot_params, 'xlabel', i) is not None:
            figure_hash(i, axis_matrix).set_xlabel(find_parameter(plot_params, 'xlabel', i))

        if find_parameter(plot_params, 'ylabel', i) is not None:
            figure_hash(i, axis_matrix).set_ylabel(find_parameter(plot_params, 'ylabel', i))

        if find_parameter(plot_params, 'annotations', i) is not None:
            for annotations in find_parameter(plot_params, 'annotations', i):
                figure_hash(i, axis_matrix).text(annotations[0], annotations[1], annotations[2])


def set_global_parameters(fig, plot_params):
    if find_parameter(plot_params, 'pad') is not None:
        fig.tight_layout(pad=find_parameter(plot_params, 'pad'))

    if find_parameter(plot_params, 'figsize') is not None:
        fig.set_size_inches(plot_params['figsize'][0], plot_params['figsize'][1])

    if find_parameter(plot_params, 'dpi') is not None:
        fig.set_dpi(plot_params['dpi'])


def plot_data(data, plot_params):
    """
        This function plots the simulation results according to the specifications

    :param file_name: pickle file resulting from simulation name
    :param module_logger: logger from simulation to store plotting info
    :param plot_json: plot_json file with plotting instructions
    """

    # Get the figure number from the list of figures
    # Add it to parameters
    try:
        figure_number = len(find_parameter(plot_params, 'figures'))
    except TypeError:
        # No figures plot only the default
        figure_number = 1
    if figure_number == 0:
        figure_number = 1

    # Create figures grid, max_lines is 2 as default and the value is defined in the function
    fig, axis_matrix = figure_hash_creation(figure_number, max_lines=find_parameter(plot_params, 'max_lines'))

    # Set parameters common to all figures, like padding and figure size.
    set_global_parameters(fig, plot_params)

    # Global call for setting global parameters
    set_figure_characteristics(axis_matrix, plot_params)

    # Now we plot
    for figure_index in range(figure_number):
        plot_curves(data, figure_hash(figure_index, axis_matrix), figure_index, plot_params)

    # Real default case
    plt.show()


if __name__ == "__main__":
    pass
