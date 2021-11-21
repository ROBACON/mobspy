import json
import os.path
import math

import numpy as np

import logging
import pickle as pkl
import logging
from pathlib import Path
import inspect
import matplotlib.pyplot as plt
import numpy

logging.basicConfig(level=logging.WARNING)
module_logging = logging.getLogger(__name__)
logging_array = [module_logging]

# Plot color global variable for cycling
__plot_color = 0


####################### PRACTICAL FUNCTIONS

def color_cycle():

    """
    :return: a simple color cycle fuction for different colors for different species
    """

    color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

    global __plot_color
    __plot_color = __plot_color + 1

    return color_list[__plot_color % len(color_list)]


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
        logging_array[0].error("Time unit not recognized. Plotting in minutes instead")
        factor = 1
        time_unit = 'min'

    # use np for speed
    time_data = list(np.array(time_data) * factor)
    return time_data, time_unit


# Hash for converting linear figure number into index
def figure_hash(current_figure, axis_list):
    """
        This function allows one to acess the figure grid with a linear input
        For instance one can access a 2x2 grid using 0, 1, 2, 3
        0 becomes 0,0
        1 becomes 1,0
        2 becomes 0,1
        3 becomes 1,1

    :param current_figure: linear number of the figure
    :param axis_list: a list with all the created axis on the multiple figure subplot
    :return: the correct axis based on the number provided
    """

    # Get the number of lines
    max_lines = len(axis_list)

    # Get the number of figures
    try:
        total_figure_number = axis_list.shape[0] * axis_list.shape[1]
    except IndexError:
        total_figure_number = axis_list.shape[0]

    # Correction for index purposes
    current_figure = current_figure - 1

    col = math.floor(current_figure / max_lines)
    if total_figure_number <= max_lines:
        return axis_list[int(current_figure % max_lines)]
    else:
        return axis_list[int(current_figure % max_lines), int(col)]


# Hash to convert total figure number into grid
def figure_hash_creation(total_figure_number, max_lines = None):
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
        column_number = int((total_figure_number + total_figure_number % max_lines)/max_lines)
        fig, axs = plt.subplots(max_lines, column_number)

    if total_figure_number == 1:
        axs = np.array([axs])

    return fig, axs


def find_parameter(params, key, index = None):
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
        except KeyError:
            return None
    # Search a parameter
    # If local return local, otherwise return global
    # If not found return nothing
    elif type(index) == int:

        try:
            return params['figures'][index][key]
        except KeyError:

            try:
                return params[key]
            except KeyError:
                return None
    # If two indexes are given
    # Look inside the plot
    # Than figure
    # Than global
    elif type(index) == tuple:

        try:
            return params['figures'][index[0]]['plots'][index[1]][key]

        except KeyError:

            try:
                return params['figures'][index[0]][key]

            except KeyError:

               try:
                   return params[key]
               except KeyError:
                   return None


################## DEALING WITH JSONS


def get_plot_params(plot_json_filename):
    """
    In:
      plot_json_filename: json file name

    Returns:
      plot parameter dictionary
    """
    with open(plot_json_filename, 'r') as file:
        try:
            json_data = json.load(file)
        except json.decoder.JSONDecodeError as e:
            logging_array[0].error(f'Error while decoding json file "{plot_json_filename}".')
            logging_array[0].error(f'{e}')
            logging_array[0].error(f'File content:\n{file.read()}', level='error')
            exit(1)

    return json_data


def handle_plot_json(params, plot_json):
    """
    :param params: simulation parameters
    :param plot_json: plot_json file name
    :return: dictionary with plot parameters
    """

    # If a plot json is not specified search the params
    if plot_json is None:
        try:
            plot_json = params['plot_json']
        except KeyError:
            plot_json = None

    # Check to see if extension is given
    # If not we add to it
    if plot_json is not None and ".json" not in plot_json:
        plot_json += ".json"

    # Default folder in params unless stated otherwise
    if plot_json is not None and "/" not in plot_json:
        plot_json = 'params/' + plot_json

    # Check json file in directory:
    local_dir = Path(inspect.getfile(inspect.currentframe()))
    sim_dir = local_dir.absolute().parents[0]
    plot_json = os.path.join(sim_dir, plot_json)

    if not os.path.exists(plot_json):
        logging_array[0].warning("No JSON file detected plotting with default parameters")

    plot_params = get_plot_params(plot_json)
    return plot_params


####################### PLOTING FUNCTIONS


def plot_curves(data, axs, figure_index, plot_params, mappings = {}):
    """

        This function plots the curves in the assigned axis
        Also it looks for the final parameters inside the plot json, at the plot level

    :param axs: axs to plot the data in
    :param data: data given to be plot from pickle file
    :param figure_index: index of the current figure analysed
    :param plot_params: parameters for plotting
    :param mappings: species mappings set by the user
    """

    # Get the plot number from the list of plots
    try:
        plot_number = len(find_parameter(plot_params, 'plots', figure_index))
    except TypeError:
        # No figures plot only the default
        plot_number = 1

    # For all plots in the figure
    # Set all parameters and plot
    for plot_index in range(plot_number):

        # Get the species from plot parameters
        if find_parameter(plot_params, key='species_to_plot', index=(figure_index, plot_index)) is not None:
            species = find_parameter(plot_params, key='species_to_plot', index=(figure_index, plot_index))
        else:
            if mappings != {}:
                species = list(mappings.keys())
            else:
                species = list(data.keys())
                species.remove('Time')

        if find_parameter(plot_params, key='runs',  index=(figure_index, plot_index)) is not None:
            runs = find_parameter(plot_params, key='runs',  index=(figure_index, plot_index))
        else:
            runs = range(len(data[species[0]]['runs']))

        for spe in species:

            # Get the parameters assigned to the species, if not assign empty for None returns
            if find_parameter(plot_params, key=spe, index=(figure_index, plot_index)) is not None:
                species_characteristics = find_parameter(plot_params, key=spe, index=(figure_index, plot_index))
            else:
                species_characteristics = {}

            # Now we search the species caracteristics for parameters
            if find_parameter(species_characteristics, key='color') is not None:
                curve_color = find_parameter(species_characteristics, key='color')
            else:
                curve_color = color_cycle()

            if find_parameter(species_characteristics, key='linestyle') is not None:
                linestyle = find_parameter(species_characteristics, key='linestyle')
            else:
                linestyle = '-'

            for run in runs:
                axs.plot(data['Time'], data[spe]['runs'][run], color=curve_color,
                         linestyle = linestyle, linewidth=None)


def set_figure_characteristics(axis_list, plot_params):
    """

        Sets the figure parameters

    :param axis_list: list of all axis in the grid
    :param params: plot parameters received
    :return: nothing axis work by register, so their methods change them outside the function scope
    """

    # Loop through all axis
    for i, axs in enumerate(axis_list):

        if find_parameter(plot_params, 'xlim', i) is not None:
            figure_hash(i, axis_list).set_xlim(find_parameter(plot_params, 'xlim', i))

        if find_parameter(plot_params, 'ylim', i) is not None:
            figure_hash(i, axis_list).set_ylim(find_parameter(plot_params, 'ylim'))

        if find_parameter(plot_params, 'logscale', i) is not None and 'X' in find_parameter(plot_params, 'logscale', i):
            figure_hash(i, axis_list).set_xscale('log')

        if find_parameter(plot_params, 'logscale', i) is not None and 'Y' in find_parameter(plot_params, 'logscale', i):
            figure_hash(i, axis_list).set_yscale('log')

        if find_parameter(plot_params, 'title', i) is not None:
            figure_hash(i, axis_list).set_title(find_parameter(plot_params, 'title', i))

        if find_parameter(plot_params, 'xlabel', i) is not None:
            figure_hash(i, axis_list).set_xlabel(find_parameter(plot_params, 'xlabel', i))

        if find_parameter(plot_params, 'ylabel', i) is not None:
            figure_hash(i, axis_list).set_ylabel(find_parameter(plot_params, 'ylabel', i))

        if find_parameter(plot_params, 'annotations', i) is not None:
            for annotations in find_parameter(plot_params, 'annotations', i):
                figure_hash(i, axis_list).text(annotations[0], annotations[1], annotations[2])


def set_global_parameters(fig, plot_params):

    if find_parameter(plot_params, 'pad') is not None:
        fig.tight_layout(pad=find_parameter(plot_params, 'pad'))

    if find_parameter(plot_params, 'figsize') is not None:
        fig.set_size_inches(plot_params['figsize'][0], plot_params['figsize'][1])

    if find_parameter(plot_params, 'dpi') is not None:
        fig.set_dpi(plot_params['dpi'])


def plot_data(file_name, plot_json=None):
    """
        This function plots the simulation results according to the specifications

    :param file_name: pickle file resulting from simulation name
    :param module_logger: logger from simulation to store plotting info
    :param plot_json: plot_json file with plotting instructions
    """

    # Extract simulation parameters and data
    simulation_data = pkl.load(open(file_name, "rb"))
    params = simulation_data['params']
    data = simulation_data['data']
    mappings = simulation_data['mappings']

    # Extract plot parameter dictionary
    plot_params = handle_plot_json(params, plot_json)

    # Get the figure number from the list of figures
    try:
        figure_number = len(find_parameter(plot_params, 'figures'))
    except TypeError:
        # No figures plot only the default
        figure_number = 1

    # Create figures grid, max_lines is 2 as default and the value is defined in the function
    fig, axis_list = figure_hash_creation(figure_number, max_lines=find_parameter(plot_params, 'max_lines'))

    # Set parameters common to all figures, like padding and figure size.
    set_global_parameters(fig, plot_params)

    # Global call for setting global parameters
    set_figure_characteristics(axis_list, plot_params)

    # Now we plot
    for figure_index in range(figure_number):
        plot_curves(data, figure_hash(figure_index, axis_list), figure_index, plot_params, mappings)

    # Real default case
    plt.show()


if __name__ == "__main__":

    plot_data("/Users/fabriciocravo/Downloads/BacSim Copasi/output/AB_2CD/test/test.pkl", 'text2.json')

