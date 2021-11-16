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
    color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

    global __plot_color
    __plot_color = __plot_color + 1

    return color_list[__plot_color % len(color_list)]


# This function only works for the standart time we are currently using
def convert_to_time_unit(time_data, unit):
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
def figure_hash(current_figure, axis_list, total_figure_number):

    max_lines = len(axis_list)

    # Correction for index purposes
    current_figure = current_figure - 1

    col = math.floor(current_figure / max_lines)
    if total_figure_number <= max_lines:
        return axis_list[int(current_figure % max_lines)]
    else:
        return axis_list[int(current_figure % max_lines), int(col)]


# Hash to convert total figure number into grid
def figure_hash_creation(total_figure_number, max_lines=2):

    if total_figure_number <= max_lines:
        fig, axs = plt.subplots(total_figure_number)
    else:
        column_number = int((total_figure_number + total_figure_number % max_lines)/max_lines)
        fig, axs = plt.subplots(max_lines, column_number)

    if total_figure_number == 1:
        axs = [axs]

    return fig, axs


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
    :param params: simulation parameters from text file to check if a plot_json was specified
    :param plot_json: plot_json file name
    :param module_logger: module logger to keep track of what is happening
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


def plot_curve(data, axs, species, curve_parameters={}):

    try:
        species = curve_parameters['species']
    except KeyError:
        if species is None:
            logging_array[0].error('Must have a species as an argument or in json file')
            exit(1)

    try:
        runs = curve_parameters['runs']
    except KeyError:
        runs = range(len(data[species[0]]['runs']))

    for specie in species:

        try:
            curve_color = curve_parameters[specie + '_color']
        except KeyError:
            curve_color = color_cycle()

        for run in runs:
            axs.plot(data['Time'], data[specie]['runs'][run], color=curve_color)


def set_figure_characteristics(data, axis_list, figure_parameters, total_figure_number, figure_list, species=None, global_call=False):

    try:
        for figure in figure_list:
            figure_hash(figure, axis_list, total_figure_number).set_xlim(figure_parameters['xlim'])
    except KeyError:
        pass

    try:
        for figure in figure_list:
            figure_hash(figure, axis_list, total_figure_number).set_xlim(figure_parameters['xlim'])
    except KeyError:
        pass

    # Count the plots
    plot_number = 1
    while True:
        try:
            assert figure_parameters['plot_' + str(plot_number)]
        except KeyError:
            plot_number = plot_number - 1
            break
        plot_number = plot_number + 1

    # For global call only set axis parameters
    if global_call and plot_number == 0:
        return 0

    if plot_number != 0:
        for plot in range(plot_number):

            curve_parameters = figure_parameters['plot_' + str(plot + 1)]

            for figure in figure_list:
                axs = figure_hash(figure, axis_list, total_figure_number)
                plot_curve(data, axs, species, curve_parameters)
    else:
        plot_curve(data,  axis_list[0], species)


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

    # This is how we deal with multiple figures
    # First we count them
    figure_number = 1
    while True:
        try:
            assert plot_params['figure_' + str(figure_number)]
        except KeyError:
            figure_number = figure_number - 1
            break
        figure_number = figure_number + 1

    # We need at least one figure for default plot
    if figure_number == 0:
        figure_number = figure_number + 1

    # TODO We add one since index != figure number
    try:
        fig, axs = figure_hash_creation(figure_number, plot_params['max_lines'])
    except KeyError:
        fig, axs = figure_hash_creation(figure_number)

    # Set properties for all graphs
    if mappings != {}:
        species = list(mappings.keys())
    else:
        species = list(data.keys())
        species.remove('Time')

    # Global call for setting global parameters
    figure_list = list(range(figure_number))
    set_figure_characteristics(data, axs, plot_params,
                               total_figure_number=figure_number, figure_list=figure_list,
                               species=species, global_call=True)

    # if no plots were asked for we plot the default
    if 'figure_1' not in plot_params and 'plot_1' not in plot_params:
        plot_curve(data, axs[0], species)
    elif 'figure_1' not in plot_params:
        pass
    else:
        for i in range(figure_number):
            set_figure_characteristics(data, axs,
                                       plot_params['figure_' + str(i + 1)],
                                       species=species, figure_list=[i],
                                       total_figure_number=figure_number)

    # Real default case
    plt.show()


if __name__ == "__main__":

    plot_data("/Users/fabriciocravo/Downloads/BacSim Copasi/output/AB_2CD/test/test.pkl", 'text2.json')

