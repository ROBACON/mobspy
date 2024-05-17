import json
import os.path
import math
import mobspy.plot_scripts.statistics_calculations as spd
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.plot_scripts.process_plot_data as ppd
import numpy as np
import pickle as pkl
from pathlib import Path
import inspect
import matplotlib.pyplot as plt
import mobspy.modules.unit_handler as uh
from pint import Quantity


####################### PRACTICAL FUNCTIONS
class Color_cycle:
    """
        This class is responsible for cycling through the different colors for different curves

        :param index: (int) Current position in the color cycle
        :param color_list: (list) = list of available colors
    """

    def __init__(self):
        self.index = 0
        self.color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    def __call__(self, n):
        """
            Call updates the index and returns the next color in the list
            Normally n is set to 1

            :param n: (int) number of positions to skip
            :return: color (str) = color name for pyplot
        """
        self.index = (self.index + n) % len(self.color_list)
        return self.color_list[self.index]


def find_species_time_series(spe, data):
    """
        There can be different time-series in MobsPy data (even experimental data, as long as it is in MobPy format)
        This function finds all the time-series the species is present in and returns it for looping thorugh all of
        them
        This function is implemented to allow for the comparison of models with experimental data or other models

        :param spe: (str) Species name
        :pram  data: (dict) Data in MobsPy format
    """
    for time_series in data:
        if spe in time_series:
            yield time_series


def get_total_figure_number(axis_matrix):
    """
        Gets the number of figures from an axis_matrix generated by pyplot. Used by the hash to place configs in the
        correct part of the axis_matrix

        :param axis_matrix: (1D or 2D numpy array) Array returned by pyplot once multiple figures are introduced into
        a single plot

        :return: total_figure_number (int) Total number of figures
    """
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

        :param current_figure: (int) linear number of the figure
        :param axis_matrix: (1D or 2D numpy array) a list with all the created axis on the multiple figure subplot

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

        This is hash used to create the figure grid automatically according to the number of figures desired by the
        user and the maximum number of lines

        For instance 4 figures with 2 as max_lines creates a 2x2 figure grid automatically

        :param total_figure_number: (int) Number of total figures to create
        :param max_lines: (int) Maximum number of lines in the grid

        :return: fig and axs = Figures grid will all the respective axis
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

        The priority for parameter search is plots => figures => global, with plot overriding others and so on

        If a parameter is defined globally it will be applied to all figures, if it defined inside a figure element
        it will only apply to that figure, if it is defined in a plot element it will only apply to that curve

        Check the readme or the tutorials for more details on the plotting structure. It is simple and versatile

        :param params: (dict) Plot parameters from python dictionary (after json conversion)
        :param key: (str) Key necessary to access the parameters
        :param index: (int) None for global search, one index for figure search, and two for figure curve search

        :return: The parameter if found, and None if nothing is found
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
    # If two indexes are given, look inside the plot, than figure, than global
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


def annotation_handling(axs, figure_index, plot_index, plot_params):
    print()

    if find_parameter(plot_params, key='annotations', index=(figure_index, plot_index)) is not None:
        annotations = find_parameter(plot_params, key='annotations', index=(figure_index, plot_index))

        if type(annotations) != list:
            simlog.warning("On plotting annotations: Annotations must a be list with dictionaries as elements")
            return 0

        for annotation_dict in annotations:
            argument_dict = {}

            if "text" not in annotation_dict:
                text = "Default Annotation"
            else:
                text = annotation_dict["text"]

            if 'coordinates' not in annotation_dict:
                coordinates = (0, 0)
            else:
                coordinates = annotation_dict['coordinates']

            arg_list = ["textcoords", "xytext", "ha", "fontsize"]
            for arg in arg_list:
                if arg in annotation_dict:
                    argument_dict[arg] = annotation_dict[arg]

            axs.annotate(text, coordinates, **argument_dict)
        else:
            return 0


####################### PLOTING FUNCTIONS
def plot_curves(data, axs, figure_index, plot_params):
    """
        This function plots the programmed curves in the assigned figure

        :param axs: (pyplot axe) axs to plot the data in
        :param data: (dict) data given in MobsPy format results['data']
        :param figure_index: (int) index of the current figure to plot curves in
        :param plot_params: parameters for plotting
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

        annotation_handling(axs, figure_index, plot_index, plot_params)

        # Get the species from plot parameters
        if find_parameter(plot_params, key='species_to_plot', index=(figure_index, plot_index)) is not None:
            species = find_parameter(plot_params, key='species_to_plot', index=(figure_index, plot_index))
        else:
            simlog.error('No species found for plotting in one of the curves or figures')

        species = sorted([spe for spe in species])

        # Get the time series to plot
        if find_parameter(plot_params, key='time_series', index=(figure_index, plot_index)) is not None:
            time_series = find_parameter(plot_params, key='time_series', index=(figure_index, plot_index))
            if type(time_series) == int:
                time_series = [time_series]
        else:
            time_series = list(range(len(data)))

        if find_parameter(plot_params, key='time_filter', index=(figure_index, plot_index)) is not None:
            low, high = find_parameter(plot_params, key='time_filter', index=(figure_index, plot_index))
        else:
            low, high = None, None

        if find_parameter(plot_params, key='y_filter', index=(figure_index, plot_index)) is not None:
            low_y, high_y = find_parameter(plot_params, key='y_filter', index=(figure_index, plot_index))
        else:
            low_y, high_y = None, None

        if find_parameter(plot_params, key='x_from', index=(figure_index, plot_index)) is not None:
            x_start, x_finish = find_parameter(plot_params, key='x_from', index=(figure_index, plot_index))
        else:
            x_start, x_finish = None, None

        if find_parameter(plot_params, key='y_from', index=(figure_index, plot_index)) is not None:
            y_start, y_finish = find_parameter(plot_params, key='y_from', index=(figure_index, plot_index))
        else:
            y_start, y_finish = None, None

        for spe in species:

            # Get the parameters assigned to the species, if not assign empty for None returns
            if find_parameter(plot_params, key=spe, index=(figure_index, plot_index)) is not None:
                species_characteristics = find_parameter(plot_params, key=spe, index=(figure_index, plot_index))
            else:
                species_characteristics = {}

            # Now we search the species characteristics for parameters
            if find_parameter(species_characteristics, key='color') is not None:
                curve_color = find_parameter(species_characteristics, key='color')
            else:
                curve_color = None

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

            if find_parameter(species_characteristics, key='ylabel') is not None:
                y_label = find_parameter(species_characteristics, key='ylabel')
                if find_parameter(plot_params, 'ylabel_fontsize', figure_index) is not None:
                    axs.set_ylabel(y_label, fontsize=plot_params['ylabel_fontsize'])
                else:
                    axs.set_ylabel(y_label)

            for ts in time_series:
                ts_time = data['Time'][ts]
                if '$' not in spe:
                    ts_data = data[spe][ts]
                else:
                    ts_data = data[ts][spe]

                if low is not None and high is not None:
                    ts_time, ts_data = ppd.time_filter_operation(low, high, ts_time, ts_data)

                if low_y is not None and high_y is not None:
                    ts_time, ts_data = ppd.y_filter_operation(low_y, high_y, ts_time, ts_data)

                if x_start is not None and x_finish is not None:
                    ref_point = ts_data[-1]
                    axs.plot([x_start, x_finish], [ref_point,  ref_point], alpha=0)

                if y_start is not None and y_finish is not None:
                    ref_point = ts_time[-1]
                    axs.plot([ref_point, ref_point], [y_start, y_finish], alpha=0)

                try:
                    if find_parameter(plot_params, key='fill_between', index=(figure_index, plot_index)) is not None \
                            and find_parameter(plot_params, key='fill_between', index=(figure_index, plot_index)):
                        try:
                            axs.fill_between(ts_time, ts_data[0], ts_data[1], color=curve_color,
                                             label=label)
                        except IndexError:
                            simlog.error('Fill_between must only have two or less runs referring to it')
                    else:
                        if curve_color is not None:
                            axs.plot(ts_time, ts_data, color=curve_color,
                                     linestyle=linestyle, linewidth=linewidth, label=label)
                        else:
                            axs.plot(ts_time, ts_data,
                                     linestyle=linestyle, linewidth=linewidth, label=label)
                        label = None
                except KeyError:
                    pass

        if find_parameter(plot_params, key='vertical_lines', index=(figure_index, plot_index)) is not None:
            unit_x = find_parameter(plot_params, key='unit_x', index=(figure_index, plot_index))

            position = find_parameter(plot_params, key='vertical_lines', index=(figure_index, plot_index))
            for p in position:
                if unit_x is not None and isinstance(p, Quantity):
                    new_p = uh.time_convert_to_other_unit(p, unit_x)
                else:
                    new_p = p
                axs.axvline(x=new_p, color='gray', linestyle='--')

    if legend_flag:
        if find_parameter(plot_params, key='prop', index=figure_index) is not None:
            prop = find_parameter(plot_params, key='prop', index=figure_index)
        else:
            prop = {'size': 10}

        if find_parameter(plot_params, key='frameon', index=figure_index) is not None:
            frameon = find_parameter(plot_params, key='frameon', index=figure_index)
        else:
            frameon = True
        axs.legend(frameon=frameon, prop=prop)


def set_figure_characteristics(axis_matrix, plot_params):
    """
        Sets the characteristics for each figure

        :param axis_matrix: (axis from pyplot subplot) Array of all axis in the grid
        :param plot_params: (dict) Plot parameters received
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
            if find_parameter(plot_params, 'title_fontsize', i) is not None:
                figure_hash(i, axis_matrix).set_title(find_parameter(plot_params, 'title', i),
                                                      fontsize=plot_params['title_fontsize'])
            else:
                figure_hash(i, axis_matrix).set_title(find_parameter(plot_params, 'title', i))

        if find_parameter(plot_params, 'xlabel', i) is not None:
            if find_parameter(plot_params, 'xlabel_fontsize', i) is not None:
                figure_hash(i, axis_matrix).set_xlabel(find_parameter(plot_params, 'xlabel', i),
                                                       fontsize=plot_params['xlabel_fontsize'])
            else:
                figure_hash(i, axis_matrix).set_xlabel(find_parameter(plot_params, 'xlabel', i))

        if find_parameter(plot_params, 'ylabel', i) is not None:
            if find_parameter(plot_params, 'ylabel_fontsize', i) is not None:
                figure_hash(i, axis_matrix).set_ylabel(find_parameter(plot_params, 'ylabel', i),
                                                       fontsize=plot_params['ylabel_fontsize'])
            else:
                figure_hash(i, axis_matrix).set_ylabel(find_parameter(plot_params, 'ylabel', i))


def set_global_parameters(fig, plot_params):
    """
        Sets the characteristics the plot window

        :param fig: (fig from pyplot subplot) Array of all axis in the grid
        :param plot_params: (dict) Plot parameters received
    """
    if find_parameter(plot_params, 'pad') is not None:
        fig.tight_layout(pad=find_parameter(plot_params, 'pad'))

    if find_parameter(plot_params, 'figsize') is not None:
        fig.set_size_inches(plot_params['figsize'][0], plot_params['figsize'][1])

    if find_parameter(plot_params, 'dpi') is not None:
        fig.set_dpi(plot_params['dpi'])

    if find_parameter(plot_params, 'suptitle') is not None:
        if find_parameter(plot_params, 'suptitle_fontsize') is not None:
            fig.suptitle(plot_params['suptitle'], fontsize=plot_params['suptitle_fontsize'])
        else:
            fig.suptitle(plot_params['suptitle'])


def plot_data(data, plot_params, return_fig_object=False):
    """
        This function plots the simulation results according to the specifications

        :param data: (dict) Data in MobsPy format
        :param plot_params: (dict) Plot parameters received
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

    tight_layout = find_parameter(plot_params, 'tight_layout')
    if tight_layout is not None:
        fig.tight_layout()

    # Now we plot
    for figure_index in range(figure_number):
        plot_curves(data, figure_hash(figure_index, axis_matrix), figure_index, plot_params)

    # Real default case
    if return_fig_object:
        return fig, axis_matrix

    save_to = find_parameter(plot_params, 'save_to')
    if save_to is None:
        plt.show()
    else:
        plt.savefig(save_to)

    return None


if __name__ == "__main__":
    pass
