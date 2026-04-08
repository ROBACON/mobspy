"""Render multi-panel hierarchical plots for parametric sweep results."""

from __future__ import annotations

import contextlib
import math
from typing import TYPE_CHECKING, Any

from mobspy.exceptions import ValidationError
from mobspy.mobspy_logging import get_logger

_logger = get_logger(__name__)
import matplotlib.pyplot as plt  # noqa: E402  # deferred import after lazy setup
import numpy as np  # noqa: E402  # deferred import after lazy setup
from pint import Quantity  # noqa: E402  # deferred import after lazy setup

import mobspy.plot.process_plot_data as ppd  # noqa: E402  # deferred import
import mobspy.units.handler as uh  # noqa: E402  # deferred import after lazy setup

if TYPE_CHECKING:
    from collections.abc import Generator


####################### PRACTICAL FUNCTIONS
class Color_cycle:  # noqa: N801  # legacy DSL public API name
    """
    This class is responsible for cycling through the
    different colors for different curves.

    Args:
        index: Current position in the color cycle.
        color_list: List of available colors.
    """

    def __init__(self) -> None:
        self.index = 0
        self.color_list = ["b", "g", "r", "c", "m", "y", "k"]

    def __call__(self, n: int) -> str:
        """
        Call updates the index and returns the next color in the list
        Normally n is set to 1

        Args:
            n: Number of positions to skip.


        Returns:
            Color name for pyplot.
        """
        self.index = (self.index + n) % len(self.color_list)
        return self.color_list[self.index]


def find_species_time_series(spe: str, data: Any) -> Generator[Any, None, None]:
    """
    There can be different time-series in MobsPy data
    (even experimental data, as long as it is in MobPy
    format). This function finds all the time-series the
    species is present in and returns it for looping
    through all of them.
    This function is implemented to allow for the comparison
    of models with experimental data or other models.

    Args:
        spe: Species name.
        data: Data in MobsPy format.
    """
    for time_series in data:
        if spe in time_series:
            yield time_series


def get_total_figure_number(axis_matrix: np.ndarray[Any, Any]) -> int:
    """
    Gets the number of figures from an axis_matrix generated
    by pyplot. Used by the hash to place configs in the
    correct part of the axis_matrix.

    Args:
        axis_matrix: Array returned by pyplot once multiple figures are introduced into
            a single plot.


    Returns:
        Total number of figures.
    """
    # Get the number of figures, using the axis_matrix
    try:
        total_figure_number = axis_matrix.shape[0] * axis_matrix.shape[1]
    except IndexError:
        total_figure_number = axis_matrix.shape[0]
    return total_figure_number


# Hash for converting linear figure number into index
def figure_hash(current_figure: int, axis_matrix: np.ndarray[Any, Any]) -> Any:
    """
    This function allows one to access the figure grid with a linear input
    For instance one can access a 2x2 grid using 0, 1, 2, 3
    0 becomes 0,0
    1 becomes 1,0
    2 becomes 0,1
    3 becomes 1,1

    Args:
        current_figure: Linear number of the figure.
        axis_matrix: A list with all the created axis on the multiple figure subplot.


    Returns:
        The correct axis based on the number provided.
    """

    # Get the number of lines
    max_lines = len(axis_matrix)
    if max_lines == 0:
        raise ValueError("axis_matrix is empty; cannot resolve figure hash")
    total_figure_number = get_total_figure_number(axis_matrix)

    col = math.floor(current_figure / max_lines)

    if total_figure_number <= max_lines:
        return axis_matrix[int(current_figure % max_lines)]
    return axis_matrix[int(current_figure % max_lines), int(col)]


# Hash to convert total figure number into grid
def figure_hash_creation(
    total_figure_number: int, max_lines: int | None = None
) -> tuple[Any, np.ndarray[Any, Any]]:
    """

    This is hash used to create the figure grid automatically
    according to the number of figures desired by the user
    and the maximum number of lines

    For instance 4 figures with 2 as max_lines creates a 2x2 figure grid automatically

    Args:
        total_figure_number: Number of total figures to create.
        max_lines: Maximum number of lines in the grid.


    Returns:
        Fig and axs = Figures grid will all the respective axis.
    """

    # Default value if nothing is set up
    if max_lines is None:
        max_lines = 2

    if total_figure_number <= max_lines:
        fig, axs = plt.subplots(total_figure_number)
    else:
        column_number = int(
            (total_figure_number + total_figure_number % max_lines) / max_lines
        )
        fig, axs = plt.subplots(max_lines, column_number)

    if total_figure_number == 1:
        axs = np.array([axs])

    return fig, axs


def find_parameter(  # noqa: PLR0911  # complex DSL method with multiple return paths
    params: dict[str, Any], key: str, index: int | tuple[int, ...] | None = None
) -> Any:
    """

    This is the heart of the plotting structure, this function
    allows one to simply set multiple characteristics for
    only one figure.

    The priority for parameter search is
    plots => figures => global, with plot overriding others
    and so on.

    If a parameter is defined globally it will be applied
    to all figures, if it defined inside a figure element
    it will only apply to that figure, if it is defined in
    a plot element it will only apply to that curve.

    Check the readme or the tutorials for more details on
    the plotting structure. It is simple and versatile.

    Args:
        params: Plot parameters from python dictionary (after json conversion).
        key: Key necessary to access the parameters.
        index: None for global search, one index for figure search, and two for figure
            curve search.


    Returns:
        The parameter if found, and None if not found.
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
    elif isinstance(index, int):
        try:
            return params["figures"][index][key]
        except (KeyError, IndexError):
            try:
                return params[key]
            except (KeyError, IndexError):
                return None
    # If two indexes are given, look inside the plot, than figure, than global
    elif isinstance(index, tuple):
        try:
            return params["figures"][index[0]]["plots"][index[1]][key]

        except (KeyError, IndexError):
            try:
                return params["figures"][index[0]][key]

            except (KeyError, IndexError):
                try:
                    return params[key]
                except (KeyError, IndexError):
                    return None
    return None


def annotation_handling(
    axs: Any, figure_index: int, plot_index: int, plot_params: dict[str, Any]
) -> int | None:
    """Apply user-defined annotations to a subplot axis."""
    if (
        find_parameter(plot_params, key="annotations", index=(figure_index, plot_index))
        is not None
    ):
        annotations = find_parameter(
            plot_params, key="annotations", index=(figure_index, plot_index)
        )

        if not isinstance(annotations, list):
            _logger.warning(
                "On plotting annotations: Annotations must "
                "be a list with dictionaries as elements"
            )
            return 0

        for annotation_dict in annotations:
            argument_dict: dict[str, Any] = {}

            if "text" not in annotation_dict:
                text = "Default Annotation"
            else:
                text = annotation_dict["text"]

            if "coordinates" not in annotation_dict:
                coordinates = (0, 0)
            else:
                coordinates = annotation_dict["coordinates"]

            arg_list = ["textcoords", "xytext", "ha", "fontsize"]
            for arg in arg_list:
                if arg in annotation_dict:
                    argument_dict[arg] = annotation_dict[arg]

            axs.annotate(text, coordinates, **argument_dict)
        return 0
    return None


####################### PLOTTING FUNCTIONS


def _get_param(
    params: dict[str, Any],
    key: str,
    index: int | tuple[int, ...] | None = None,
    default: Any = None,
) -> Any:
    """Retrieve a plot parameter with a default fallback."""
    val = find_parameter(params, key=key, index=index)
    return val if val is not None else default


def _get_plot_filters(
    plot_params: dict[str, Any],
    figure_index: int,
    plot_index: int,
) -> tuple[Any, Any, Any, Any, Any, Any, Any, Any]:
    """Extract time/y/x/y range filters for a single plot."""
    idx = (figure_index, plot_index)
    time_filter = find_parameter(plot_params, key="time_filter", index=idx)
    low, high = time_filter if time_filter is not None else (None, None)

    y_filter = find_parameter(plot_params, key="y_filter", index=idx)
    low_y, high_y = y_filter if y_filter is not None else (None, None)

    x_from = find_parameter(plot_params, key="x_from", index=idx)
    x_start, x_finish = x_from if x_from is not None else (None, None)

    y_from = find_parameter(plot_params, key="y_from", index=idx)
    y_start, y_finish = y_from if y_from is not None else (None, None)

    return low, high, low_y, high_y, x_start, x_finish, y_start, y_finish


def _get_species_style(
    species_characteristics: dict[str, Any],
) -> tuple[Any, str, Any, Any]:
    """Extract curve style attributes from species characteristics."""
    curve_color = _get_param(species_characteristics, "color")
    linestyle = _get_param(species_characteristics, "linestyle", default="-")
    linewidth = _get_param(species_characteristics, "linewidth")
    label = _get_param(species_characteristics, "label")
    return curve_color, linestyle, linewidth, label


def _plot_single_curve(  # noqa: PLR0913  # complex function signature
    axs: Any,
    ts_time: Any,
    ts_data: Any,
    curve_color: Any,
    linestyle: str,
    linewidth: Any,
    label: Any,
    fill_between: bool,
) -> Any:
    """Plot a single time series curve or fill_between region. Returns updated label."""
    if fill_between:
        try:
            axs.fill_between(
                ts_time,
                ts_data[0],
                ts_data[1],
                color=curve_color,
                label=label,
            )
        except IndexError as e:
            raise ValidationError(
                "Fill_between must only have two or less runs referring to it"
            ) from e
    else:
        plot_kwargs: dict[str, Any] = {
            "linestyle": linestyle,
            "linewidth": linewidth,
            "label": label,
        }
        if curve_color is not None:
            plot_kwargs["color"] = curve_color
        axs.plot(ts_time, ts_data, **plot_kwargs)
        label = None
    return label


def _apply_filters_and_ranges(  # noqa: PLR0913  # complex function signature
    axs: Any,
    ts_time: Any,
    ts_data: Any,
    low: Any,
    high: Any,
    low_y: Any,
    high_y: Any,
    x_start: Any,
    x_finish: Any,
    y_start: Any,
    y_finish: Any,
) -> tuple[Any, Any]:
    """Apply time/y filters and invisible range markers."""
    if low is not None and high is not None:
        ts_time, ts_data = ppd.time_filter_operation(low, high, ts_time, ts_data)
    if low_y is not None and high_y is not None:
        ts_time, ts_data = ppd.y_filter_operation(low_y, high_y, ts_time, ts_data)
    if x_start is not None and x_finish is not None:
        axs.plot([x_start, x_finish], [ts_data[-1], ts_data[-1]], alpha=0)
    if y_start is not None and y_finish is not None:
        axs.plot([ts_time[-1], ts_time[-1]], [y_start, y_finish], alpha=0)
    return ts_time, ts_data


def _plot_species_curves(  # noqa: PLR0913  # complex function signature
    data: Any,
    axs: Any,
    species: list[Any],
    time_series: list[int],
    plot_params: dict[str, Any],
    figure_index: int,
    plot_index: int,
    filters: tuple[Any, ...],
    fill_between: bool,
) -> bool:
    """Plot all species curves for one plot. Returns True if any label was set."""
    idx = (figure_index, plot_index)
    low, high, low_y, high_y, x_start, x_finish, y_start, y_finish = filters
    legend_flag = False

    for spe in species:
        spe_chars = _get_param(plot_params, spe, idx, default={})
        curve_color, linestyle, linewidth, label = _get_species_style(spe_chars)
        if label is not None:
            legend_flag = True

        ylabel_val = _get_param(spe_chars, "ylabel")
        if ylabel_val is not None:
            fontsize = _get_param(plot_params, "ylabel_fontsize", figure_index)
            _apply_label_with_fontsize(axs, "set_ylabel", ylabel_val, fontsize)

        for ts in time_series:
            ts_time = data["Time"][ts]
            ts_data = data[spe][ts] if "$" not in spe else data[ts][spe]
            ts_time, ts_data = _apply_filters_and_ranges(
                axs,
                ts_time,
                ts_data,
                low,
                high,
                low_y,
                high_y,
                x_start,
                x_finish,
                y_start,
                y_finish,
            )
            with contextlib.suppress(KeyError):
                label = _plot_single_curve(
                    axs,
                    ts_time,
                    ts_data,
                    curve_color,
                    linestyle,
                    linewidth,
                    label,
                    fill_between,
                )

    return legend_flag


def _draw_vertical_lines(
    axs: Any,
    plot_params: dict[str, Any],
    idx: tuple[int, int],
) -> None:
    """Draw vertical guide lines if configured."""
    vlines = _get_param(plot_params, "vertical_lines", idx)
    if vlines is not None:
        unit_x = _get_param(plot_params, "unit_x", idx)
        for p in vlines:
            new_p = (
                uh.time_convert_to_other_unit(p, unit_x)
                if unit_x is not None and isinstance(p, Quantity)
                else p
            )
            axs.axvline(x=new_p, color="gray", linestyle="--")


def plot_curves(
    data: Any, axs: Any, figure_index: int, plot_params: dict[str, Any]
) -> None:
    """
    This function plots the programmed curves in the assigned figure

    Args:
        axs: Axs to plot the data in.
        data: Data given in MobsPy format results['data'].
        figure_index: Index of the current figure to plot curves in.
        plot_params: Parameters for plotting.
    """
    try:
        plot_number = len(find_parameter(plot_params, "plots", figure_index))
    except TypeError:
        plot_number = 1
    if plot_number == 0:
        plot_number = 1

    legend_flag = False
    for plot_index in range(plot_number):
        annotation_handling(axs, figure_index, plot_index, plot_params)
        idx = (figure_index, plot_index)

        species = find_parameter(plot_params, key="species_to_plot", index=idx)
        if species is None:
            raise ValidationError(
                "No species found for plotting in one of the curves or figures"
            )
        species = sorted(species)

        time_series = _get_param(plot_params, "time_series", idx)
        if time_series is None:
            time_series = list(range(len(data)))
        elif isinstance(time_series, int):
            time_series = [time_series]

        filters = _get_plot_filters(plot_params, figure_index, plot_index)
        fill_between = bool(_get_param(plot_params, "fill_between", idx))

        if _plot_species_curves(
            data,
            axs,
            species,
            time_series,
            plot_params,
            figure_index,
            plot_index,
            filters,
            fill_between,
        ):
            legend_flag = True

        _draw_vertical_lines(axs, plot_params, idx)

    if legend_flag:
        prop = _get_param(plot_params, "prop", figure_index, default={"size": 10})
        frameon = _get_param(plot_params, "frameon", figure_index, default=True)
        axs.legend(frameon=frameon, prop=prop)


def _apply_label_with_fontsize(
    ax: Any,
    setter_name: str,
    value: Any,
    fontsize: Any,
) -> None:
    """Call an axis setter method, optionally passing fontsize."""
    setter = getattr(ax, setter_name)
    if fontsize is not None:
        setter(value, fontsize=fontsize)
    else:
        setter(value)


def set_figure_characteristics(
    axis_matrix: np.ndarray[Any, Any], plot_params: dict[str, Any]
) -> None:
    """
    Sets the characteristics for each figure

    Args:
        axis_matrix: Array of all axis in the grid.
        plot_params: Plot parameters received.
    """
    total_figure_number = get_total_figure_number(axis_matrix)
    for i in range(total_figure_number):
        ax = figure_hash(i, axis_matrix)

        xlim = _get_param(plot_params, "xlim", i)
        if xlim is not None:
            ax.set_xlim(xlim)

        ylim = _get_param(plot_params, "ylim", i)
        if ylim is not None:
            ax.set_ylim(ylim)

        logscale = _get_param(plot_params, "logscale", i)
        if logscale is not None:
            if "X" in logscale:
                ax.set_xscale("log")
            if "Y" in logscale:
                ax.set_yscale("log")

        label_configs: list[tuple[str, str, str]] = [
            ("title", "title_fontsize", "set_title"),
            ("xlabel", "xlabel_fontsize", "set_xlabel"),
            ("ylabel", "ylabel_fontsize", "set_ylabel"),
        ]
        for key, fontsize_key, setter_name in label_configs:
            val = _get_param(plot_params, key, i)
            if val is not None:
                fs = _get_param(plot_params, fontsize_key, i)
                _apply_label_with_fontsize(ax, setter_name, val, fs)


def set_global_parameters(fig: Any, plot_params: dict[str, Any]) -> None:
    """
    Sets the characteristics the plot window

    Args:
        fig: Array of all axis in the grid.
        plot_params: Plot parameters received.
    """
    if find_parameter(plot_params, "pad") is not None:
        fig.tight_layout(pad=find_parameter(plot_params, "pad"))

    if find_parameter(plot_params, "figsize") is not None:
        fig.set_size_inches(plot_params["figsize"][0], plot_params["figsize"][1])

    if find_parameter(plot_params, "dpi") is not None:
        fig.set_dpi(plot_params["dpi"])

    if find_parameter(plot_params, "suptitle") is not None:
        if find_parameter(plot_params, "suptitle_fontsize") is not None:
            fig.suptitle(
                plot_params["suptitle"], fontsize=plot_params["suptitle_fontsize"]
            )
        else:
            fig.suptitle(plot_params["suptitle"])


def plot_data(
    data: Any, plot_params: dict[str, Any], return_fig_object: bool = False
) -> Any:
    """
    This function plots the simulation results according to the specifications

    Args:
        data: Data in MobsPy format.
        plot_params: Plot parameters received.
    """
    # Get the figure number from the list of figures
    # Add it to parameters
    try:
        figure_number = len(find_parameter(plot_params, "figures"))
    except TypeError:
        # No figures plot only the default
        figure_number = 1
    if figure_number == 0:
        figure_number = 1

    # Create figures grid, max_lines is 2 as default
    # and the value is defined in the function
    fig, axis_matrix = figure_hash_creation(
        figure_number, max_lines=find_parameter(plot_params, "max_lines")
    )

    # Set parameters common to all figures, like padding and figure size.
    set_global_parameters(fig, plot_params)

    # Global call for setting global parameters
    set_figure_characteristics(axis_matrix, plot_params)

    tight_layout = find_parameter(plot_params, "tight_layout")
    if tight_layout is not None:
        fig.tight_layout()

    # Now we plot
    for figure_index in range(figure_number):
        plot_curves(
            data, figure_hash(figure_index, axis_matrix), figure_index, plot_params
        )

    # Real default case
    if return_fig_object:
        return fig, axis_matrix

    save_to = find_parameter(plot_params, "save_to")
    if save_to is None:
        plt.show()
    else:
        plt.savefig(save_to)

    return None


if __name__ == "__main__":
    pass
