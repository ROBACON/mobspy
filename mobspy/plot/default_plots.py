"""Generate matplotlib plots from simulation time-series data."""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from typing import Any

from pint import Quantity

import mobspy.plot.hierarchical_plot as hp
import mobspy.plot.process_plot_data as ppd
import mobspy.plot.statistics as sc
from mobspy.constants import AVERAGE_SUFFIX, DEVIATION_SUFFIX
from mobspy.exceptions import ValidationError

# Thresholds for legend label font size selection
_SMALL_LABEL_THRESHOLD = 5
_MEDIUM_LABEL_THRESHOLD = 10


def read_plot_json(plot_json_filename: str) -> dict[str, Any]:
    """
    This function converts a plot_json file into a dictionary

    Args:
        plot_json_filename: JSON file name.


    Returns:
        Converted JSON as dictionary.
    """
    with Path(plot_json_filename).open(encoding="utf-8") as file:
        try:
            json_data = json.load(file)
        except (json.JSONDecodeError, ValueError) as e:
            raise ValidationError(
                "The following error happened while "
                "decoding json file "
                f'"{plot_json_filename}":\n' + str(e)
            ) from e

    return json_data  # type: ignore[no-any-return]  # dynamic dispatch


def set_plot_units(new_plot_params: dict[str, Any]) -> None:
    """
    Sets the plot labels to the unit names by adding to the xlabel and ylabel

    Args:
        new_plot_params: Plot parameters after some changes.
    """
    if "xlabel" not in new_plot_params:
        new_plot_params["xlabel"] = "Time"

        if new_plot_params["unit_x"] is not None and not (
            new_plot_params.get("ignore_unit_label_x")
        ):
            if not isinstance(new_plot_params["unit_x"], Quantity):
                new_plot_params["xlabel"] += f" ({new_plot_params['unit_x']}s)"
            else:
                new_plot_params["xlabel"] += f" ({new_plot_params['unit_x'].units}s)"

    if "ylabel" not in new_plot_params:
        if new_plot_params["output_concentration"]:
            new_plot_params["ylabel"] = "Conc."
        else:
            new_plot_params["ylabel"] = "Counts"

        if new_plot_params["unit_y"] is not None and not new_plot_params.get(
            "ignore_unit_label_y"
        ):
            if not isinstance(new_plot_params["unit_y"], Quantity):
                new_plot_params["ylabel"] += f" ({new_plot_params['unit_y']})"
            else:
                new_plot_params["ylabel"] += f" ({new_plot_params['unit_y'].units})"


def stochastic_plot(
    species: set[str] | list[str],
    data: Any,
    plot_params: dict[str, Any],
) -> Any:
    """
    Design default stochastic plot using MobsPy plotting
    hierarchy. It then passes the parameters for plotting
    in the hierarchical plot module.

    Args:
        species: List of species names in MobsPy str format (queries are performed with
            the query plot data function).
        data: Data in MobsPy format Simulation.results['data'].
        plot_params: Dictionary with the plot parameters supplied by the user before the
            modifications from this function.
    """

    # Data Handling - Copy data object to not interfere with simulation data
    species, data = ppd.query_plot_data(species, data)
    ppd.check_plot_parameters(species, plot_params)

    data_to_plot = data

    try:
        new_plot_params = deepcopy(plot_params)
    except (TypeError, AttributeError):
        new_plot_params = plot_params
    set_plot_units(new_plot_params)

    # Plot config
    new_plot_params["frameon"] = False
    new_plot_params["figures"] = []
    new_plot_params["pad"] = 1.5
    color_cycler = hp.Color_cycle()
    for spe in species:
        # We define new 'mappings' with the resulting
        # runs for the statistics for the plot structure
        try:
            plots_for_spe_i: list[dict[str, Any]] = []
            plots_for_spe_i_sta: list[dict[str, Any]] = []

            processed_runs = sc.average_plus_standard_deviation(spe, data_to_plot)

            key_average = spe + AVERAGE_SUFFIX
            key_dev = spe + DEVIATION_SUFFIX

            data_to_plot[0][key_average] = processed_runs[0]
            data_to_plot[0][key_dev] = [processed_runs[1], processed_runs[2]]

            # We define the standard plot for the average and deviation
            plots_for_spe_i.append({"species_to_plot": [spe]})

            plots_for_spe_i_sta.append(
                {"species_to_plot": [key_average], "time_series": [0]}
            )
            plots_for_spe_i_sta.append(
                {"species_to_plot": [key_dev], "fill_between": True, "time_series": [0]}
            )

        except ValueError as e:
            raise ValidationError(f"{spe} species not found in data") from e
        new_plot_params["figures"].append(
            {"ylabel": spe + " " + new_plot_params["ylabel"], "plots": plots_for_spe_i}
        )
        new_plot_params["figures"].append(
            {
                "ylabel": spe + " " + new_plot_params["ylabel"],
                "plots": plots_for_spe_i_sta,
            }
        )

    # Setting species parameters
    for spe in species:
        key_average = spe + AVERAGE_SUFFIX
        key_dev = spe + DEVIATION_SUFFIX

        color = color_cycler(1)

        new_plot_params[key_dev] = {
            "color": (0.8, 0.8, 0.8),
            "linestyle": ":",
            "label": "std. dev",
        }
        if spe not in new_plot_params:
            new_plot_params[spe] = {"color": color, "label": spe}
            new_plot_params[key_average] = {
                "color": color,
                "linestyle": "-",
                "label": "mean",
            }
        else:
            new_plot_params[key_average] = {}

            for par in new_plot_params[spe]:
                new_plot_params[key_average][par] = new_plot_params[spe][par]
                new_plot_params[key_dev][par] = new_plot_params[spe][par]

            if "label" not in new_plot_params[spe]:
                new_plot_params[spe]["label"] = spe

            if "color" not in new_plot_params[spe]:
                new_plot_params[spe]["color"] = color
                new_plot_params[key_average]["color"] = color
                new_plot_params[key_average]["linestyle"] = "-"
                new_plot_params[key_average]["label"] = "mean"
            else:
                new_plot_params[key_average]["color"] = new_plot_params[spe]["color"]
                new_plot_params[key_average]["linestyle"] = "-"
                new_plot_params[key_average]["label"] = "mean"

    return hp.plot_data(data_to_plot, new_plot_params)


def deterministic_plot(
    species: set[str] | list[str],
    data: Any,
    plot_params: dict[str, Any],
) -> Any:
    """
    Design default deterministic plot using MobsPy
    plotting hierarchy. It then passes the parameters
    for plotting in the hierarchical plot module.

    Args:
        species: List of species names in MobsPy str format (queries are performed with
            the query plot data function).
        data: Data in MobsPy format Simulation.results['data'].
        plot_params: Dictionary with the plot parameters supplied by the user before the
            modifications from this function.
    """
    # Data Handling
    species, data = ppd.query_plot_data(species, data)
    ppd.check_plot_parameters(species, plot_params)

    try:
        new_plot_params = deepcopy(plot_params)
    except (TypeError, AttributeError):
        new_plot_params = plot_params
    set_plot_units(new_plot_params)

    #  Plot Config
    new_plot_params["frameon"] = False
    new_plot_params["species_to_plot"] = species
    color_cycler = hp.Color_cycle()
    for spe in species:
        if spe not in new_plot_params:
            new_plot_params[spe] = {"label": spe, "color": color_cycler(1)}
        else:
            if "label" not in new_plot_params[spe]:
                new_plot_params[spe]["label"] = spe
            if "color" not in new_plot_params[spe]:
                new_plot_params[spe]["color"] = color_cycler(1)

    return hp.plot_data(data, new_plot_params)


def parametric_plot(
    species: set[str] | list[str],
    data: Any,
    plot_params: dict[str, Any],
) -> Any:
    """Generate a plot overlaying curves for each parameter sweep."""
    max_labels = 15
    current_labels = 0

    try:
        new_plot_params = deepcopy(plot_params)
    except (TypeError, AttributeError):
        new_plot_params = plot_params
    set_plot_units(new_plot_params)

    new_plot_params["frameon"] = False
    new_plot_params["figures"] = []
    new_plot_params["pad"] = 1.5

    previous_parameter = data.ts_model_parameters[0]

    # Update plot to add new curve
    def update_plot(  # noqa: PLR0913  # complex function signature
        spe: str,
        temp_ts: list[int],
        i: int,
        p: Any,
        previous_parameter: Any,
        current_labels: int,
    ) -> tuple[dict[str, Any], list[int], Any]:
        """Append a new curve entry for one species/parameter pair."""
        label = (
            str(spe) + " " + str(previous_parameter)
            if current_labels < max_labels
            else None
        )
        temp_dict = {
            "species_to_plot": spe,
            "time_series": temp_ts,
            spe: {"label": label},
        }
        return temp_dict, [i], p

    # Extracting time-series per parameter in sweep
    for spe in species:
        temp_ts: list[int] = []
        plots: list[dict[str, Any]] = []
        for i, p in enumerate(data.ts_model_parameters):
            if str(p) == str(previous_parameter):
                temp_ts.append(i)
            else:
                current_labels += 1
                temp_dict, temp_ts, previous_parameter = update_plot(
                    spe, temp_ts, i, p, previous_parameter, current_labels
                )
                plots.append(temp_dict)

            # For the final value
            if i == len(data) - 1:
                if str(p) != str(previous_parameter):
                    current_labels += 1
                    temp_dict, temp_ts, previous_parameter = update_plot(
                        spe,
                        [len(data) - 1],
                        len(data) - 1,
                        p,
                        previous_parameter,
                        current_labels,
                    )
                else:
                    current_labels += 1
                    temp_dict, temp_ts, previous_parameter = update_plot(
                        spe,
                        temp_ts,
                        len(data) - 1,
                        p,
                        previous_parameter,
                        current_labels,
                    )

                current_labels = 0
                plots.append(temp_dict)

        new_plot_params["figures"].append({"plots": plots})

    # Adjusting for label size
    if len(data.ts_model_parameters) < _SMALL_LABEL_THRESHOLD:
        prop = {"size": 10}
    elif len(data.ts_model_parameters) < _MEDIUM_LABEL_THRESHOLD:
        prop = {"size": 8}
    else:
        prop = {"size": 6}

    new_plot_params["prop"] = prop
    return hp.plot_data(data, new_plot_params)


def raw_plot(
    data: Any,
    parameters_or_file: dict[str, Any] | str,
    return_fig: bool = False,
) -> Any:
    """
    Plots data from a json or parameter dictionary
    configured according to the hierarchical plot
    structure. Does not accept parameters from a
    Simulation object, it must be given the parameters
    in its entirety.

    Args:
        data: Data in MobsPy format Simulation.results['data'].
        parameters_or_file: Dictionary originated from a JSON or JSON file name.
        return_fig: Return figure instead of plotting.
    """
    if isinstance(parameters_or_file, str) and parameters_or_file[-5:] == ".json":
        plot_params = read_plot_json(parameters_or_file)
    elif isinstance(parameters_or_file, dict):
        plot_params = parameters_or_file
    else:
        raise ValidationError(
            "Raw plot only takes json files or parameters for configuration"
        )

    species = list(data.ts_data[0].keys())
    ppd.check_plot_parameters(species, plot_params)

    # Here we add some base data just in case the user has not supplied it
    return hp.plot_data(data, plot_params, return_fig_object=return_fig)


if __name__ == "__main__":
    raw_plot({}, "default_parameters.json")
