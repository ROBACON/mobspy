from __future__ import annotations

from copy import deepcopy
from typing import Any

from mobspy.exceptions import ValidationError
from mobspy.mobspy_logging import get_logger

simlog = get_logger(__name__)
import mobspy.plot_params.example_plot_reader as epr  # noqa: E402


def query_plot_data(species: set[str] | list[str], data: Any) -> tuple[list[str], Any]:
    """
    Performs a query of the plot data, when one wishes to plot
    species with characteristics. It creates a new data structure
    with the results from the query added to it for the plotting
    structure

    :param species: (str) Species name in string format
    :param data: (dict) Data in MobsPy dictionary format
    """
    new_data = deepcopy(data)

    species_to_plot: set[str] = set()
    for time_series in data.ts_data:
        for key in time_series.keys():  # noqa: SIM118
            if key in species:
                species_to_plot.add(key)

    for spe in species:
        if spe in species_to_plot:
            continue

        for i, _ in enumerate(data.ts_data):
            new_data[i][spe] = data[spe][i]

        species_to_plot.add(spe)

    species_to_plot_list = sorted(species)
    return species_to_plot_list, new_data


def check_plot_parameters(species: list[str], plot_params: dict[str, Any]) -> None:
    """
    Performs a check of the plot_parameters given. To see if the
    parameters are correctly named

    :param species: (str) Species in str format
    :param plot_params: (dict) Plot parameter dictionary
    """
    dictionary = epr.get_example_plot_parameters()

    if "Time" in plot_params:
        raise ValidationError("Time must not be a plot parameter name")

    for spe in species:
        if spe in dictionary.keys():  # noqa: SIM118
            raise ValidationError(
                f"Plotting is impossible, species {spe} is a parameter name"
            )

    # Check if parameters are valid
    validated_keys: set[str] = set()
    for key in plot_params:
        if key not in dictionary and key not in species:
            continue
        validated_keys.add(key)

    # Check if query is present
    for key in plot_params:
        if key in validated_keys:
            continue
        spe_name = key.split(".")[0]
        if spe_name not in species:
            simlog.warning(f"Parameter {key} not supported")
        validated_keys.add(key)


def time_filter_operation(
    low: float, high: float, time_data: list[float], data: list[float]
) -> tuple[list[float], list[float]]:
    new_time_data: list[float] = []
    new_data: list[float] = []

    for t, d in zip(time_data, data):
        if t < low:
            continue
        if low < t < high:
            new_time_data.append(t)
            new_data.append(d)
        elif t > high:
            break

    return new_time_data, new_data


def y_filter_operation(
    low_y: float, high_y: float, time_data: list[float], data: list[float]
) -> tuple[list[float], list[float]]:
    new_time_data: list[float] = []
    new_data: list[float] = []

    for t, d in zip(time_data, data):
        if d < low_y:
            continue
        if low_y <= d <= high_y:
            new_time_data.append(t)
            new_data.append(d)
        elif d > high_y:
            continue

    return new_time_data, new_data
