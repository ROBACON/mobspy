"""
mobspy.data_handler.process_result_data.py

Handles converting the output data from a simulation into desired-units or concentration
"""

from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING, Any

from scipy.constants import N_A

from mobspy.mobspy_logging import get_logger
from mobspy.modules.unit_registry import u

if TYPE_CHECKING:
    from mobspy.types import SimulationParameters

_logger = get_logger(__name__)


def extract_time_and_volume_list(
    list_of_params: list[SimulationParameters],
) -> tuple[list[float], list[float], bool]:
    """
    This function extracts the list of durations from all
    concatenated simulations (or one for single simulation).
    It also extracts the respective volumes at each
    simulation. It does not work if there is a change in
    volume and a simulation without fixed duration.

    Conditional simulations have a parameter called '_end_condition' in their dictionary

    Args:
        list_of_params: List of parameters of all concatenated simulations.
    """
    no_fixed_volume = False
    no_fixed_dur = False
    initial_volume: float = list_of_params[0]["volume"]
    for par in list_of_params:
        if par["_end_condition"] is not None:
            no_fixed_dur = True

        if par["volume"] != initial_volume:
            no_fixed_volume = True

    flag_concentration = True
    if no_fixed_dur and len(list_of_params) > 1 and no_fixed_volume:
        flag_concentration = False
        _logger.warning(
            "Could not resolve simulation volume due to "
            "multiple simulations with at least one "
            "with a conditional duration. The output "
            "will be printed in counts instead"
        )

    volume_list: list[float] = []
    previous_time: float = list_of_params[0]["duration"]
    sim_time_list: list[float] = [previous_time]
    if no_fixed_volume:
        for i, par in enumerate(list_of_params):
            volume_list.append(par["volume"])

            if i == 1:
                continue
            current_time: float = previous_time + par["duration"]
            sim_time_list.append(current_time)
            previous_time = current_time
    else:
        volume_list = [initial_volume]
        sim_time_list = [sum([par["duration"] for par in list_of_params])]

    return volume_list, sim_time_list, flag_concentration


def convert_data_to_desired_unit(  # noqa: PLR0913
    data: dict[str, list[float]],
    time_list: list[float],
    volume_list: list[float],
    unit_x: str | None = None,
    unit_y: str | None = None,
    output_concentration: bool = False,
    model_context: Any = None,
) -> dict[str, list[float]]:
    """Converts the simulation output data from the MobsPy standard units
    to the desired units specified by the user

    Args:
        data: Resulting data from a MobsPy simulation execution.
        time_list: List of times where the volume changes (single simulation: only one).
        volume_list: List of volumes changes (single simulation: only one).
        unit_x: Unit that the user desires the time in, defaults to dimensionless
            (None).
        unit_y: Unit that the user desires the y axis to be in (Concentration or
            counts).
        output_concentration: Decide if output should be a concentration or count.


    Returns:
        Input data converted to the desired units.
    """
    ur = u.unit_registry_object
    converted_data = deepcopy(data)

    if unit_x is not None:
        time_unit = model_context.time_unit if model_context is not None else ur.seconds
        converted_data["Time"] = [
            (time * time_unit).to(unit_x).magnitude  # pyright: ignore[reportAttributeAccessIssue]
            for time in data["Time"]
        ]

    if output_concentration:
        converted_data = convert_to_concentration(
            data, converted_data, volume_list, time_list
        )

    if unit_y is not None:
        _apply_unit_y_conversion(
            converted_data,
            unit_y,
            output_concentration,
            model_context,
            ur,
        )

    return converted_data


def _multiply_data_by_factor(
    converted_data: dict[str, list[float]],
    factor: float,
) -> None:
    """Scale all non-Time species data by a constant factor."""
    for key in list(converted_data):
        if key == "Time":
            continue
        converted_data[key] = [count * factor for count in converted_data[key]]


def _apply_unit_y_conversion(
    converted_data: dict[str, list[float]],
    unit_y: str,
    output_concentration: bool,
    model_context: Any,
    ur: Any,
) -> None:
    """Apply substance/concentration unit conversion on the y-axis data."""
    _substance_is_molar = model_context is not None and model_context.substance_is_molar

    if "mol" in str(unit_y):
        if not _substance_is_molar:
            _multiply_data_by_factor(converted_data, N_A**-1)
        if output_concentration:
            if _substance_is_molar:
                source = 1 * model_context.substance_unit / model_context.volume_unit
            else:
                source = 1 * ur.molar
            _multiply_data_by_factor(
                converted_data,
                source.to(unit_y).magnitude,  # pyright: ignore[reportAttributeAccessIssue]
            )
        else:
            if _substance_is_molar:
                source = 1 * model_context.substance_unit
            else:
                source = 1 * ur.moles
            _multiply_data_by_factor(
                converted_data,
                source.to(unit_y).magnitude,  # pyright: ignore[reportAttributeAccessIssue]
            )
    elif output_concentration:
        source = 1 / model_context.volume_unit if _substance_is_molar else 1 / ur.l
        _multiply_data_by_factor(
            converted_data,
            source.to(unit_y).magnitude,  # pyright: ignore[reportAttributeAccessIssue]
        )


def convert_to_concentration(
    data: dict[str, list[float]],
    converted_data: dict[str, list[float]],
    volume_list: list[float],
    time_list: list[float],
) -> dict[str, Any]:
    """
    Converts output data from counts to concentration according to simulation volume

    Args:
        data: Simulation data.
        converted_data: Data converted to requested units.
        volume_list: List of volumes of all simulations (more than one if concatenated).
        time_list: List of durations of each simulation (to check for respective volume
            in results).
    """
    new_data: dict[str, Any] = {}
    for key in data:
        new_data[key] = []

    current_volume = volume_list[0]
    k = 0
    for i in range(len(data["Time"])):
        if data["Time"][i] > time_list[k] and len(volume_list) > 1:
            k = k + 1
            current_volume = volume_list[k]

        for key in data:
            if key == "Time":
                continue
            new_data[key].append(data[key][i] / current_volume)

    new_data["Time"] = converted_data["Time"]

    return new_data
