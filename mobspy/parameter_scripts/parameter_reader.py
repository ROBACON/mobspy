"""
This module is responsible for processing the
parameters given to MobsPy before starting the
simulation
"""

from __future__ import annotations

import json
from datetime import datetime
from typing import Any

from pint import Quantity

import mobspy.modules.unit_handler as uh
from mobspy.exceptions import ParameterError
from mobspy.modules.mobspy_expressions import u


def read_json(json_file_name: str) -> Any:
    """
    Reads json file

    :param plot_json_filename: json file name

    :raise simlog.error: If the file was not able to be read

    :return: plot parameter dictionary
    """
    with open(json_file_name) as file:
        try:
            json_data = json.load(file)
        except json.decoder.JSONDecodeError as e:
            raise ParameterError("Error reading file") from e

    return json_data


def name_output_file(params: dict[str, Any]) -> None:
    """
    Gives a name to the output file - just date
    time in case the user has not specified one

    :param params: (dict) Dictionary with simulation parameters
    """

    file_name = "r_"
    file_name += str(datetime.now()) + ".json"
    params["absolute_output_file"] = params["output_dir"] + file_name


def check_stochastic_repetitions_seeds(params: dict[str, Any]) -> None:
    """
    The list of seeds must be equal to the number of repetitions specified

    :param params (dict) = Dictionary with simulation parameters
    """
    if "seeds" in params:
        try:
            if params["repetitions"] != len(params["seeds"]):
                raise ParameterError("Seeds must be equal to the number of repetitions")
        except Exception as e:
            raise ParameterError("Parameter seeds must be a list") from e


def convert_parameters_for_COPASI(params: dict[str, Any]) -> None:
    """
    Converts parameters units to MobsPy standard units
    (basiCO needs seconds for simulation duration)

    :param params: (dict) = Dictionary with simulation parameters
    """
    for key, p in params.items():
        if (
            key == "duration"
            and isinstance(p, Quantity)
            and p.dimensionality != "[time]"
        ):
            raise ParameterError(
                "The duration of the simulation is not in units of time"
            )

        if (
            isinstance(p, Quantity)
            and (key != "unit_x" and key != "unit_y")
            and str(p.dimensionality) == "[time]"
        ):
            params[key] = p.convert("second").magnitude
            continue


def convert_unit_parameters(params: dict[str, Any]) -> None:
    units = ["unit_x", "unit_y"]

    for un in units:
        if un in params and params[un] is not None:
            if isinstance(params[un], Quantity):
                params[un] = params[un]
            else:
                try:
                    params[un] = u.unit_registry_object(params[un])
                except Exception as e:
                    raise ParameterError(
                        f"The unit in parameter {un} did not parse"
                    ) from e


def convert_time_parameters_after_compilation(
    value: int | float | Quantity,
    model_context: Any = None,
) -> int | float | Quantity:
    """
    This function converts the duration if the model was already compiled
    """
    if isinstance(value, Quantity) and str(value.dimensionality) == "[time]":
        if model_context is not None:
            value = model_context.convert_time(value)
        else:
            value = value.convert("second").magnitude
    return value


def convert_volume_after_compilation(
    dimension: int | None,
    parameters_for_sbml: dict[str, Any],
    value: int | float | Quantity,
    model_context: Any = None,
) -> int | float:
    if isinstance(value, Quantity):
        message = (
            f"Error converting volume parameter (value={value}, dimension={dimension})."
            "\n The dimension is set to three at the moment of the compilation if"
            " not specified beforehand.\n Please set a volume in the correct"
            " dimension before compilation."
        )
        uh.extract_length_dimension(
            str(value.dimensionality), dimension, context=message
        )

    value = uh.convert_volume(value, dimension, model_context=model_context)
    if model_context is not None:
        vol_unit_id = model_context.get_sbml_volume_units_id()
    else:
        vol_unit_id = "dimensionless"
    parameters_for_sbml["volume"] = (value, vol_unit_id)
    return value


"""
all basico methods listed here for future reference - if more need to be added
methods = {
        'deterministic': COPASI.CTaskEnum.Method_deterministic,
        'lsoda': COPASI.CTaskEnum.Method_deterministic,
        'hybrid': COPASI.CTaskEnum.Method_hybrid,
        'hybridode45': COPASI.CTaskEnum.Method_hybridODE45,
        'hybridlsoda': COPASI.CTaskEnum.Method_hybridLSODA,
        'adaptivesa': COPASI.CTaskEnum.Method_adaptiveSA,
        'tauleap': COPASI.CTaskEnum.Method_tauLeap,
        'stochastic': COPASI.CTaskEnum.Method_stochastic,
        'directmethod': COPASI.CTaskEnum.Method_directMethod,
        'radau5': COPASI.CTaskEnum.Method_RADAU5,
        'sde': COPASI.CTaskEnum.Method_stochasticRunkeKuttaRI5,
    }
"""


def check_method_parameter(params: dict[str, Any]) -> None:
    # Method takes preference from the user assignment,
    # but the code was made for 'simulation_method
    if params["method"] is not None:
        params["simulation_method"] = params["method"]
    params["simulation_method"] = params["simulation_method"].lower()

    valid_basiCO_deterministic = ["deterministic", "lsoda"]
    valid_basiCO_stochastic = [
        "hybrid",
        "hybridode45",
        "hybridlsoda",
        "tauleap",
        "stochastic",
        "directmethod",
        "sde",
    ]

    if (
        params["simulation_method"]
        not in valid_basiCO_deterministic + valid_basiCO_stochastic
    ):
        raise ParameterError(
            "The simulation method "
            f"{params['simulation_method']}"
            " is not compatible with MobsPy"
        )

    if params["simulation_method"] in valid_basiCO_deterministic:
        if params["rate_type"] is None:
            params["rate_type"] = "deterministic"

        if params["plot_type"] is None:
            params["plot_type"] = "deterministic"
    else:
        if params["rate_type"] is None:
            params["rate_type"] = "stochastic"

        if params["plot_type"] is None and params["repetitions"] == 1:
            params["plot_type"] = "deterministic"
        else:
            params["plot_type"] = "stochastic"


def check_duration_unit(params: dict[str, Any]) -> None:
    if isinstance(params["duration"], Quantity) and params["unit_x"] is None:
        params["unit_x"] = 1 * params["duration"].units


def parameter_process(params: dict[str, Any]) -> None:
    check_duration_unit(params)
    convert_unit_parameters(params)
    name_output_file(params)
    check_stochastic_repetitions_seeds(params)
    convert_parameters_for_COPASI(params)
    check_method_parameter(params)


def manually_process_each_parameter(
    simulation_object: Any,
    duration: int | float | Quantity | None,
    volume: int | float | Quantity | None,
    dimension: int | None,
    repetitions: int | None,
    level: int | None,
    simulation_method: str | None,
    start_time: int | float | None,
    r_tol: float | None,
    a_tol: float | None,
    seeds: list[int] | None,
    step_size: int | float | None,
    jobs: int | None,
    unit_x: str | Quantity | None,
    unit_y: str | Quantity | None,
    output_concentration: bool | None,
    output_event: bool | None,
    output_file: str | None,
    save_data: bool | None,
    plot_data: bool | None,
    rate_type: str | None,
    plot_type: str | None,
) -> None:
    if duration is not None:
        simulation_object.duration = duration

    if volume is not None:
        simulation_object.volume = volume

    if dimension is not None:
        simulation_object.dimension = dimension

    if repetitions is not None:
        simulation_object.repetitions = repetitions

    if level is not None:
        simulation_object.level = level

    if simulation_method is not None:
        simulation_object.simulation_method = simulation_method

    if start_time is not None:
        simulation_object.start_time = start_time

    if r_tol is not None:
        simulation_object.r_tol = r_tol

    if a_tol is not None:
        simulation_object.a_tol = a_tol

    if seeds is not None:
        simulation_object.seeds = seeds

    if step_size is not None:
        simulation_object.step_size = step_size

    if jobs is not None:
        simulation_object.jobs = jobs

    if unit_x is not None:
        simulation_object.unit_x = unit_x

    if unit_y is not None:
        simulation_object.unit_y = unit_y

    if output_concentration is not None:
        simulation_object.output_concentration = output_concentration

    if output_event is not None:
        simulation_object.output_event = output_event

    if output_file is not None:
        simulation_object.output_file = output_file

    if save_data is not None:
        simulation_object.save_data = save_data

    if plot_data is not None:
        simulation_object.plot_data = plot_data

    if rate_type is not None:
        simulation_object.rate_type = rate_type

    if plot_type is not None:
        simulation_object.plot_type = plot_type
