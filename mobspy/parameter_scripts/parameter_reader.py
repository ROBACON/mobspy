"""
This module is responsible for processing the
parameters given to MobsPy before starting the
simulation
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any

from pint import Quantity

import mobspy.modules.unit_handler as uh
from mobspy.exceptions import ParameterError
from mobspy.modules.unit_registry import u


def read_json(json_file_name: str) -> Any:
    """
    Reads json file

    Args:
        plot_json_filename: Json file name.


    Returns:
        Plot parameter dictionary.
    """
    with Path(json_file_name).open(encoding="utf-8") as file:
        try:
            json_data = json.load(file)
        except json.decoder.JSONDecodeError as e:
            raise ParameterError("Error reading file") from e

    return json_data


def name_output_file(params: dict[str, Any]) -> None:
    """
    Gives a name to the output file - just date
    time in case the user has not specified one

    Args:
        params: Dictionary with simulation parameters.
    """

    file_name = "r_"
    file_name += str(datetime.now()) + ".json"
    params["absolute_output_file"] = params["output_dir"] + file_name


def check_stochastic_repetitions_seeds(params: dict[str, Any]) -> None:
    """
    The list of seeds must be equal to the number of repetitions specified

    """
    if "seeds" in params and params["seeds"] is not None:
        try:
            if params["repetitions"] != len(params["seeds"]):
                raise ParameterError("Seeds must be equal to the number of repetitions")
        except TypeError as e:
            raise ParameterError("Parameter seeds must be a list") from e


def convert_parameters_for_COPASI(params: dict[str, Any]) -> None:  # noqa: N802
    """
    Converts parameters units to MobsPy standard units
    (basiCO needs seconds for simulation duration)

    Args:
        params: Dictionary with simulation parameters.
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
            and (key not in {"unit_x", "unit_y"})
            and str(p.dimensionality) == "[time]"
        ):
            params[key] = p.convert("second").magnitude
            continue


def convert_unit_parameters(params: dict[str, Any]) -> None:
    """Parse unit_x/unit_y strings into Pint Quantity objects."""
    units = ["unit_x", "unit_y"]

    for un in units:
        if un in params and params[un] is not None:
            if isinstance(params[un], Quantity):
                params[un] = params[un]
            else:
                try:
                    params[un] = u.unit_registry_object(params[un])
                except (ValueError, AttributeError) as e:
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
    """Convert and store a volume value in SBML parameters."""
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
    """Resolve simulation method and infer rate/plot types."""
    if params["method"] is not None:
        params["simulation_method"] = params["method"]
    params["simulation_method"] = params["simulation_method"].lower()

    valid_basiCO_deterministic = ["deterministic", "lsoda"]  # noqa: N806
    valid_basiCO_stochastic = [  # noqa: N806
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
    """Auto-set unit_x from duration when it carries a Pint unit."""
    if isinstance(params["duration"], Quantity) and params["unit_x"] is None:
        params["unit_x"] = 1 * params["duration"].units


def parameter_process(params: dict[str, Any]) -> None:
    """Run all parameter validation and conversion steps."""
    check_duration_unit(params)
    convert_unit_parameters(params)
    name_output_file(params)
    check_stochastic_repetitions_seeds(params)
    convert_parameters_for_COPASI(params)
    check_method_parameter(params)


_VALID_RUN_PARAMS: frozenset[str] = frozenset(
    {
        "duration",
        "volume",
        "dimension",
        "repetitions",
        "level",
        "simulation_method",
        "start_time",
        "r_tol",
        "a_tol",
        "seeds",
        "step_size",
        "jobs",
        "unit_x",
        "unit_y",
        "output_concentration",
        "output_event",
        "output_file",
        "save_data",
        "plot_data",
        "rate_type",
        "plot_type",
    }
)


def manually_process_each_parameter(
    simulation_object: Any,
    **kwargs: Any,
) -> None:
    """Apply non-None keyword arguments onto a simulation object."""
    for attr_name, attr_value in kwargs.items():
        if attr_name not in _VALID_RUN_PARAMS:
            msg = f"Unknown run parameter: {attr_name!r}"
            raise TypeError(msg)
        if attr_value is not None:
            setattr(simulation_object, attr_name, attr_value)
