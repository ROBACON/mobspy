"""Expand parametric sweep definitions into concrete model variants."""

from __future__ import annotations

import contextlib
import itertools
import logging
from copy import deepcopy
from typing import TYPE_CHECKING, Any

from mobspy.constants import DOT_SEPARATOR, SBML_LOCATION

_logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from mobspy.types import CompiledModelDict, ParameterSweepList, ParameterUsedInfo


def assign_values_to_model(
    parameter_name: str,
    parameter_value: int | float,
    models: list[CompiledModelDict],
    locations: set[str],
) -> None:
    """Write a parameter value into compiled model dicts at each location."""
    for location in locations:
        if location == SBML_LOCATION:
            for model in models:
                with contextlib.suppress(KeyError):
                    existing = model.parameters_for_sbml.get(parameter_name)
                    unit_str = existing[1] if existing else "dimensionless"
                    model.parameters_for_sbml[parameter_name] = (
                        parameter_value,
                        unit_str,
                    )
        else:
            for model in models:
                try:
                    model.species_for_sbml[location] = parameter_value
                    dot_loc = location.replace(DOT_SEPARATOR, ".")
                    model.species_not_mapped[dot_loc] = parameter_value
                except KeyError:
                    pass


def generate_all_sbml_models(
    model_parameters: dict[str, ParameterUsedInfo],
    list_of_models: list[CompiledModelDict],
) -> tuple[ParameterSweepList, list[dict[str, int | float]]]:
    """Build model copies for every combination of swept parameters."""
    names: list[str] = []
    used_in: list[set[str]] = []
    values: list[list[int | float]] = []

    to_return: ParameterSweepList = []

    if not model_parameters:
        return [list_of_models], []

    keys = sorted(model_parameters.keys())
    for key in keys:
        item = model_parameters[key]

        names.append(item.name)
        used_in.append(item.used_in)
        if isinstance(item.values, list):
            if len(item.values):
                values.append(item.values)
        else:
            values.append([item.values])

    parameter_list_of_dic: list[dict[str, int | float]] = []
    for v in itertools.product(*values):
        parameter_dic: dict[str, int | float] = {}
        models_copy = deepcopy(list_of_models)
        for i, name in enumerate(names):
            assign_values_to_model(name, v[i], models_copy, used_in[i])
            parameter_dic[name] = v[i]

        parameter_list_of_dic.append(parameter_dic)
        to_return.append(models_copy)

    return to_return, parameter_list_of_dic


def unite_parameter_dictionaries(
    dict_1: dict[str, Any],
    dict_2: dict[str, Any],
) -> dict[str, Any]:
    """Merge two parameter-usage dicts, unifying ``used_in`` sets."""
    for key in dict_2:  # noqa: PLC0206  # iterating dict keys directly
        if key not in dict_1:
            dict_1[key] = dict_2[key]
        else:
            if dict_1[key].values != dict_2[key].values:
                _logger.warning(
                    "Parameter '%s' has different values in composed simulations "
                    "(%s vs %s). Using the first value.",
                    key,
                    dict_1[key].values,
                    dict_2[key].values,
                )
            new_used_in = dict_1[key].used_in.union(dict_2[key].used_in)
            dict_1[key].used_in = new_used_in

    return dict_1
