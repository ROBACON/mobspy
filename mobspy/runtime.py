"""Build execution plans without modifying source models or configuration."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace
from itertools import product
from typing import Any

from mobspy.constants import SBML_LOCATION
from mobspy.exceptions import ParameterError
from mobspy.types import ConcreteModel, ExecutionPlan, RunSettings, SpeciesDict
from mobspy.units.parameters import parameter_values


def build_execution_plan(
    models: Sequence[ConcreteModel], parameters: Sequence[dict[str, Any]]
) -> ExecutionPlan:
    """Expand shared parameters into independent compiled model variants."""
    values: dict[str, list[int | float]] = {}
    for model in models:
        for name, info in model.parameters_used.items():
            raw = info.object.value if info.object is not None else info.values
            options = raw if isinstance(raw, list) else [raw]
            if not options:
                raise ParameterError(f"Parameter {name!r} has an empty sweep")
            if name in values and values[name] != options:
                raise ParameterError(
                    f"Shared parameter {name!r} has different values across stages"
                )
            values[name] = list(options)
    names = sorted(values)
    selections = tuple(product(*(range(len(values[name])) for name in names)))
    combinations = tuple(
        {name: values[name][index] for name, index in zip(names, selected, strict=True)}
        for selected in selections
    )
    chains = tuple(
        tuple(
            _with_parameters(model, dict(zip(names, selected, strict=True)))
            for model in models
        )
        for selected in selections
    )
    return ExecutionPlan(
        models=chains,
        settings=tuple(RunSettings.from_parameters(p) for p in parameters),
        parameter_values=combinations,
    )


def _with_parameters(model: ConcreteModel, selected: dict[str, int]) -> ConcreteModel:
    species = SpeciesDict(model.species)
    parameters = dict(model.parameters)
    for name, info in model.parameters_used.items():
        values = info.values if isinstance(info.values, list) else [info.values]
        value = values[selected[name]]
        for location in info.used_in:
            if location == SBML_LOCATION:
                unit = parameters.get(name, (0, "dimensionless"))[1]
                parameters[name] = (value, unit)
            else:
                species[location] = (
                    parameter_values(info.object, model.unit_context, counts=True)[
                        selected[name]
                    ]
                    if info.object is not None
                    else value
                )
    return replace(model, species=species, parameters=parameters)
