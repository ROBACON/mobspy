"""Convert backend results without mutating simulation or model state."""

from __future__ import annotations

from typing import Any

from mobspy.exceptions import SimulationError
from mobspy.results.process_data import (
    convert_data_to_desired_unit,
    extract_time_and_volume_list,
)
from mobspy.results.time_series import MobsPyTimeSeries, SimulationResults
from mobspy.types import BackendResults, ExecutionPlan, TimeSeriesDataDict


def process_results(
    raw_results: BackendResults,
    plan: ExecutionPlan,
    parameters: list[dict[str, Any]],
) -> tuple[SimulationResults, SimulationResults]:
    """Convert a backend's native output using the plan's first-stage units."""
    if len(raw_results) != len(plan.models):
        raise SimulationError("Backend result count does not match the execution plan")
    output = dict(parameters[0])
    first_model = plan.models[0][0]
    context = first_model.unit_context
    resolved = [dict(p) for p in parameters]
    if context is not None:
        for params, model in zip(resolved, plan.models[0], strict=True):
            source = model.unit_context
            if source is not None:
                params["volume"] *= (
                    (1 * source.volume_unit).to(context.volume_unit).magnitude
                )
                params["duration"] *= (
                    (1 * source.time_unit).to(context.time_unit).magnitude
                )
    volumes, times, concentration_available = extract_time_and_volume_list(resolved)
    unit_y = output["unit_y"]
    if not concentration_available or (
        unit_y is not None and "[length]" not in unit_y.dimensionality
    ):
        output["output_concentration"] = False
    models = [
        model.to_compiled_model(dict(model.species), dict(model.mappings))
        for model in plan.models[0]
    ]
    parameter_objects = {
        name: obj
        for model in plan.models[0]
        for name, obj in model.parameter_objects.items()
    }
    combinations = plan.parameter_values or tuple({} for _ in plan.models)
    series: list[MobsPyTimeSeries] = []
    for runs, values in zip(raw_results, combinations, strict=True):
        if len(runs) != plan.settings[0].repetitions:
            raise SimulationError(
                "Backend repetition count does not match the execution plan"
            )
        for data in runs:
            converted = convert_data_to_desired_unit(
                data,
                times,
                volumes,
                output["unit_x"],
                output["unit_y"],
                output["output_concentration"],
                model_context=context,
            )
            series.append(
                MobsPyTimeSeries(TimeSeriesDataDict(converted, output, models), values)
            )
    return SimulationResults(series, parameter_objects), SimulationResults(
        series[:1], None, True
    )
