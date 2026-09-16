"""Unit conversion between independently compiled simulation stages."""

from __future__ import annotations

from scipy.constants import N_A

from mobspy.constants import END_FLAG_SPECIES_NAME
from mobspy.types import RawTimeSeries
from mobspy.units.model_context import ModelUnitContext


def amount_factor(source: ModelUnitContext, target: ModelUnitContext) -> float:
    """Number of target substance units per source substance unit."""
    source_items = (
        float((1 * source.substance_unit).to("mole").magnitude) * N_A
        if source.substance_unit is not None
        else 1.0
    )
    target_items = (
        float((1 * target.substance_unit).to("mole").magnitude) * N_A
        if target.substance_unit is not None
        else 1.0
    )
    return source_items / target_items


def normalize_time_series(
    data: RawTimeSeries,
    source: ModelUnitContext | None,
    target: ModelUnitContext | None,
) -> RawTimeSeries:
    """Express a stage's native output in the first stage's unit system."""
    if source is None or target is None:
        return {name: list(values) for name, values in data.items()}
    time_factor = float((1 * source.time_unit).to(target.time_unit).magnitude)
    concentration_factor = amount_factor(source, target) * float(
        (1 * target.volume_unit).to(source.volume_unit).magnitude
    )
    return {
        name: [
            value * (time_factor if name == "Time" else concentration_factor)
            for value in values
        ]
        if name != END_FLAG_SPECIES_NAME
        else list(values)
        for name, values in data.items()
    }
