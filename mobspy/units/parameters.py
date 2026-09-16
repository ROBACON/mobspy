"""Resolve symbolic parameter values without changing their declarations."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pint import Quantity
from scipy.constants import N_A

from mobspy.exceptions import ParameterError
from mobspy.units.model_context import ModelUnitContext
from mobspy.units.validation import real_magnitude

if TYPE_CHECKING:
    from mobspy.dsl.mobspy_parameters import Internal_Parameter_Constructor


def parameter_values(
    parameter: Internal_Parameter_Constructor,
    context: ModelUnitContext | None,
    *,
    counts: bool = False,
) -> list[float]:
    """Resolve a sweep into model units, or initial amounts for count bindings."""
    values = parameter.value if isinstance(parameter.value, list) else [parameter.value]
    if not values:
        raise ParameterError(f"Parameter {parameter.name!r} has an empty sweep")
    if not parameter.has_units():
        return [real_magnitude(value) for value in values]
    ctx = context or ModelUnitContext()
    quantities = [
        Quantity(value / parameter.conversion_factor, parameter.original_unit)
        for value in values
    ]
    if counts:
        return [
            real_magnitude(ctx.convert_counts(value, ctx.resolved_volume_magnitude))
            for value in quantities
        ]
    return [_convert_parameter(value, ctx) for value in quantities]


def _convert_parameter(value: Quantity, context: ModelUnitContext) -> float:
    dimensions = dict(value.dimensionality)
    if set(dimensions) - {"[time]", "[length]", "[substance]"}:
        raise ParameterError(f"Unsupported parameter dimensions: {value}")
    length = real_magnitude(dimensions.get("[length]", 0))
    substance = real_magnitude(dimensions.get("[substance]", 0))
    time = real_magnitude(dimensions.get("[time]", 0))
    if length and not context.dimension:
        raise ParameterError("A spatial parameter requires a spatial model")
    target = (
        context.time_unit**time
        * context.volume_unit ** (length / (context.dimension or 1))
        * (context.substance_unit or context.time_unit._REGISTRY.mole) ** substance
    )
    magnitude = real_magnitude(value.to(target).magnitude)
    return magnitude * (N_A**substance if context.substance_unit is None else 1)
