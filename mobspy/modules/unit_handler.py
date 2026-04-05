"""Convert Pint quantities for rates and counts into dimensionless values.

All conversion logic lives in :class:`ModelUnitContext`.  The public
functions in this module are thin wrappers that create a default
context when none is provided, ensuring a single code path for
unit conversion.

Dimension extraction and validation are delegated to
:class:`ModelUnitContext` static methods.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from pint import Quantity

from mobspy.modules.mobspy_expressions import OverrideQuantity
from mobspy.modules.unit_registry import u

if TYPE_CHECKING:
    from mobspy.modules.model_unit_context import ModelUnitContext


def _default_context(dimension: int | None = None) -> ModelUnitContext:
    """Create a default ModelUnitContext for legacy callers."""
    from mobspy.modules.model_unit_context import ModelUnitContext  # noqa: PLC0415

    dim = dimension if dimension is not None else 3
    return ModelUnitContext(dimension=dim)


def convert_rate(
    quantity: int | float | Quantity,  # type: ignore[type-arg]
    reaction_order: int,
    dimension: int | None,
    model_context: ModelUnitContext | None = None,
) -> tuple[float | int | Any, int | None, bool]:
    """Convert a rate from user units to model-standard units.

    Delegates entirely to ``ModelUnitContext.convert_rate``.
    """
    ctx = model_context or _default_context(dimension)
    return ctx.convert_rate(quantity, reaction_order, dimension)


def convert_counts(
    quantity: int | float | Quantity | Any,  # type: ignore[type-arg]
    volume: int | float,
    dimension: int | None,
    model_context: ModelUnitContext | None = None,
) -> Any:
    """Convert counts from user units to model-standard units.

    Delegates entirely to ``ModelUnitContext.convert_counts``.
    """
    ctx = model_context or _default_context(dimension)
    return ctx.convert_counts(quantity, volume)


def check_dimension(
    dimension: int | None,
    value: int | float | str,
    error_context: bool | str = False,
) -> int:
    """Check dimension consistency.

    Delegates to ``ModelUnitContext.check_dimension``.
    """
    from mobspy.modules.model_unit_context import ModelUnitContext  # noqa: PLC0415

    return ModelUnitContext.check_dimension(dimension, value, error_context)


def extract_length_dimension(
    unit_string: str,
    dimension: int | None,
    reaction_order: int | None = None,
    context: bool | str = False,
) -> int | bool:
    """Extract the volume dimension from a Pint dimensionality.

    Delegates to ``ModelUnitContext.extract_length_dimension``.
    """
    from mobspy.modules.model_unit_context import ModelUnitContext  # noqa: PLC0415

    return ModelUnitContext.extract_length_dimension(
        unit_string, dimension, reaction_order, context
    )


def convert_volume(
    volume: int | float | Quantity,  # type: ignore[type-arg]
    dimension: int | None,
    model_context: ModelUnitContext | None = None,
) -> int | float:
    """Convert volume to model volume units.

    Delegates entirely to ``ModelUnitContext.convert_volume``.
    """
    ctx = model_context or _default_context(dimension)
    return ctx.convert_volume(volume)


def convert_time(
    time: int | float | Quantity,  # type: ignore[type-arg]
    model_context: ModelUnitContext | None = None,
) -> int | float | None:
    """Convert time to model time units.

    Delegates entirely to ``ModelUnitContext.convert_time``.
    """
    ctx = model_context or _default_context()
    return ctx.convert_time(time)


def time_convert_to_other_unit(
    time: int | float | Quantity,  # type: ignore[type-arg]
    other_unit: str,
) -> int | float | None:
    """Convert time to a specified unit.

    Args:
        time: Any time used.
        other_unit: Target unit string.
    """
    if isinstance(time, Quantity):
        dim = dict(time.dimensionality)
        if dim.get("[time]") and len(dim) == 1:
            return time.convert(other_unit).magnitude  # type: ignore[no-any-return]
    else:
        return time  # pyright: ignore[reportReturnType]
    return None


def deep_copy_quantities(quantity: Any) -> Any:
    """Custom deepcopy for Pint Quantity objects.

    Args:
        quantity: The Quantity object to copy.

    Returns:
        A new Quantity object with the same value and units.
    """
    if isinstance(quantity, Quantity):
        value = quantity.magnitude
        unit = str(quantity.units)
        Q = value * u.unit_registry_object.__getattr__(unit)  # noqa: N806
        return OverrideQuantity(Q)
    return quantity
