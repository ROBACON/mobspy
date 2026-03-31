"""Convert Pint quantities for rates and counts into dimensionless values."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

from pint import DimensionalityError, Quantity
from scipy.constants import N_A

from mobspy.exceptions import UnitError
from mobspy.modules.mobspy_expressions import OverrideQuantity, u

if TYPE_CHECKING:
    from mobspy.modules.model_unit_context import ModelUnitContext


def convert_rate(
    quantity: int | float | Quantity,  # type: ignore[type-arg]
    reaction_order: int,
    dimension: int | None,
    model_context: ModelUnitContext | None = None,
) -> tuple[float | int | Any, int | None, bool]:
    """
    This function converts the rate from the users given unit to MobsPy standard units

    Args:
        quantity: If it is a quantity object convert, otherwise it remains the same.
        reaction_order: Number of reactants in the reaction, to check if the rate is in
            the correct unit.
        dimension: Model's dimension (1D, 2D, 3D, ... ).
    """

    if model_context is not None:
        return model_context.convert_rate(quantity, reaction_order, dimension)

    warnings.warn(
        "Legacy unit conversion without model_context is deprecated",
        DeprecationWarning,
        stacklevel=2,
    )
    volume_power = reaction_order - 1
    # For objects that cannot be deep copied
    converted_quantity = deep_copy_quantities(quantity)

    # Check to see if rate dimension is valid
    if dimension is None and isinstance(quantity, Quantity) and reaction_order > 1:
        dimension = extract_length_dimension(
            str(quantity.dimensionality), dimension, reaction_order
        )

    if isinstance(quantity, Quantity):
        dim = dict(quantity.dimensionality)
        has_time = "[time]" in dim
        has_substance = "[substance]" in dim
        has_length = "[length]" in dim

        try:
            if has_time and not has_substance and not has_length:
                # Pure 1/[time] rate (count-based, order 0 or 1)
                converted_quantity = converted_quantity.convert("1/seconds")
                return converted_quantity.magnitude, dimension, True
            if has_substance and not has_length:
                # [substance]/[time] rate (e.g. mol/s)
                converted_quantity = converted_quantity.convert("moles/seconds")
                return converted_quantity.magnitude * N_A, dimension, True
            if has_substance:
                # [length]^n/([substance]^m*[time]) concentration rate with moles
                if dimension is None:
                    raise UnitError(
                        "dimension is required for concentration "
                        "rate conversion with substance units"
                    )
                converted_quantity = converted_quantity.convert(
                    f"decimeters ** {dimension * volume_power}"
                    f"/(moles ** {volume_power} * seconds)"
                )
                return (
                    converted_quantity.magnitude / (N_A**volume_power),
                    dimension,
                    False,
                )
            # [length]^n/[time] concentration rate
            if dimension is None:
                raise UnitError(
                    "dimension is required for concentration rate conversion"
                )
            converted_quantity = converted_quantity.convert(
                f"decimeters ** {dimension * volume_power}/seconds"
            )
            return converted_quantity.magnitude, dimension, False
        except (DimensionalityError, TypeError, ValueError) as e:
            raise UnitError(
                str(e) + "\n" + f"Problem converting rate {quantity} \n"
                f"Is the rate in the form [volume]**{volume_power}/[time]?"
            ) from e
    else:
        return quantity, dimension, False


def convert_counts(
    quantity: int | float | Quantity | Any,  # type: ignore[type-arg]
    volume: int | float,
    dimension: int | None,
    model_context: ModelUnitContext | None = None,
) -> Any:
    """
    This function converts the counts from the users
    given unit to MobsPy standard units. It also
    converts concentrations into counts.

    Args:
        quantity: If it is a quantity object convert, otherwise it remains the same.
        volume: Volume in liters (converted beforehand).
        dimension: Model's dimension (1D, 2D, 3D, ... ).
        model_context: Optional ModelUnitContext for user-unit conversion.


    Returns:
        Converted unit into MobsPy standard units.
    """
    if model_context is not None:
        return model_context.convert_counts(quantity, volume)

    warnings.warn(
        "Legacy unit conversion without model_context is deprecated",
        DeprecationWarning,
        stacklevel=2,
    )
    converted_quantity = deep_copy_quantities(quantity)

    if isinstance(quantity, Quantity):
        dim = dict(quantity.dimensionality)
        has_length = "[length]" in dim
        has_substance = "[substance]" in dim
        is_dimensionless = len(dim) == 0

        if not has_length and not has_substance and not is_dimensionless:
            raise UnitError(
                f"The assigned quantity {quantity} is"
                " neither a count nor a concentration"
            )
        if is_dimensionless:
            return quantity.magnitude

        try:
            if has_substance:
                if has_length:
                    converted_quantity = converted_quantity.convert(
                        f"moles/(decimeter ** {dimension})"
                    )
                    converted_quantity = converted_quantity * volume

                converted_quantity = converted_quantity.magnitude * N_A
            else:
                if has_length:
                    converted_quantity = converted_quantity.convert(
                        f"1/(decimeter ** {dimension})"
                    )
                    converted_quantity = converted_quantity * volume
                converted_quantity = converted_quantity.magnitude
        except (DimensionalityError, TypeError, ValueError) as e:
            raise UnitError(
                str(e) + "\n" + f"Problem converting count {quantity} \n"
                "Is it really a count or concentration?"
            ) from e
    return converted_quantity


def check_dimension(
    dimension: int | None,
    value: int | float | str,
    error_context: bool | str = False,
) -> int:
    """
    Checks for dimension consistency. It "stores" the
    first dimension it was given by returning it.

    Args:
        dimension: Model's dimension (1D, 2D, 3D ...).
        value: Dimension value being analysed.
        error_context: Context of the error if dimensions are not consistent.


    Returns:
        Model's dimension (1D, 2D, 3D ...).
    """
    if dimension is None:
        dimension = int(value)
    elif dimension != int(value):
        message = (
            "The dimensions are not consistent. "
            "There are at least two units given "
            "for different dimension models."
        )
        if error_context:
            message = message + "\n " + str(error_context)
        raise UnitError(message)
    return dimension


def extract_length_dimension(
    unit_string: str,
    dimension: int | None,
    reaction_order: int | None = None,
    context: bool | str = False,
) -> int | bool:
    """Extract the volume dimension from a Pint dimensionality.

    Uses Pint's dimensionality dict API instead of string parsing
    for robustness across Pint versions.

    Args:
        unit_string: Unit dimensionality in str format, or a Pint UnitsContainer-like
            object coerced to str.
        dimension: Model's dimension (1D, 2D, 3D ...).
        reaction_order: Number of reactants in a reaction (for dimensional consistency
            in rates).
        context: Context of the error if dimensions are not consistent.
    """
    # Parse the length exponent from the dimensionality string
    # using Pint's UnitRegistry to get the dict representation
    length_power = _extract_length_power(unit_string)

    if length_power is None:
        return False

    if reaction_order is None:
        dimension = check_dimension(dimension, length_power, context)
    elif reaction_order == 1:
        if length_power != 0:
            raise UnitError(
                "Unimolecular reaction (order 1) should not have "
                f"[length] in rate units, but got exponent {length_power}"
            )
        dimension = check_dimension(dimension, 0, context)
    else:
        volume_dim = int(length_power / (reaction_order - 1))
        dimension = check_dimension(dimension, volume_dim, context)

    return dimension


def _extract_length_power(unit_string: str) -> int | None:
    """Extract the [length] exponent from a Pint dimensionality string.

    Parses the dimensionality dict via the unit registry for
    robustness, with a string fallback for edge cases.
    """
    try:
        # Use Pint's API to parse dimensionality
        ur = u.unit_registry_object
        dim_container = ur.parse_expression(unit_string).dimensionality
        length_exp = dim_container.get("[length]", 0)
        if length_exp == 0:
            return None
        return int(length_exp)
    except (DimensionalityError, TypeError, ValueError, AttributeError):
        # Fallback: parse from string representation
        if "[length]" not in unit_string:
            return None
        parts = unit_string.split()
        try:
            pos = parts.index("[length]")
            if pos + 1 < len(parts) and parts[pos + 1] == "**":
                return int(parts[pos + 2])
            return 1
        except (ValueError, IndexError):
            return None


def convert_volume(
    volume: int | float | Quantity,  # type: ignore[type-arg]
    dimension: int | None,
    model_context: ModelUnitContext | None = None,
) -> int | float:
    """
    Converts volume to model volume units (or decimetre**dimension in legacy mode).

    Args:
        volume: Volume used in simulation.
        model_context: Optional ModelUnitContext for user-unit conversion.


    Returns:
        The converted volume in model units.
    """
    if model_context is not None:
        return model_context.convert_volume(volume)

    warnings.warn(
        "Legacy unit conversion without model_context is deprecated",
        DeprecationWarning,
        stacklevel=2,
    )
    if isinstance(volume, Quantity):
        dimension = extract_length_dimension(str(volume.dimensionality), dimension)
        converted: int | float = volume.convert(f"decimeter ** {dimension}").magnitude  # type: ignore[assignment]
        return converted
    return volume  # pyright: ignore[reportReturnType]


def convert_time(
    time: int | float | Quantity,  # type: ignore[type-arg]
    model_context: ModelUnitContext | None = None,
) -> int | float | None:
    """
    Converts time to model time units (or seconds in legacy mode).

    Args:
        time: Any time used.
        model_context: Optional ModelUnitContext for user-unit conversion.
    """
    if model_context is not None:
        return model_context.convert_time(time)

    warnings.warn(
        "Legacy unit conversion without model_context is deprecated",
        DeprecationWarning,
        stacklevel=2,
    )
    if isinstance(time, Quantity):
        dim = dict(time.dimensionality)
        if dim.get("[time]") and len(dim) == 1:
            return time.convert("second").magnitude  # type: ignore[no-any-return]
    else:
        return time  # pyright: ignore[reportReturnType]
    return None


def time_convert_to_other_unit(
    time: int | float | Quantity,  # type: ignore[type-arg]
    other_unit: str,
) -> int | float | None:
    """Converts time to a specified unit.

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
    """
    Custom deepcopy function for Pint Quantity objects.

    Args:
        quantity: The Quantity object to copy.

    Returns:
        A new Quantity object with the same value and units.
    """
    # Extract the value and unit from the original quantity
    if isinstance(quantity, Quantity):
        value = quantity.magnitude
        unit = str(quantity.units)  # Convert to string to ensure unit compatibility

        # Create a new Quantity using the OverrideUnitRegistry object (u)
        Q = value * u.unit_registry_object.__getattr__(unit)  # noqa: N806
        return OverrideQuantity(Q)
    return quantity
