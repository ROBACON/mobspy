from __future__ import annotations

from typing import TYPE_CHECKING, Any

from pint import Quantity
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

    :param quantity: (int, float, Quantity) If it is a
        quantity object convert, otherwise it remains
        the same
    :param reaction_order: (int) number of reactants in
        the reaction, to check if the rate is in the
        correct unit
    :param dimension: (int) model's dimension
        (1D, 2D, 3D, ... )

    :param quantity: (int, float) converted unit into
        MobsPy standard units
    """

    if model_context is not None:
        return model_context.convert_rate(quantity, reaction_order, dimension)

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
            elif has_substance and not has_length:
                # [substance]/[time] rate (e.g. mol/s)
                converted_quantity = converted_quantity.convert("moles/seconds")
                return converted_quantity.magnitude * N_A, dimension, True
            elif has_substance:
                # [length]^n/([substance]^m*[time]) concentration rate with moles
                converted_quantity = converted_quantity.convert(
                    f"decimeters ** {dimension * volume_power}"
                    f"/(moles ** {volume_power} * seconds)"
                )
                return (
                    converted_quantity.magnitude / (N_A**volume_power),
                    dimension,
                    False,
                )
            else:
                # [length]^n/[time] concentration rate
                converted_quantity = converted_quantity.convert(
                    f"decimeters ** {dimension * volume_power}/seconds"
                )
                return converted_quantity.magnitude, dimension, False
        except Exception as e:
            raise UnitError(
                str(e) + "\n" + f"Problem converting rate {quantity} \n"
                f"Is the rate in the form [volume]**{volume_power}/[time]?"
            ) from e
    else:
        return quantity, dimension, False


def convert_counts(
    quantity: int | float | Quantity | Any,  # type: ignore[type-arg]
    volume: int | float,
    dimension: int,
    model_context: ModelUnitContext | None = None,
) -> Any:
    """
    This function converts the counts from the users
    given unit to MobsPy standard units. It also
    converts concentrations into counts.

    :param quantity: (int, float, Quantity) If it is a
        quantity object convert, otherwise it remains
        the same
    :param volume: (int, float) volume in liters (converted beforehand)
    :param dimension: (int) model's dimension (1D, 2D, 3D, ... )
    :param model_context: optional ModelUnitContext for user-unit conversion

    :return: converted_quantity (int, float) = converted unit into MobsPy standard units
    """
    if model_context is not None:
        return model_context.convert_counts(quantity, volume)

    converted_quantity = deep_copy_quantities(quantity)

    if isinstance(quantity, Quantity):
        dim = dict(quantity.dimensionality)
        has_length = "[length]" in dim
        has_substance = "[substance]" in dim
        is_dimensionless = len(dim) == 0

        if not has_length and not has_substance and not is_dimensionless:
            raise UnitError(
                f"The assigned quantity {quantity} is neither a count or concentration"
            )
        elif is_dimensionless:
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
        except Exception as e:
            raise UnitError(
                str(e) + "\n" + f"Problem converting rate {quantity} \n"
                f"Is it really a count or concentration?"
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

    :param dimension: (int) model's dimension
        (1D, 2D, 3D ...)
    :param value: (int) dimension value being analysed
    :param error_context: (bool or str) context of the
        error if dimensions are not consistent

    :raise simlog.error: If dimensions are not consistent
        through the given units
        (units in 1D with 2D mixed)

    :return: dimension (int) = model's dimension (1D, 2D, 3D ...)
    """
    if dimension is None:
        dimension = int(value)
    else:
        if dimension != int(value):
            message = (
                "The dimensions are not consistent. "
                "There are at least two units given "
                "for different dimension models."
            )
            if error_context:
                message = message + "\n " + error_context
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

    :param unit_string: (str) unit dimensionality in str format,
        or a Pint UnitsContainer-like object coerced to str
    :param dimension: (int) model's dimension (1D, 2D, 3D ...)
    :param reaction_order: (int) number of reactants in
        a reaction (for dimensional consistency in rates)
    :param context: (bool or str) context of the error if
        dimensions are not consistent
    """
    # Parse the length exponent from the dimensionality string
    # using Pint's UnitRegistry to get the dict representation
    length_power = _extract_length_power(unit_string)

    if length_power is None:
        return False

    if reaction_order is None:
        dimension = check_dimension(dimension, length_power, context)
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
    except Exception:
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

    :param volume: (int, float, Quantity) volume used in simulation
    :param model_context: optional ModelUnitContext for user-unit conversion

    :return: the converted volume in model units
    """
    if model_context is not None:
        return model_context.convert_volume(volume)

    if isinstance(volume, Quantity):
        dimension = extract_length_dimension(str(volume.dimensionality), dimension)
        volume = volume.convert(f"decimeter ** {dimension}").magnitude
        return volume
    else:
        return volume


def convert_time(
    time: int | float | Quantity,  # type: ignore[type-arg]
    model_context: ModelUnitContext | None = None,
) -> int | float | None:
    """
    Converts time to model time units (or seconds in legacy mode).

    :param time: (int, float, Quantity) any time used
    :param model_context: optional ModelUnitContext for user-unit conversion
    """
    if model_context is not None:
        return model_context.convert_time(time)

    if isinstance(time, Quantity):
        dim = dict(time.dimensionality)
        if dim.get("[time]") and len(dim) == 1:
            return time.convert("second").magnitude
    else:
        return time
    return None


def time_convert_to_other_unit(
    time: int | float | Quantity,  # type: ignore[type-arg]
    other_unit: str,
) -> int | float | None:
    """Converts time to a specified unit.

    :param time: (int, float, Quantity) any time used
    :param other_unit: target unit string
    """
    if isinstance(time, Quantity):
        dim = dict(time.dimensionality)
        if dim.get("[time]") and len(dim) == 1:
            return time.convert(other_unit).magnitude
    else:
        return time
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
        Q = value * u.unit_registry_object.__getattr__(unit)
        return OverrideQuantity(Q)
    else:
        return quantity


if __name__ == "__main__":
    A = 5 * u.meters
    B = deep_copy_quantities(A)
    print(B)
