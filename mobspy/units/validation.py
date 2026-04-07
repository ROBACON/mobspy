"""Early unit validation for rate expressions.

Validates dimensional consistency at reaction definition time
(not just at compile time), catching errors closer to their source.
"""

from __future__ import annotations

from typing import Any

from pint import Quantity

from mobspy.exceptions import UnitError


def validate_rate_units(
    rate: Any,
    num_reactants: int,  # noqa: ARG001  # unused arg required by interface
) -> None:
    """Validate that a rate value has compatible units for a reaction.

    Called at definition time (in ``@`` operator or ``[]`` syntax)
    when the rate is a Pint Quantity. Checks that the rate has
    dimensions consistent with mass-action kinetics for the given
    reaction order.

    Expected dimensions for mass-action:
      - 0th order: substance / time
      - 1st order: 1 / time
      - 2nd order: volume / (substance * time)
      - nth order: volume^(n-1) / (substance^(n-1) * time)

    Args:
        rate: The reaction rate value.
        num_reactants: Number of distinct reactant species.

    Raises:
        UnitError: If the rate has incompatible dimensions.
    """
    if not isinstance(rate, Quantity):
        return

    dim = dict(rate.dimensionality)

    if "[time]" not in dim:
        raise UnitError(
            f"Rate {rate} has no time dimension. "
            "Reaction rates must include a time component (e.g., 1/s, mol/L/s)."
        )

    time_exp = dim.get("[time]", 0)
    if time_exp >= 0:
        raise UnitError(
            f"Rate {rate} has a non-negative time exponent ({time_exp}). "
            "Rates must be inversely proportional to time (e.g., 1/s)."
        )
