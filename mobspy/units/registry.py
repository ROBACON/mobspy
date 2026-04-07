"""MobsPy unit registry: wraps pint to return OverrideQuantity instances.

This module breaks the import cycle where ``model_unit_context``,
``unit_handler``, and other modules previously imported ``u``
from the heavyweight ``mobspy_expressions`` module.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from pint import UnitRegistry

if TYPE_CHECKING:
    from mobspy.expressions.evaluation import OverrideQuantity


def _make_override_quantity(q_object: Any) -> OverrideQuantity:
    """Lazy import to avoid circular dependency with expression_definer."""
    from mobspy.expressions.evaluation import OverrideQuantity  # noqa: PLC0415

    return OverrideQuantity(q_object)  # pyright: ignore[reportReturnType]


class OverrideUnitRegistry:
    """Pint UnitRegistry wrapper that returns OverrideQuantity instances.

    This allows quantities created via ``u.hour`` or ``u(1, "liter")``
    to participate in MobsPy expression building when ``_ms_active``
    is set.
    """

    def __init__(self) -> None:
        self.unit_registry_object = UnitRegistry()

    def __call__(self, *args: Any, **kwargs: Any) -> OverrideQuantity:
        q_object = self.unit_registry_object(*args, **kwargs)
        return _make_override_quantity(q_object)

    def __getattr__(self, item: str) -> OverrideQuantity:
        if item == "h":
            item = "hour"

        q_object = 1 * self.unit_registry_object.__getattr__(item)
        return _make_override_quantity(q_object)


u = OverrideUnitRegistry()
