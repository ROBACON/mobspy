from __future__ import annotations

from typing import Any

from mobspy.modules.mobspy_expressions import (
    ExpressionDefiner as me_ExpressionDefiner,
)
from mobspy.modules.mobspy_expressions import (
    u,
)


def _iter_live_expressions() -> list[me_ExpressionDefiner]:
    """Return all live ExpressionDefiner instances from the registry."""
    live = []
    for ref in me_ExpressionDefiner._registry.values():
        obj = ref()
        if obj is not None:
            live.append(obj)
    return live


def expression_compilation_initiation() -> None:
    """Activate expression context on all live ExpressionDefiner instances.

    Uses the class-level registry instead of walking the call stack.
    """
    u._ms_active = True
    for expr in _iter_live_expressions():
        expr._ms_active = True


def expression_compilation_finish() -> None:
    """Deactivate expression context on all live ExpressionDefiner instances."""
    u._ms_active = False
    for expr in _iter_live_expressions():
        expr._ms_active = False


class Unit_Context_Setter:
    def __enter__(self) -> None:
        expression_compilation_initiation()

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: Any,
    ) -> None:
        expression_compilation_finish()
