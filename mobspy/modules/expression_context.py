"""Manage the global expression-compilation context for unit-aware rate building."""

from __future__ import annotations

from typing import Any

from mobspy.modules.mobspy_expressions import _ms_active_ctx


def expression_compilation_initiation() -> None:
    """Activate expression context globally via context variable."""
    _ms_active_ctx.set(True)


def expression_compilation_finish() -> None:
    """Deactivate expression context globally via context variable."""
    _ms_active_ctx.set(False)


class Unit_Context_Setter:  # noqa: N801
    """Context manager for unit-aware compilation.

    Activates expression-building mode.
    """

    def __enter__(self) -> None:
        expression_compilation_initiation()

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: Any,
    ) -> None:
        expression_compilation_finish()
