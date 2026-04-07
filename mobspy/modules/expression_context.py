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


def set_compilation_model_context(model_context: Any) -> None:
    """Store the model context on the session for expression unit conversion."""
    from mobspy.modules.session_context import get_session  # noqa: PLC0415

    get_session().model_context = model_context


def get_compilation_model_context() -> Any:
    """Retrieve the model context from the session."""
    from mobspy.modules.session_context import get_session  # noqa: PLC0415

    return get_session().model_context


class Unit_Context_Setter:  # noqa: N801
    """Context manager for unit-aware compilation.

    Activates expression-building mode and stores the model context.
    """

    def __init__(self, model_context: Any = None) -> None:
        self._model_context = model_context

    def __enter__(self) -> None:
        expression_compilation_initiation()
        if self._model_context is not None:
            set_compilation_model_context(self._model_context)

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: Any,
    ) -> None:
        expression_compilation_finish()
        set_compilation_model_context(None)
