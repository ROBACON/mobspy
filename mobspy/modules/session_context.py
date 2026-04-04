"""Consolidated thread-local session state for the MobsPy DSL.

All per-thread context (registry, rate storage, event context, etc.)
lives in a single :class:`SessionContext` held in one ``ContextVar``.
Individual modules access their field through accessor functions
defined here.
"""

from __future__ import annotations

from contextvars import ContextVar
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from mobspy.modules.declarations import ModelRegistry


@dataclass
class SessionContext:
    """All per-thread DSL state in one place.

    Each field corresponds to a formerly independent ``ContextVar``.
    """

    # ModelRegistry (was _registry_cv in declarations.py)
    registry: ModelRegistry | None = field(default=None)

    # Auto-naming counter (was _entity_counter_cv in reactions.py)
    entity_counter: int = 0

    # Legacy rate buffering for [] syntax (was _last_rate_cv in reactions.py)
    last_rate: Any = None

    # Active Simulation for event context (was _simulation_context_cv in species.py)
    simulation_context: Any = None

    # Any context characteristics (was _meta_any_ctx_cv in species.py)
    meta_any_ctx: set[str] = field(default_factory=set)

    # Assignment mode flag (was _asg_context_cv in assignments_implementation.py)
    asg_context: bool = False

    # Any characteristic accumulation (was _any_chars_cv in any_species.py)
    any_chars: set[str] = field(default_factory=set)

    # Nested Any context stack (was _any_stack_cv in any_species.py)
    any_stack: list[set[str]] = field(default_factory=list)

    # Any context entry flag (was _any_building_cv in any_species.py)
    any_building: bool = False

    # Expression building mode (was _ms_active_ctx in species_operators.py)
    ms_active: bool = False

    # Parameter registry (was class-level parameter_stack with threading.Lock)
    parameter_stack: dict[str, Any] = field(default_factory=dict)


_session_cv: ContextVar[SessionContext] = ContextVar("_session_cv")


def get_session() -> SessionContext:
    """Return the thread-local SessionContext, creating one if needed."""
    try:
        return _session_cv.get()
    except LookupError:
        ctx = SessionContext()
        _session_cv.set(ctx)
        return ctx


def reset_session() -> None:
    """Reset all session state for the current thread.

    Used by test fixtures to prevent state leakage between tests.
    """
    try:
        _session_cv.get()
    except LookupError:
        return
    _session_cv.set(SessionContext())
