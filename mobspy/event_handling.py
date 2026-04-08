"""
Event handling module: context managers and event compilation for MobsPy simulations.

Provides the EventHandlingMixin class that adds event_time(), event_condition(),
and related event management capabilities to the Simulation class.
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import TYPE_CHECKING, Any

from pint import Quantity

from mobspy.constants import STD_CHAR
from mobspy.dsl.reactions import Reacting_Species
from mobspy.dsl.species import Species
from mobspy.exceptions import EventError, ValidationError
from mobspy.types import SimulationEventData

if TYPE_CHECKING:
    from collections.abc import Generator


class EventHandlingMixin:
    """Mixin providing event context managers and compilation for Simulation."""

    # These attributes are provided by the Simulation class
    _event_time: float | int
    current_event_count_data: list[Any]
    total_packed_events: list[SimulationEventData]
    number_of_context_comparisons: int
    pre_number_of_context_comparisons: int
    _context_not_active: bool
    _conditional_event: bool

    @classmethod
    def event_compilation_error(cls) -> None:
        """Raise an EventError describing the expected condition format."""
        raise EventError(
            "The event condition did not compile.\n"
            "Please make sure it follows the following format:\n"
            "For simple conditions - if C1 \n"
            "For and based condition - if (C1) & (C2)\n"
            "For or based conditions - if (C1) | (C2)\n"
            "Please include the parentheses"
        )

    def event_context_finish(self) -> None:
        """Remove the context in all meta-species and reset variables."""
        self._event_time = 0
        Species.reset_simulation_context()
        self._context_not_active = True

    def event_context_add(self, time: float | int | Any, trigger: str) -> None:
        """Add an event to the event context.

        Args:
            trigger: Condition that triggers the event when fulfilled.
            time: Time to wait before triggering the event.
        """
        event_data = SimulationEventData(
            event_time=time,
            event_counts=list(self.current_event_count_data),
            trigger=trigger,
        )

        self.current_event_count_data = []
        self.pre_number_of_context_comparisons = self.number_of_context_comparisons
        self.number_of_context_comparisons = 0

        if event_data.event_counts:
            self.total_packed_events.append(event_data)

        self.event_context_finish()

    def event_context_initiator(self) -> None:
        """Set the context in all meta-species."""
        Species.set_simulation_context(self)  # type: ignore[arg-type]

    def _event_handler(self) -> None:
        """Activate the current context, checking it is the only one active."""
        if self._context_not_active:
            self._context_not_active = False
            self._set_parameter("_with_event", True)  # type: ignore[attr-defined]
            self.event_context_initiator()
        else:
            raise EventError("MobsPy does not support multiple context calls")

    @contextmanager
    def event_condition(
        self, trigger: str, delay: float | int | Quantity = 0
    ) -> Generator[int, None, None]:
        """Context manager for condition events.

        Used in ``with Simulation.event_condition(trigger):`` format to define
        events that trigger when a condition is met.

        Args:
            trigger: Condition string that triggers the event when fulfilled.
            delay: Time to wait before triggering (can be int, float, or Quantity).

        Yields:
            Always yields 0 for context manager compatibility.

        Raises:
            EventError: If invalid trigger syntax is used.
            ValidationError: If invalid trigger type is provided.
        """
        try:
            if isinstance(trigger, (bool, float, int)):
                raise ValidationError(
                    f"MobsPy has received an invalid trigger type: {type(trigger)}. "
                    "Please make sure you are not using the operator == "
                    "for creating event conditions"
                )

            self._conditional_event = True
            self._event_handler()
            yield 0
        finally:
            self._conditional_event = False
            delay_val = delay.magnitude if isinstance(delay, Quantity) else delay
            self.event_context_add(delay_val, trigger)

    @contextmanager
    def event_time(self, time: float | int | Quantity) -> Generator[int, None, None]:
        """Context manager for time events.

        Used in ``with Simulation.event_time(time):`` format to define
        events that trigger after a specified time.

        Args:
            time: Time delay before event triggers.

        Yields:
            Always yields 0 for context manager compatibility.
        """
        try:
            self._event_handler()
            yield 0
        finally:
            time_val = time.magnitude if isinstance(time, Quantity) else time
            self.event_context_add(time_val, "true")

    def at(
        self,
        time: float | int | Quantity,
        assignments: dict[Species | Reacting_Species, Any],
    ) -> None:
        """Define a time event without a context manager.

        Equivalent to::

            with S.event_time(t):
                A(count)
                B(count)

        Args:
            time: Time at which the event fires.
            assignments: ``{species_or_reacting: count}`` pairs.
        """
        self._set_parameter("_with_event", True)  # type: ignore[attr-defined]
        time_val = time.magnitude if isinstance(time, Quantity) else time
        event_counts = _resolve_event_assignments(assignments)
        event_data = SimulationEventData(
            event_time=time_val,
            event_counts=event_counts,
            trigger="true",
        )
        if event_counts:
            self.total_packed_events.append(event_data)

    def when(
        self,
        condition: Any,
        assignments: dict[Species | Reacting_Species, Any],
        delay: float | int | Quantity = 0,
    ) -> None:
        """Define a condition event without a context manager.

        Equivalent to::

            with S.event_condition(trigger, delay):
                A(count)

        Args:
            condition: Trigger condition (e.g. ``A <= 10``).
            assignments: ``{species_or_reacting: count}`` pairs.
            delay: Delay after condition becomes true.
        """
        if isinstance(condition, (bool, float, int)):
            raise ValidationError(
                f"Invalid trigger type: {type(condition)}. "
                "Do not use == for event conditions."
            )
        self._set_parameter("_with_event", True)  # type: ignore[attr-defined]
        delay_val = delay.magnitude if isinstance(delay, Quantity) else delay
        event_counts = _resolve_event_assignments(assignments)
        event_data = SimulationEventData(
            event_time=delay_val,
            event_counts=event_counts,
            trigger=condition,
        )
        if event_counts:
            self.total_packed_events.append(event_data)


def _resolve_event_assignments(
    assignments: dict[Species | Reacting_Species, Any],
) -> list[dict[str, Any]]:
    """Convert ``{species: count}`` dict into internal event count format."""
    event_counts: list[dict[str, Any]] = []
    for target, quantity in assignments.items():
        if isinstance(target, Reacting_Species):
            species_obj = target.list_of_reactants[0]["object"]
            chars = target.list_of_reactants[0]["characteristics"]
        elif isinstance(target, Species):
            species_obj = target
            chars = STD_CHAR
        else:
            raise ValidationError(
                f"Event assignment keys must be Species or Reacting_Species, "
                f"got {type(target)}"
            )
        event_counts.append(
            {
                "species": species_obj,
                "characteristics": chars,
                "quantity": quantity,
            }
        )
    return event_counts
