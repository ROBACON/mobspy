"""Immutable declaration types produced by DSL operators.

These types decouple the DSL surface (operator overloading) from the
internal compilation pipeline.  Operators produce declarations that
accumulate in a thread-local ModelRegistry; the compiler reads
declarations instead of mutating Species objects.

ContextVar usage summary
------------------------
MobsPy uses ``ContextVar`` for thread-safe DSL state.  The refactored
architecture reduces reliance on implicit context:

=============================  =========  ==============================
ContextVar                     Status     Notes
=============================  =========  ==============================
``_last_rate_cv``              Legacy     Only needed for ``[]`` syntax.
                                          The ``@`` operator bypasses it
                                          via ``RatedProduct``.
``_entity_counter_cv``         Active     Generates unique species names.
``_simulation_context_cv``     Active     Tracks active event context.
                                          ``S.at()``/``S.when()`` bypass
                                          it entirely.
``_meta_any_ctx_cv``           Active     ``with Any.char:`` context.
``_asg_context_cv``            Active     Assignment mode flag.
``_registry_cv``               New        Thread-local ModelRegistry.
=============================  =========  ==============================
"""

from __future__ import annotations

from contextvars import ContextVar
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from mobspy.modules.species import Species


# ---------------------------------------------------------------------------
# Declaration value objects (immutable)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReactantRef:
    """Reference to a species on one side of a reaction."""

    species: Species
    characteristics: frozenset[str]
    stoichiometry: int | float = 1
    label: int | float | str | None = None


@dataclass(frozen=True)
class ReactionDecl:
    """A single reaction declaration produced by ``>>``.

    Stored in the ModelRegistry and later read by the compiler.
    """

    reactants: tuple[ReactantRef, ...]
    products: tuple[ReactantRef, ...]
    rate: Any = None
    is_reversible: bool = False
    reverse_rate: Any = None


@dataclass(frozen=True)
class CountAssignment:
    """Initial count assignment produced by ``Species(count)``."""

    species: Species
    characteristics: frozenset[str] | str
    quantity: Any


@dataclass(frozen=True)
class EventDecl:
    """Event declaration produced by ``S.at()`` or ``S.when()``."""

    trigger: str
    delay: float | int | Any = 0
    assignments: tuple[EventAssignmentDecl, ...] = ()


@dataclass(frozen=True)
class EventAssignmentDecl:
    """Single assignment within an event."""

    species: Species
    characteristics: frozenset[str] | str
    quantity: Any


# ---------------------------------------------------------------------------
# RatedProduct: intermediate produced by the ``@`` operator
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RatedProduct:
    """Carries a rate alongside product species info.

    Created by ``B @ rate`` (via ``__rmatmul__``), consumed by
    ``A >> rated_product`` (via ``__rshift__``).
    """

    products: list[dict[str, Any]]
    rate: Any
    is_reversible: bool = False
    reverse_rate: Any = None


# ---------------------------------------------------------------------------
# Model registry (thread-local accumulator)
# ---------------------------------------------------------------------------


class ModelRegistry:
    """Thread-local accumulator for DSL declarations.

    Operators append declarations here.  ``Simulation.__init__``
    snapshots and clears the registry.
    """

    def __init__(self) -> None:
        self.reactions: list[ReactionDecl] = []
        self.reaction_objects: list[Any] = []
        self._counts: dict[tuple[int, Any], CountAssignment] = {}
        self.events: list[EventDecl] = []

    @property
    def counts(self) -> list[CountAssignment]:
        """Return count assignments as an ordered list."""
        return list(self._counts.values())

    def add_reaction(
        self,
        decl: ReactionDecl,
        reaction_obj: Any = None,
    ) -> None:
        """Register a reaction declaration and its Reactions object."""
        self.reactions.append(decl)
        if reaction_obj is not None:
            self.reaction_objects.append(reaction_obj)

    def add_count(self, assignment: CountAssignment) -> None:
        """Register or overwrite a count assignment.

        Keyed by ``(id(species), characteristics)`` so that
        reassigning a count updates the existing entry.
        """
        key = (id(assignment.species), assignment.characteristics)
        self._counts[key] = assignment

    def remove_counts_for(self, species: Any) -> None:
        """Remove all count assignments for a species.

        Called by ``Species.reset_quantities()`` to keep the
        registry in sync.
        """
        species_id = id(species)
        to_remove = [k for k in self._counts if k[0] == species_id]
        for k in to_remove:
            del self._counts[k]

    def add_event(self, event: EventDecl) -> None:
        """Register an event declaration."""
        self.events.append(event)

    def snapshot(self) -> ModelRegistry:
        """Return a shallow copy of the current state."""
        snap = ModelRegistry()
        snap.reactions = list(self.reactions)
        snap.reaction_objects = list(self.reaction_objects)
        snap._counts = dict(self._counts)
        snap.events = list(self.events)
        return snap

    def clear(self) -> None:
        """Reset all accumulated declarations."""
        self.reactions.clear()
        self.reaction_objects.clear()
        self._counts.clear()
        self.events.clear()

    def snapshot_and_clear(self) -> ModelRegistry:
        """Snapshot then clear.  Used by Simulation.__init__."""
        snap = self.snapshot()
        self.clear()
        return snap


_registry_cv: ContextVar[ModelRegistry] = ContextVar("_registry_cv")


def get_registry() -> ModelRegistry:
    """Return the thread-local ModelRegistry, creating one if needed."""
    try:
        return _registry_cv.get()
    except LookupError:
        reg = ModelRegistry()
        _registry_cv.set(reg)
        return reg


def snapshot_registry() -> ModelRegistry:
    """Snapshot the current registry (does NOT clear it)."""
    return get_registry().snapshot()


def snapshot_and_clear_registry() -> ModelRegistry:
    """Snapshot and clear the current registry."""
    return get_registry().snapshot_and_clear()
