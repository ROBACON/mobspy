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

from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, TypeAlias

if TYPE_CHECKING:
    from pint import Quantity

    from mobspy.dsl.reactions import Reactions
    from mobspy.dsl.species import Species
    from mobspy.expressions.evaluation import OverrideQuantity
    from mobspy.expressions.nodes import ExprNode
    from mobspy.expressions.rate_builder import RateExpression

    # Rate values: numeric, string, callable (lambda replay), expression AST,
    # Quantity/OverrideQuantity, or tuple for reversible rates (k_fwd, k_rev).
    RateValue: TypeAlias = (
        int
        | float
        | str
        | Callable[..., Any]
        | RateExpression
        | ExprNode
        | Quantity
        | OverrideQuantity
        | tuple[Any, ...]
        | None
    )


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
    rate: RateValue = None
    is_reversible: bool = False
    reverse_rate: RateValue = None


@dataclass(frozen=True)
class CountAssignment:
    """Initial count assignment produced by ``Species(count)``."""

    species: Species
    characteristics: frozenset[str] | str
    quantity: int | float | Quantity


@dataclass(frozen=True)
class EventDecl:
    """Event declaration produced by ``S.at()`` or ``S.when()``."""

    trigger: str
    delay: float | int = 0
    assignments: tuple[EventAssignmentDecl, ...] = ()


@dataclass(frozen=True)
class EventAssignmentDecl:
    """Single assignment within an event."""

    species: Species
    characteristics: frozenset[str] | str
    quantity: int | float | Quantity


# ---------------------------------------------------------------------------
# RatedProduct: intermediate produced by the ``@`` operator
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RatedProduct:
    """Carries a rate alongside product species info.

    Created by ``B @ rate`` (via ``__matmul__``), consumed by
    ``A >> rated_product`` (via ``__rshift__``).

    Supports ``B + (C @ rate)`` via ``__radd__``: the left-hand
    species is prepended to the product list, returning a new
    ``RatedProduct``.  This makes ``@`` work naturally in
    multi-product reactions without requiring parentheses.
    """

    products: list[dict[str, Any]]
    rate: RateValue
    is_reversible: bool = False
    reverse_rate: RateValue = None

    def __radd__(self, other: Any) -> RatedProduct:
        """Support ``Species + RatedProduct``."""
        from mobspy.dsl.reactions import Reacting_Species  # noqa: PLC0415
        from mobspy.dsl.species import Species  # noqa: PLC0415  # circular import

        if isinstance(other, Species):
            other = Reacting_Species(other, set())
        if isinstance(other, Reacting_Species):
            merged = list(other.list_of_reactants) + self.products
            return RatedProduct(
                products=merged,
                rate=self.rate,
                is_reversible=self.is_reversible,
                reverse_rate=self.reverse_rate,
            )
        return NotImplemented

    def _adjust_rate(self, other: Any, op: str) -> RatedProduct:
        """Apply an arithmetic operation to the rate.

        Handles the Python precedence issue where ``B @ 0.01 / u.hour``
        parses as ``(B @ 0.01) / u.hour``. By supporting arithmetic on
        RatedProduct, both ``B @ (0.01 / u.hour)`` and
        ``B @ 0.01 / u.hour`` produce the same result.
        """
        if op == "mul":
            new_rate: Any = self.rate * other
        elif op == "truediv":
            new_rate = self.rate / other
        else:
            return NotImplemented  # type: ignore[no-any-return]  # dynamic dispatch
        rev = self.reverse_rate
        if rev is not None:
            if op == "mul":
                rev = rev * other
            elif op == "truediv":
                rev = rev / other
        return RatedProduct(
            products=self.products,
            rate=new_rate,
            is_reversible=self.is_reversible,
            reverse_rate=rev,
        )

    def __mul__(self, other: Any) -> RatedProduct:
        """Support ``(B @ rate) * scalar``."""
        return self._adjust_rate(other, "mul")

    def __rmul__(self, other: Any) -> RatedProduct:
        """Support ``scalar * (B @ rate)``."""
        return self._adjust_rate(other, "mul")

    def __truediv__(self, other: Any) -> RatedProduct:
        """Support ``(B @ rate) / unit``."""
        return self._adjust_rate(other, "truediv")


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
        self.reaction_objects: list[Reactions] = []
        self._counts: dict[tuple[int, frozenset[str] | str], CountAssignment] = {}
        self.events: list[EventDecl] = []
        self._reaction_species: dict[int, frozenset[int]] = {}

    @property
    def counts(self) -> list[CountAssignment]:
        """Return count assignments as an ordered list."""
        return list(self._counts.values())

    def add_reaction(
        self,
        decl: ReactionDecl,
        reaction_obj: Reactions | None = None,
    ) -> None:
        """Register a reaction declaration and its Reactions object."""
        self.reactions.append(decl)
        if reaction_obj is not None:
            self.reaction_objects.append(reaction_obj)
            species_ids: set[int] = set()
            for ref in decl.reactants:
                species_ids.add(id(ref.species))
            for ref in decl.products:
                species_ids.add(id(ref.species))
            self._reaction_species[id(reaction_obj)] = frozenset(species_ids)

    def add_count(self, assignment: CountAssignment) -> None:
        """Register or overwrite a count assignment.

        Keyed by ``(id(species), characteristics)`` so that
        reassigning a count updates the existing entry.
        """
        key = (id(assignment.species), assignment.characteristics)
        self._counts[key] = assignment

    def remove_counts_for(self, species: Species) -> None:
        """Remove all count assignments for a species.

        Called by ``Species.reset_quantities()`` to keep the
        registry in sync.
        """
        species_id = id(species)
        to_remove = [k for k in self._counts if k[0] == species_id]
        for k in to_remove:
            del self._counts[k]

    def remove_reactions_for(self, species: Species) -> None:
        """Remove all reactions where the given species participates.

        Called by ``Species.reset_reactions()`` to keep the
        registry in sync.
        """
        species_id = id(species)
        to_remove_ids: set[int] = set()
        for rxn_id, spe_ids in self._reaction_species.items():
            if species_id in spe_ids:
                to_remove_ids.add(rxn_id)

        # Filter both parallel lists together
        new_reactions: list[ReactionDecl] = []
        new_objects: list[Reactions] = []
        for decl, obj in zip(self.reactions, self.reaction_objects, strict=False):
            if id(obj) not in to_remove_ids:
                new_reactions.append(decl)
                new_objects.append(obj)
        self.reactions = new_reactions
        self.reaction_objects = new_objects
        for rxn_id in to_remove_ids:
            self._reaction_species.pop(rxn_id, None)

    def reactions_for_species(
        self,
        model_species_ids: frozenset[int],
    ) -> set[Reactions]:
        """Return reaction objects where at least one participant is in the set.

        Args:
            model_species_ids: Set of ``id(species)`` for all species in the
                model (including their references).
        """
        result: set[Reactions] = set()
        for rxn_obj in self.reaction_objects:
            spe_ids = self._reaction_species.get(id(rxn_obj), frozenset())
            if spe_ids & model_species_ids:
                result.add(rxn_obj)
        return result

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
        snap._reaction_species = dict(self._reaction_species)
        return snap

    def clear(self) -> None:
        """Reset all accumulated declarations."""
        self.reactions.clear()
        self.reaction_objects.clear()
        self._counts.clear()
        self.events.clear()
        self._reaction_species.clear()

    def snapshot_and_clear(self) -> ModelRegistry:
        """Snapshot then clear.  Used by Simulation.__init__."""
        snap = self.snapshot()
        self.clear()
        return snap


def get_registry() -> ModelRegistry:
    """Return the thread-local ModelRegistry, creating one if needed."""
    from mobspy.dsl.session_context import get_session  # noqa: PLC0415

    session = get_session()
    if session.registry is None:
        session.registry = ModelRegistry()
    return session.registry


def snapshot_registry() -> ModelRegistry:
    """Snapshot the current registry (does NOT clear it)."""
    return get_registry().snapshot()


def snapshot_and_clear_registry() -> ModelRegistry:
    """Snapshot and clear the current registry."""
    return get_registry().snapshot_and_clear()
