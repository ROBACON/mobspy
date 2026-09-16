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
from weakref import ReferenceType, WeakValueDictionary, ref

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
    """Weak discovery index for declarations owned by their species.

    Captured models use owning snapshots; the ambient session cannot keep an
    otherwise unreachable model alive. Removing a simulation never edits this
    index or the declarations of another simulation.
    """

    def __init__(self, *, owning: bool = False) -> None:
        self._objects: dict[int, Reactions] | WeakValueDictionary[int, Reactions] = (
            {} if owning else WeakValueDictionary()
        )
        self._counts: dict[
            tuple[int, frozenset[str] | str],
            tuple[ReferenceType[Species], Any],
        ] = {}
        self.events: list[EventDecl] = []
        self._owners: list[Species] | None = [] if owning else None

    @property
    def reaction_objects(self) -> list[Reactions]:
        return list(self._objects.values())

    @property
    def reactions(self) -> list[ReactionDecl]:
        return [reaction.declaration for reaction in self._objects.values()]

    @property
    def counts(self) -> list[CountAssignment]:
        result = []
        for (_, chars), (reference, quantity) in list(self._counts.items()):
            species = reference()
            if species is not None:
                result.append(CountAssignment(species, chars, quantity))
        return result

    def add_reaction(self, decl: ReactionDecl, reaction_obj: Reactions) -> None:
        reaction_obj.declaration = decl
        self._objects[id(reaction_obj)] = reaction_obj

    def add_count(self, assignment: CountAssignment) -> None:
        key = (id(assignment.species), assignment.characteristics)
        self._counts[key] = (
            ref(assignment.species, lambda _: self._counts.pop(key, None)),
            assignment.quantity,
        )
        if self._owners is not None:
            self._owners.append(assignment.species)

    def remove_counts_for(self, species: Species) -> None:
        for key in list(self._counts):
            if key[0] == id(species):
                del self._counts[key]

    def remove_reactions_for(self, species: Species) -> None:
        for reaction in species.get_reactions():
            self._objects.pop(id(reaction), None)
            for participant in (
                *reaction.declaration.reactants,
                *reaction.declaration.products,
            ):
                participant.species._reactions.discard(reaction)

    def reactions_for_species(
        self, model_species_ids: frozenset[int]
    ) -> set[Reactions]:
        return {
            reaction
            for reaction in self._objects.values()
            if any(
                id(participant.species) in model_species_ids
                for participant in (
                    *reaction.declaration.reactants,
                    *reaction.declaration.products,
                )
            )
        }

    def add_event(self, event: EventDecl) -> None:
        self.events.append(event)

    def snapshot(self) -> ModelRegistry:
        snapshot = ModelRegistry(owning=True)
        for reaction in self.reaction_objects:
            snapshot.add_reaction(reaction.declaration, reaction)
        for assignment in self.counts:
            snapshot.add_count(assignment)
        snapshot.events = list(self.events)
        return snapshot

    def snapshot_for_species(self, species_ids: frozenset[int]) -> ModelRegistry:
        snapshot = ModelRegistry(owning=True)
        for reaction in self.reactions_for_species(species_ids):
            snapshot.add_reaction(reaction.declaration, reaction)
        for assignment in self.counts:
            if id(assignment.species) in species_ids:
                snapshot.add_count(assignment)
        return snapshot

    def clear(self) -> None:
        self._objects.clear()
        self._counts.clear()
        self.events.clear()
        if self._owners is not None:
            self._owners.clear()

    def snapshot_and_clear(self) -> ModelRegistry:
        snapshot = self.snapshot()
        self.clear()
        return snapshot


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
