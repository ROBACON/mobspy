"""Owned snapshots of DSL declarations, independent of simulation execution."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from mobspy.dsl.declarations import get_registry
from mobspy.dsl.list_species import List_Species
from mobspy.dsl.mobspy_parameters import Internal_Parameter_Constructor
from mobspy.dsl.species import Species
from mobspy.exceptions import ValidationError

if TYPE_CHECKING:
    from mobspy.dsl.reactions import Reactions


class Model:
    """Capture the declarations for a set of species.

    ``Model(A | B)`` can be reused by several simulations, including in another
    thread. Subsequent DSL declarations do not change its reaction/count
    snapshot. Species characteristics should be defined before capturing it.
    """

    def __init__(
        self,
        species: Species | List_Species,
        reactions: set[Reactions] | None = None,
    ) -> None:
        if not isinstance(species, (Species, List_Species)):
            raise ValidationError("Model requires Species joined with the | operator")
        selected = List_Species(species)  # type: ignore[arg-type]
        expanded = set(selected)
        for item in selected:
            expanded.update(item._linked_species)
        self.species = List_Species(expanded)
        ids = frozenset(
            id(reference)
            for item in self.species
            for reference in item.get_references()
        )
        self.declarations = get_registry().snapshot_for_species(ids)
        self.reactions = frozenset(
            reactions
            if reactions is not None
            else {
                reaction
                for item in self.species
                for reference in item.get_references()
                for reaction in reference.get_reactions()
            }
        )
        self.counts: list[dict[str, Any]] = [
            {
                "object": item,
                "characteristics": count["characteristics"].copy()
                if isinstance(count["characteristics"], set)
                else count["characteristics"],
                "quantity": count["quantity"],
            }
            for item in self.species
            for count in item.get_quantities()
        ]
        # String rates resolve names against the namespace at model creation,
        # rather than whichever thread happens to compile the model later.
        self.parameters = dict(Internal_Parameter_Constructor.parameter_stack)
