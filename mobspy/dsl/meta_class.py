"""Meta-species, meta-reactions, and related classes for MobsPy.

.. deprecated::
    Import directly from the sub-modules instead::

        from mobspy.dsl.species import Species
        from mobspy.dsl.reactions import Reactions

This module re-exports everything from the split sub-modules so that
existing imports like ``from mobspy.dsl.meta_class import Species``
continue to work unchanged.

The actual implementations live in:
- ``mobspy.dsl.species`` -- Species class
- ``mobspy.dsl.reactions`` -- Reactions, Reacting_Species, _Last_rate_storage
- ``mobspy.dsl.list_species`` -- List_Species
- ``mobspy.dsl.species_constructors`` -- BaseSpecies, New, Zero, One, EndFlagSpecies
"""

from __future__ import annotations

from mobspy.dsl.declarations import RatedProduct

# Re-export everything for backward compatibility
from mobspy.dsl.list_species import List_Species
from mobspy.dsl.reactions import (
    Assignment_Opp_Imp,
    Reacting_Species,
    Reactions,
    _Last_rate_storage,
    _methods_Reacting_Species,
)
from mobspy.dsl.species import (
    Species,
    _create_reaction_from_rated,
    _methods_Species,
    clean_species_name,
)
from mobspy.dsl.species_constructors import (
    BaseSpecies,
    EndFlagSpecies,
    ListSpecies,
    New,
    One,
    Zero,
    _Create_Species,
    compile_species_number_line,
)

__all__ = [
    "Assignment_Opp_Imp",
    "BaseSpecies",
    "EndFlagSpecies",
    "ListSpecies",
    "List_Species",
    "New",
    "One",
    "RatedProduct",
    "Reacting_Species",
    "Reactions",
    "Species",
    "Zero",
    "_Create_Species",
    "_Last_rate_storage",
    "_create_reaction_from_rated",
    "_methods_Reacting_Species",
    "_methods_Species",
    "clean_species_name",
    "compile_species_number_line",
]
