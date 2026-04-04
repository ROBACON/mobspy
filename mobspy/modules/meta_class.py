"""Meta-species, meta-reactions, and related classes for MobsPy.

.. deprecated::
    Import directly from the sub-modules instead::

        from mobspy.modules.species import Species
        from mobspy.modules.reactions import Reactions

This module re-exports everything from the split sub-modules so that
existing imports like ``from mobspy.modules.meta_class import Species``
continue to work unchanged.

The actual implementations live in:
- ``mobspy.modules.species`` -- Species class
- ``mobspy.modules.reactions`` -- Reactions, Reacting_Species, _Last_rate_storage
- ``mobspy.modules.list_species`` -- List_Species
- ``mobspy.modules.species_constructors`` -- BaseSpecies, New, Zero, One, EndFlagSpecies
"""

from __future__ import annotations

from mobspy.modules.declarations import RatedProduct

# Re-export everything for backward compatibility
from mobspy.modules.list_species import List_Species
from mobspy.modules.reactions import (
    Assignment_Opp_Imp,
    Reacting_Species,
    Reactions,
    _Last_rate_storage,
    _methods_Reacting_Species,
)
from mobspy.modules.species import (
    Species,
    _create_reaction_from_rated,
    _methods_Species,
    clean_species_name,
)
from mobspy.modules.species_constructors import (
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
