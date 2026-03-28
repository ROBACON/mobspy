"""Meta-species, meta-reactions, and related classes for MobsPy.

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

# Re-export everything for backward compatibility
from mobspy.modules.list_species import List_Species as List_Species
from mobspy.modules.reactions import (
    Assignment_Opp_Imp as Assignment_Opp_Imp,
    Reactions as Reactions,
    Reacting_Species as Reacting_Species,
    _Last_rate_storage as _Last_rate_storage,
    _methods_Reacting_Species as _methods_Reacting_Species,
)
from mobspy.modules.species import (
    Species as Species,
    _get_multiline_code_context as _get_multiline_code_context,
    _methods_Species as _methods_Species,
    clean_species_name as clean_species_name,
)
from mobspy.modules.species_constructors import (
    BaseSpecies as BaseSpecies,
    EndFlagSpecies as EndFlagSpecies,
    ListSpecies as ListSpecies,
    New as New,
    One as One,
    Zero as Zero,
    _Create_Species as _Create_Species,
    compile_species_number_line as compile_species_number_line,
)

if __name__ == "__main__":
    pass
