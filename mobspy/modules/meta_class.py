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
)
from mobspy.modules.reactions import (
    Reacting_Species as Reacting_Species,
)
from mobspy.modules.reactions import (
    Reactions as Reactions,
)
from mobspy.modules.reactions import (
    _Last_rate_storage as _Last_rate_storage,
)
from mobspy.modules.reactions import (
    _methods_Reacting_Species as _methods_Reacting_Species,
)
from mobspy.modules.species import (
    Species as Species,
)
from mobspy.modules.species import (
    _get_multiline_code_context as _get_multiline_code_context,
)
from mobspy.modules.species import (
    _methods_Species as _methods_Species,
)
from mobspy.modules.species import (
    clean_species_name as clean_species_name,
)
from mobspy.modules.species_constructors import (
    BaseSpecies as BaseSpecies,
)
from mobspy.modules.species_constructors import (
    EndFlagSpecies as EndFlagSpecies,
)
from mobspy.modules.species_constructors import (
    ListSpecies as ListSpecies,
)
from mobspy.modules.species_constructors import (
    New as New,
)
from mobspy.modules.species_constructors import (
    One as One,
)
from mobspy.modules.species_constructors import (
    Zero as Zero,
)
from mobspy.modules.species_constructors import (
    _Create_Species as _Create_Species,
)
from mobspy.modules.species_constructors import (
    compile_species_number_line as compile_species_number_line,
)
