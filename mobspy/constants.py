"""Internal sentinel strings and magic constants used throughout MobsPy.

These constants replace hard-coded magic strings to improve maintainability,
enable IDE navigation, and prevent typo-related bugs.

Examples:
    >>> from mobspy.constants import DOT_SEPARATOR, ALL_CHAR
    >>> DOT_SEPARATOR
    '_dot_'
    >>> ALL_CHAR
    'all$'
"""

from __future__ import annotations

# Separator used in species names for SBML compatibility (dots are invalid in SBML IDs)
DOT_SEPARATOR: str = "_dot_"

# Characteristic markers used in the compiler and reaction system
ALL_CHAR: str = "all$"
STD_CHAR: str = "std$"
NOT_CHAR: str = "not$"

# Expression mode prefixes for species references in rate expressions
COUNT_PREFIX: str = "$count$"
CONCENTRATION_PREFIX: str = "$concentration$"
ASSIGNMENT_PREFIX: str = "$asg_"

# Null species placeholder in expressions
NULL_SPECIES: str = "$Null"

# SBML parameter tracking location marker
SBML_LOCATION: str = "$sbml"

# Pre-defined singleton species names
ZERO_SPECIES_NAME: str = "_S0"
ONE_SPECIES_NAME: str = "_S1"
END_FLAG_SPECIES_NAME: str = "_End_Flag_MetaSpecies"

# Context meta-species markers
CONTEXT_SPECIES_NAME: str = "Context_MetaSpecies"
CONTEXT_ANY_SPECIES_NAME: str = "Context_Any_MetaSpecies"

# Auto-generated name prefixes
RATE_NAME_PREFIX: str = "N$"
SET_SPECIES_PREFIX: str = "_set_spe_"
POSITION_PREFIX: str = "$_pos_"

# Statistical data suffixes for plot processing
AVERAGE_SUFFIX: str = "$average"
DEVIATION_SUFFIX: str = "$deviation"
