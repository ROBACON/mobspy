from mobspy import *
"""
    MobsPy can also add characteristics to species 
    and construct reactions based on this change of characteristics by a querying mechanism
    In this example we describe the Query Mechanism in detail 
"""

# First we define 2 BaseSpecies
Age, Color = BaseSpecies(2)

# We add characteristics to them using the dot operator
# Characteristics can be added explicitly - For instance the Color Meta-Species
# Or they can be added inside a reaction - For instance the Age Meta-Species
Color.blue, Color.red, Color.yellow
Age.young >> Age.old [1]


# We multiply Age and Color
# The new Meta-Species formed now inherits the characteristics from Color and Age
# But it inherits them as a two dimensional space where the characteristics sets are orthogonal
# Dummy now represents all the possible combinations of characteristics from the meta-species
# it inherits from separated by a dot
# As an example Dummy now represents the following species
# => Dummy.blue.young, Dummy.blue.old, Dummy.red.young, Dummy.red.old, ....
Dummy = Age*Color

# We can now query between the meta-species that Dummy represents
# This allows one to define a reactions for only a subset
# In the meta-reaction bellow only the species from Dummy with the old characteristic
# are assigned the reaction
Dummy.old >> Zero [1]

# We can also query over multiple dimensions, order is irrelevant
Dummy.young.red >> Dummy.blue [1]

# Generate the compile text
# Results also in the 3_Queries.txt
MySim = Simulation(Dummy)
print(MySim.compile())
