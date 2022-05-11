from mobspy import *
"""
    MobsPy can also add characteristics to species 
    and construct reactions based on this change of characteristics by a querying mechanism
    In this example we describe the Query Mechanism in detail 
"""

# First we define 2 BaseSpecies.
Age, Color = BaseSpecies(2)

# Characteristics:
#   We add characteristics to them using the dot operator.
#   Characteristics can be added explicitly or implicitly.
#   For instance the Color  meta-species uses explicit adding.
Color.blue, Color.red, Color.yellow
#   Alternatively, they can be added implicitly within the reaction definition.
#   An example for the latter is the Age  meta-species.
Age.young >> Age.old [1]


# We next multiply Age and Color.
# The new  meta-species inherits the characteristics from Color and Age.
# It inherits them as a two dimensional space where the characteristics sets are orthogonal.
# The Dummy meta-species represents all the possible combinations of characteristics from the meta-species
# it inherits from, separated by a dot.
# As an example, Dummy now represents the following species:
#   Dummy.blue.young, Dummy.blue.old, Dummy.red.young, Dummy.red.old, ....
Dummy = Age*Color

# Queries:
#   We can now query the meta-species that Dummy represents.
#   This allows one to define a reactions for a subset of species in Dummy
#   that match the query.
#   In the meta-reaction bellow, only the species from Dummy with the old characteristic
#   are assigned the reaction
Dummy.old >> Zero [1]

# We can also query over multiple characteristics with the order being irrelevant.
Dummy.young.red >> Dummy.blue [1]

# To generate the compiled string, we execute
MySim = Simulation(Dummy)
print(MySim.compile())
# Results are available in 3_Queries.txt
