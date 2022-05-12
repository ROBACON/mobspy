from mobspy import *
"""
    This model is designed to explain the count assignment process in detail
"""

# Start with the base species
Age, Size = BaseSpecies(2)

# Define some simple reactions
Age.young >> Age.old [1]
Size.small >> Size.medium [1]
Size.medium >> Size.big [1]

# Define a species through multiplication
Thing = Age*Size

# Now the meta-species Thing has 6 species
# If we use the call operator in Things it will assign the count to the default state
# The default state is formed by only the first characteristic that added to each meta-species
# So the count will be assigned to Thing.small.young
Thing(100)

# If one wants to assign a count to another species, just use the dot operator
# Every characteristic not specified will be replaced by the default
Thing.medium(100)
Thing.big.old(100)

MySim = Simulation(Thing)
print(MySim.compile())


