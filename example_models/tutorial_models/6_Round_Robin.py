from mobspy import *

"""
    In this example we will construct a model to explain the round-robin structure in detail
    The round-robin is a cyclic assignment structure that helps MobsPy deal with non-linear reactions
    We will use some examples to show it's importance
"""

# We define the species
Age = BaseSpecies(1)
Animal = New(Age)

# We give it characteristics
Animal.a1, Animal.a2

# We assign a age reaction
Age.young >> Age.old [1]

# Here we have a competition reaction
# The first impression is that the old animal kills the young animal
# HOWEVER, MobsPy uses a round-robin system where the first species in the product
# becomes the first-species in the reactants and so it continues in a cyclic fashion
# until all products have been assigned a species
# Here we actually have the young animal killing the old and becoming old afterwards
Animal.young + Animal.old >> Animal.old [1]

# To solve this we can either change the order or use .label()
Animal.young + Animal.old.label(1) >> Animal.old.label(1) [1]

# The cyclic round robin style is especially useful in reproduction reaction
# As the species keep the characteristics of the parent
Animal >> 2*Animal [1]

# We recommend commenting one of the competition reactions and running the model
# to visualize the round-robin system
MySim = Simulation(Animal)
print(MySim.compile())






