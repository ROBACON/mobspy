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

# The Round-Robin mechanism can be seen through the reaction bellow:
#   In this reaction we want to define a competition reaction between an young and an old animal
#   We would like the old animal to survive the competition
#
#   However, how does MobsPy distinguish this case from the case where the
#   young animal survives and becomes old?
#
#   In MobsPy the order of the meta-species in the reaction is important
#   With the usage of a rond-robin system the first reactant becomes the first-species in the product
#   and it continues in a cyclic fashion until all products have been assigned a species
Animal.old + Animal.young >> Animal.old [1]

# The .label() can be used to circumvent this if necessary
# In this reaction we have the young Animal killing the old and becoming old himself
Animal.old + Animal.young.label(1) >> Animal.old.label(1) [1]

# The cyclic round robin style is especially useful in reproduction reaction
# Since the assignment cycles in a round-robin fashion the first product will be assigned the first reactant
# and the second product will cycle through the reactants and be assigned the first reactant as well
# Thus, allowing the products to have the same characteristics as the reactant
Animal >> 2*Animal [1]

# We recommend commenting one of the competition reactions and running the model
# to visualize the round-robin system
MySim = Simulation(Animal)
print(MySim.compile())






