from mobspy import *

"""
    In this example we will construct a model to explain the round-robin structure in detail.
    The round-robin is a cyclic assignment structure that helps MobsPy deal with non-linear reactions.
    We will use some examples to show its importance.
"""

# We define the species:
Age = BaseSpecies(1)
Animal = New(Age)

# We assign characteristics to a species:
Animal.a1, Animal.a2

# We define an age reaction:
Age.young >> Age.old [1]

# The Round-Robin mechanism is demonstrated with the reaction bellow:
#   In this reaction we want to define a competition reaction between a young and an old animal.
#   We would like the old animal to survive the competition.
#
#   However, how does MobsPy distinguish this case from the case where the
#   young animal survives and becomes old?
#
#   In MobsPy the order of the meta-species in the reaction is important.
#   With the usage of a rond-robin system the first reactant becomes the first-species in the product
#   and it continues in a cyclic fashion until all products have been mapped to a species.
Animal.old + Animal.young >> Animal.old [1]
# Here, the old animal is the one that survives.

# The .label() can be used to circumvent this default behavior, if required.
# In the next reaction, we have the young Animal killing the old and becoming old itself.
Animal.old + Animal.young.label(1) >> Animal.old.label(1) [1]
# Here, the young animal is the one that survives.

# While arguably, the fact that once the old animal survives and once the young one
# does not make a difference in our setting, it does if additional characteristics are
# added to the Animal, e.g., it having a color.

# The cyclic round robin style is especially useful in reproduction reaction.
# Since the assignment cycles in a round-robin fashion, the first product will be assigned the first reactant
# and the second product will cycle through the reactants and be assigned the first reactant as well.
# This results, in the example below, the two products to have the same characteristics as the single reactant.
Animal >> 2*Animal [1]

# We recommend commenting out one of the competition reactions and running the model
# to play with the round-robin system.
MySim = Simulation(Animal)
print(MySim.compile())






