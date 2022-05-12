from mobspy import *

"""
    In this model we discuss rate assignment possibilities in detail
"""

# Start with the base species
# And Add some characteristics
Position, HandShaker, Lifter, Box = BaseSpecies(4)
Lifter.weak, Lifter.strong
Position.p1, Position.p2, Position.p3

# The rate function arguments work in the following manner:
#   The rate function will return a different value based on the characteristics of the reactant
#   Since meta-species are actually a set of species
#   A meta-reaction involves several reactions between all combinations of those species
#   And at every defined reaction, one of the species from the set will be used
#   Finally, the species that are currently being used to construct one reaction from the meta-reaction
#   will be passed as arguments to a rate function
#
# Now we want to define a box breaking reaction:
#   Here we wish to say that strong lifters break the box faster
#   In this case r1 will receive a species from the Lifter meta-species as it is the first
#   to appear in the reaction
#   If another argument was passed (r2) to the function MobsPy would pass species from the Box meta-species
#   The arguments can have any name you wish
#   With the dot notation they return true if the species has that characteristic and false otherwise
#   In this way only the species from the meta-species lifter that have the strong characteristics will receive
#   the higher rate
Lifter + Box >> Lifter [lambda r1: 10 if r1.strong else 1]

# The function can return strings as well
# This is the case where the user wishes to use something different than mass action kinetics
Lifter >> Lifter [lambda r1: f'1/(1+{r1})']

# Finally one may wish to check the value of a certain characteristic to program their rates
# Here when two Position species are at the same place to become a Handshake
# This can be done by using the call argument with the meta-species the characteristic has been added to
Position + Position >> HandShaker [lambda r1, r2: 100 if Position(r1) == Position(r2) else 0]

MySim = Simulation(Position | Lifter | Box | HandShaker)
print(MySim.compile())





