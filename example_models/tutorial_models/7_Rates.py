from mobspy import *

"""
    In this model we discuss rate assignment possibilities in detail
"""

# Start with the base species
# And Add some characteristics
Position, HandShaker, Lifter, Box = BaseSpecies(4)
Lifter.weak, Lifter.strong
Position.p1, Position.p2, Position.p3

# Firstly we want to define a box breaking reaction
# But we wish to say that strong lifters break the box faster
# To do so we define a function
# MobsPy will pass the species from each meta-species as an argument for the function
# To check if the species passed has the characteristic wanted just use the dot notation
# It will return True when it does and False otherwise
# The arguments can be named anything but the first argument will receive the first reactant and so on
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





