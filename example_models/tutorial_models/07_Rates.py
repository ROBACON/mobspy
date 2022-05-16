from mobspy import *

"""
    In this model we discuss rate assignment possibilities in detail.
"""

# We start with the base species
# and add some characteristics:
Position, HandShaker, Lifter, Box = BaseSpecies(4)
Lifter.weak, Lifter.strong
Position.p1, Position.p2, Position.p3

# The rate function arguments work in the following manner:
#   The rate function returns a value based on the characteristics of the reactants.
#   The value may thus differ for the species in the meta-species.
#   A meta-reaction creates a reaction for all combinations of species.
#   The species that are currently being used to construct a reaction
#   are passed as arguments to the rate function.
#
# Now we want to define a box breaking reaction:
#   Here, we wish to say that strong lifters break the box faster.
#   In this case r1 will receive a species from the Lifter meta-species as it is the first
#   to appear in the reaction.
#   If another argument (say, r2) were passed to the function, MobsPy would pass species from the Box meta-species.
#   This is optional, however.
#   Arguments can have any name.
#   The rate function may filter on characteristics: The dot notation returns true iff the species has that
#   characteristic, and false otherwise.
#   We use this to make sure that only the species from the meta-species lifter that have the strong characteristics
#   will receive the higher rate.
Lifter + Box >> Lifter [lambda r1: 10 if r1.strong else 1]
# The value returned by the rate function is interpreted as the reaction rate constant of a
# rate that follows mass-action kinetics.

# The rate function can return strings as well.
# This can be used if kinetics different from mass-action kinetics are to be specified.
Lifter >> Lifter [lambda r1: f'1/(1+{r1})']

# Finally, one may wish to check the value of a certain characteristic to define a rate.
# For example, if two Position species are at the same place, we want them to become a Handshake.
# This can be done by using the call argument with the meta-species the characteristic has been added to.
Position + Position >> HandShaker [lambda r1, r2: 100 if Position(r1) == Position(r2) else 0]

MySim = Simulation(Position | Lifter | Box | HandShaker)
print(MySim.compile())





