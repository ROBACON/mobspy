from mobspy import *
"""
    Restraints
    MobsPy has some restraints for modeling regarding the characteristics combinations
    In this model we will discuss them
"""

try:
    # Here we define two base Species
    A, C = BaseSpecies(2)

    # We them add two characteristics to characteristics to each of them
    A.aaa, A.bbb
    C.aaa, C.ccc

    # And add them to the simulation object
    MySim = Simulation(A | C)
    MySim.level = 0
    MySim.compile()
except SystemExit:
    pass

# This model results in an error
# Each species contains a set of the characteristics added to it directly (Without considering inheritance)
# If two species share an element inside this set MobsPy will throw an error
# And specify the two species that have the common characteristic
# Characteristics must be shared through inheritance not by directly adding them
# This keeps the characteristic state-space independent

