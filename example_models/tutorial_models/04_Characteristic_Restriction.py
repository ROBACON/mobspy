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

# This model results in an error:
#   Each species contains a set of the characteristics added to it directly,
#   which means that is not obtained by inheritance
#   If two different meta-species in a simulation have an element in common in this set
#   the compiler will throw an error
#   In the error the name of the meta-species and the characteristic sets will be displayed
#
# To avoid this use inheritance to share characteristics or simply rename the characteristics
#

