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
#   Each species contains a set of the characteristics added to it directly
#   (i.e., characteristics are not obtained by inheritance).
#   If in a simulation two different meta-species contain the same element,
#   the compiler will throw an error.
#   The error message contains the name of the meta-species and the sets of species they cprrespond to.
#
# To avoid this error, use inheritance to share characteristics, or simply rename the characteristics.
#

