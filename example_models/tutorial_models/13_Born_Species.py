from mobspy import *
""" 
    In this model we will explain the born species. This is how we refer to species that appear in the products
    but do not appear in the reactants
"""

# Let's look at the following mobspy model:
A = BaseSpecies(1)
B, C = New(A, 2)
B.b1, B.b2, C.c1, C.c2


Zero >> 2 * A[1]
# In the reaction above we add the meta-species A from which B, C inherit from
# However, A does not appear in the reactants
# MobsPy them does not have a paring, so no meta-species to pick from the reactant
# Here only two reactions are defined from this meta-reaction
# Which are the default states of all meta-species

MySim = Simulation(B | C)
print(MySim.compile())

# If one wants the born species meta-reaction to create reactions for all the species in the set,
# he can use the All default order in the following way:

A = BaseSpecies(1)
B, C = New(A, 2)
B.b1, B.b2, C.c1, C.c2
Zero >> 2 * A[1]
MySim = Simulation(B | C)
MySim.default_order = All
print(MySim.compile())