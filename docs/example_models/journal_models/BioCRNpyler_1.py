from mobspy import *

A, B, C, D = BaseSpecies()

A >> 2*B [3]
B >> C + D [1.4]

S = Simulation(A | B | C | D)
print(S.compile())

