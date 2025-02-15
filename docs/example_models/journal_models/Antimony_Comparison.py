from mobspy import *

Replicator, Mortal = BaseSpecies()
Replicator >> 2*Replicator [lambda r: 2*(100 - r)*r]
Mortal >> Zero [1]

A, B, C = New(Replicator*Mortal)
A + B >> C [1]

S = Simulation(A | B | C)
print(S.generate_antimony()[0])


