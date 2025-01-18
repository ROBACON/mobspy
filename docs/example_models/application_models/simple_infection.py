from mobspy import *

Age, Mortal, Infectable, Virus = BaseSpecies()
Reproducer = New(Age)
V1, V2 = New(Virus)

Age.young >> Age.old [1]

Reproducer >> 2*Reproducer.young[0.1]


def infection_rate(r1, r2):
    factor = 0.01
    factor = 2*factor if r1.old else 1*factor
    factor = 2*factor if r2.is_a(V2) else 1*factor
    return factor


Infectable.not_infected + Virus >> Infectable.infected [infection_rate]
Mortal >> Zero [lambda r1: 2 if r1.infected else 0.01]
Cell = Infectable * Mortal * Age * Reproducer

Cell(100)
V1(20), V2(25)
S = Simulation(Cell | V1 | V2)
S.duration = 5
S.method = 'stochastic'
S.repetitions = 5
S.run()

