from mobspy import *
from tellurium import loada as te_load_anti

if __name__ == '__main__':

    A, TestSpe = BaseSpecies()
    a = ModelParameters([1, 2])

    R1 = A + TestSpe >> Zero [a]
    R2 = A >> 2*A [0.01]

    A(2), TestSpe(1)
    S = Simulation(A | TestSpe)

    S = S - R1
    S.duration = 10
    print(S.compile())


