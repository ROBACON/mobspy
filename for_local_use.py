from mobspy import *
from tellurium import loada as te_load_anti

if __name__ == '__main__':

    A, TestSpe = BaseSpecies()
    a = ModelParameters([1, 2])

    A + TestSpe >> Zero [a]
    A >> 2*A [0.01]

    A(2), TestSpe(1)
    S = Simulation(A | TestSpe)
    S.duration = 10

    anti_model = S.generate_antimony()[0]
    r = te_load_anti(anti_model)
    results = r.simulate(0, 10, 3)
    assert results[0][1] == 2
    assert results[-1][1] > 1
    assert results[-1][2] < 0.01


