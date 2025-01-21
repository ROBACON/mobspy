from mobspy import *
from mobspy import modules

if __name__ == '__main__':

    k1 = ModelParameters(1)
    A, B, C, D = BaseSpecies()

    k1 = ModelParameters([1, 2, 3])

    A >> Zero [2*k1]

    A(100)
    S1 = Simulation(A)
    S1.duration = 10
    # print(S.compile())
    # print(S.generate_sbml()[0])

    A.reset_reactions()
    B >> Zero [1]

    B(200)
    S2 = Simulation(A | B)
    S2.duration = 5

    S = S1 + S2
    S.compile()
    # S.run()
    S1.update_model([k1, 1])
    S.run()

    # TODO - Fix 2*k1 in B later
    # TODO - Fix multiple compilation issue
