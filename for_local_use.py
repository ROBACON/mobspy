from mobspy import *
from mobspy import modules

if __name__ == '__main__':

    # Replace parameters using units
    A = BaseSpecies()
    A.a1, A.a2
    k1 = ModelParameters(1)

    A >> Zero [k1]

    A(100)
    S = Simulation(A)
    S.compile()

    S.update_model([A, 200])
    print(S.generate_sbml()[0])


