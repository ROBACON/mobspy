from mobspy import *
from mobspy import modules

if __name__ == '__main__':

    # Replace parameters using units
    A = BaseSpecies()
    k1 = ModelParameters(1/u.h)

    A >> Zero [k1]

    S = Simulation(A)
    S.compile()

    S.update_model([k1, 1/u.s])
    print(S.generate_sbml()[0])


