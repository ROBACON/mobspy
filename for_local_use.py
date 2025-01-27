from mobspy import *
from mobspy import modules

if __name__ == '__main__':

    # TODO: Ask Matthias about removed tests

    # Replace parameters using units
    A = BaseSpecies()
    k1, k2 = ModelParameters(1, 1)

    A >> Zero [1, 1]
    A >> Zero [1/(k1 + k2), k1]

    A >> 2*A [10]

    S = Simulation(A)
    print(S.compile())

    # print(A.get_characteristics())


