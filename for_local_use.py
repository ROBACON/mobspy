from mobspy import *
from mobspy import modules

if __name__ == '__main__':

    # TODO: Ask Matthias about removed tests

    # Replace parameters using units
    A = BaseSpecies()
    A.a1, A.a2
    k1 = ModelParameters(1)

    A >> Zero [k1]

    A(100)
    S = Simulation(A)
    S.level=-1
    S.compile()

    S.update_model([A, 200/u.l])
    print(S.generate_sbml()[0])

    # print(A.get_characteristics())


