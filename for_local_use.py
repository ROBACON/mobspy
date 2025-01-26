from mobspy import *
from mobspy import modules

if __name__ == '__main__':

    # TODO: Ask Matthias about removed tests

    # Replace parameters using units
    A = BaseSpecies()
    A.a1, A.a2
    B = New(A)
    B.b1, B.b2
    k1 = ModelParameters(1)

    B >> Zero [k1]

    B(100), B.b2(100)
    S = Simulation(B)
    S.compile()

    S.update_model([B, 200/u.l],[B.b2, 300/u.l])

    print(S.generate_sbml()[0])

    # print(A.get_characteristics())


