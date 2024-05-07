from mobspy import *


if __name__ == '__main__':

    # TODO: Fix All in Assignments

    A, B, P = BaseSpecies()

    B >> 2*B [lambda r: (100-r)*r]
    B >> Zero [1]

    A.assign(B**2/(B**2 + 50))

    A >> A + P [1]
    P >> Zero [1]

    S = Simulation(A | B | P)
    S.run()
    print(S.compile())


