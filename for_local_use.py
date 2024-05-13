from mobspy import *


if __name__ == '__main__':

    # C = BaseSpecies()
    # C.c1, C.c2, C.c3

    #def idk(r):
    #    print(r.get_state())
    #    return 1

    #C >> Zero [idk]

    #C(100)
    #S = Simulation(C)
    #S.compile()

    A, B, C = BaseSpecies()

    with Assign:
        r = (A + B)/B
        # C.c1(r)

    S = Simulation(A | B | C)
    print(S.compile())
