from mobspy import *


if __name__ == '__main__':

    A, B, C, D = BaseSpecies()
    B.b1, B.b2, B.b3
    #C.c1, C.c2
    #D.d1, D.d2

    with Assign:
        All[B](5)
    # exit()

    S = Simulation(A | B | C | D)
    print(S.compile())
    exit()


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

    A.assign(2*B.b1)

    S = Simulation(A | B)
    print(S.compile())
    exit()

    with Assign:
        print('Hi')
        r = (A + B)/B
        C.c1(r)
        print('Hey')
    exit()

    S = Simulation(A | B | C)
    print(S.compile())
