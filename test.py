from mobspy import *
import os


if __name__ == '__main__':

    A, B = BaseSpecies(2)
    A.young, A.old
    A >> B [lambda r: 5 if r.young else 1]
    Simulation(A | B).compile()
    exit()

    A, B, C, D = BaseSpecies(4)
    Test = New(A, 1)

    Test.b(30) >> B [lambda: f'{Test}']
    Test.c(30) >> C [lambda: f'{Test}']
    Test.d(30) >> D [lambda: f'{Test}']
    A = Simulation(Test | B | C | D)
    exit()


