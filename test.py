from mobspy import *
import os


if __name__ == '__main__':

    A, B, C, D = BaseSpecies(4)
    Test = New(A, 1)

    Test.b(30) >> B [lambda: f'{Test}']
    Test.c(30) >> C [lambda: f'{Test}']
    Test.d(30) >> D [lambda: f'{Test}']
    Simulation(Test | B | C | D).run()
    exit()


