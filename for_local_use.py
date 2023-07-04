from mobspy import *
import basico
import matplotlib.pyplot as plt
import os
import inspect
from mobspy.sbml_simulator.builder import build
import timeit


if __name__ == '__main__':

    A = BaseSpecies()

    Zero >> A[lambda r1, r2, r3: 20]

    S = Simulation(A)
    S.compile()
    exit()

    '''
    start_time = timeit.default_timer()
    B = BaseSpecies()

    for i in range(5):
        B >> 2 * B [lambda r1: (i / u.h) * 1 / (1 + 10 / r1)]

    S = Simulation(B)
    S.level = -1
    S.compile()

    elapsed_2 = timeit.default_timer() - start_time

    start_time = timeit.default_timer()
    A = BaseSpecies()

    for i in range(5):
        A >> 2*A [f'{i}*1/(1 + 10/{A})']

    S = Simulation(A)
    S.level = -1
    S.compile()

    elapsed_1 = timeit.default_timer() - start_time

    print(elapsed_1, elapsed_2)
    exit()
    '''












