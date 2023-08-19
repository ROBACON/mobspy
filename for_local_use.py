from mobspy import *
import numpy as np
import basico
import matplotlib.pyplot as plt
import os
import inspect
from mobspy.sbml_simulator.builder import build
import timeit


if __name__ == '__main__':

    from mobspy import *

    A, B = BaseSpecies()

    A >> 2 * A[1.02]
    B >> 2 * B[1]

    A(1), B(1)
    S1 = Simulation(A | B)
    S1.duration = 3

    A + B >> Zero[0.1]

    S2 = Simulation(A | B)
    S2.duration = (A <= 0) | (B <= 0)
    S2.method = 'stochastic'

    S = S1 + S2
    S.plot_config.suptitle, S.plot_config.suptitle_fontsize = 'A vs B Algorithm', 16
    S.repetitions = 5
    S.run()

    A = [(0, 0), (0, 0), (0, 0)]
    B = A[1:]

    for (a, b), (c, d) in zip(A, B):
        print(a, b, c, d)

    exit()

    # Test script
    x = MobsPyExpression('A', None, dimension=3,
                         count_in_model=True,
                         concentration_in_model=False,
                         count_in_expression=False,
                         concentration_in_expression=False)

    r = x ** 2.8

























