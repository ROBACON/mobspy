from mobspy import *
import basico
import matplotlib.pyplot as plt
import os
import inspect
from mobspy.sbml_simulator.builder import build

if __name__ == '__main__':

    A, B = BaseSpecies()
    a, b, c, d = ModelParameters([1, 2, 3], 1, 10, 1)

    A >> 2*A ['a + b']

    A(b)
    S = Simulation(A | B)

    with S.event_time(c):
        B(d)

    S.duration = 5
    S.compile()













