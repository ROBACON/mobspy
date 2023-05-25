from mobspy import *
import basico
import matplotlib.pyplot as plt
import os
import inspect
from mobspy.sbml_simulator.builder import build

if __name__ == '__main__':

    A = BaseSpecies()
    A.a1, A.a2
    a, b = MSParameters(1, 2)

    A >> 2*A [b]

    All[A](a)
    S = Simulation(A)
    S.duration = 5
    S.compile()






