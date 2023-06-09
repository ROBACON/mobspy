from mobspy import *
import basico
import matplotlib.pyplot as plt
import os
import inspect
from mobspy.sbml_simulator.builder import build

if __name__ == '__main__':


    exit()

    A = BaseSpecies()
    A.a1, A.a2
    a, b, c, d, f, h = MSParameters([1, 2], [1, 2], [1, 2], [1, 2], [1, 2], [1, 2])

    A >> 2*A [lambda: f'5*(b + c)/10']

    All[A](a)
    S1 = Simulation(A)

    with S1.event_time(0):
        A.a2(d)

    with S1.event_time(2):
        A.a1('a + b')

    with S1.event_time(f):
        A.a1(d)

    S1.duration = 3

    B = BaseSpecies()

    B >> 2*B [h]

    B(a)
    S2 = Simulation(A | B)
    S2.duration = 2

    S = S1 + S2
    S.run()






