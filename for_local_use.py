from mobspy import *
import matplotlib.pyplot as plt

if __name__ == '__main__':

    A = BaseSpecies()

    A >> 2*A [lambda r: (1000 - r)*r]
    A >> Zero [1]

    A(1)
    S = Simulation(A)
    S.duration = 10
    S.run()


