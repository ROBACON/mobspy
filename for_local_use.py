from mobspy import *
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':

    A = BaseSpecies()

    Zero >> A[42]
    A >> Zero[1]

    S = Simulation(A)
    S.volume = 1*u.milliliter
    S.run()
    assert int(S.results[A][-1]) == 42

