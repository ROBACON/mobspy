from mobspy import *
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':

    A = BaseSpecies()
    A >> Zero [1]

    S = Simulation(A)

    with S.event_time(0):
        Zero >> A [1]







