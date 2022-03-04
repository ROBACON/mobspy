from mobspy import *
import os


if __name__ == '__main__':

    A, B, C, D = BaseSpecies(4)
    A(20) + B(10) >> 2*C + D [0.1]
    MySim = Simulation(A | B | C | D)
    MySim.duration = 100
    MySim.save_data = False
    MySim.plot_data = False
    MySim.run()
    MySim.plot_deterministic(A, C)
    exit()


