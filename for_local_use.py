from mobspy import *
import numpy as np

if __name__ == '__main__':

    A = BaseSpecies()
    A(10000)

    A >> Zero [1]

    S = Simulation(A)
    S.duration = 1000
    # S.plot_data = False
    S.plot_config.logscale = ["Y"]
    S.plot_config.y_filter = [0.01, 1e50]
    S.plot_config.x_from = [0, 2000]
    S.run()
    print(S.fres[A][-1])

