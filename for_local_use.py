from mobspy import *
import matplotlib.pyplot as plt


if __name__ == '__main__':
    A, B = BaseSpecies(2)
    A + B >> Zero[0.01]

    A(50), B(30)
    S2 = Simulation(A | B)
    S2.method = 'stochastic'
    S2.repetitions = 3
    S2.plot_data = False
    S2.duration = (A <= 0) | (B <= 0)
    S2.run()



