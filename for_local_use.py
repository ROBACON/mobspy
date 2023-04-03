from mobspy import *
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':

    A, B = BaseSpecies(2)
    A >> 2 * A[1]

    A(1)
    S1 = Simulation(A | B)
    S1.save_data = False
    S1.plot_data = True
    S1.duration = 3
    S1.level = -1

    A.reset_reactions()
    A + B >> Zero[0.01]

    B(50)
    S2 = Simulation(A | B)
    S2.duration = (A <= 0.5) | (B <= 0.5)
    S2.level = -1

    Sim = S1 + S2
    Sim.run()
    print(Sim.results[A][-1])
