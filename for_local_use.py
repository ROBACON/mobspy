from mobspy import *
import matplotlib.pyplot as plt


if __name__ == '__main__':

    A, B = BaseSpecies(2)
    A >> 2 * A[1]

    A(1)
    S1 = Simulation(A)
    S1.save_data = False
    # S1.plot_data = False
    S1.duration = A >= 70

    A.reset_reactions()
    A + B >> Zero[0.01]

    B(50)
    S2 = Simulation(A | B)
    S2.method = 'stochastic'
    S2.duration = (A <= 0) | (B <= 0)

    Sim = S1 + S2
    Sim.run()
    print(Sim.results[A][-1])
    print(Sim.results[B][-1])
    print(Sim.results[A][-1] == 0 or Sim.results[B][-1] == 0)





