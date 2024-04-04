from mobspy import *
import matplotlib.pyplot as plt

if __name__ == '__main__':

    A, B = BaseSpecies()

    A >> 2 * A[1.05 / u.s]
    B >> 2 * B[1 / u.s]

    A(1), B(1)
    S1 = Simulation(A | B)
    S1.duration = 3 * u.s

    A + B >> Zero[0.1 / u.s]

    S2 = Simulation(A | B)
    S2.duration = (A <= 0) | (B <= 0)
    S2.method = 'stochastic'

    S = S1 + S2
    S.plot_config.suptitle, S.plot_config.suptitle_fontsize = 'Mutual Annihilation Protocol', 18
    S.plot_config.xlabel_fontsize, S.plot_config.ylabel_fontsize = 14, 14
    S.plot_config.tight_layout = True
    S.plot_config.save_to = "testing_tight_layout.png"
    S.repetitions = 10 * u.s
    S.run()
    exit()

    A = BaseSpecies()

    A >> Zero [1]

    A(100)
    S = Simulation(A)
    S.duration = 10
    S.run(plot_data=False)
    fig, axis_matrix = S.plot_raw({}, return_fig=True)
    plt.plot()

