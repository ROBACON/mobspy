from mobspy import *
import os

A, B = BaseSpecies()

# Replication reactions
A >> 2 * A [1.05 / u.h]
B >> 2 * B [1 / u.h]

# Initial counts
A(1 / u.ml), B(1 / u.ml)

# First phase
S1 = Simulation(A | B)
S1.duration = 3 * u.h
S1.volume = 1 * u.ml

A + B >> Zero [0.1 / u.h]

S2 = Simulation(A | B)
S2.duration = (A <= 0) | (B <= 0)
S2.method = 'stochastic'
S2.volume = 1*u.ml

S = S1 + S2
S.unit_x, S.unit_y = u.h, 1 / u.ml
S.add_plot_params(A={'ylabel': "A (mL⁻¹)"}, B={'ylabel': "B (mL⁻¹)"},
                  vertical_lines=[3], tight_layout=True,
                  xlabel_fontsize=14, ylabel_fontsize=14,
                  suptitle='Mutual Annihilation Protocol', suptitle_fontsize=18)
S.plot_config.save_to = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + \
                        '/images/Mutual_Annihilation/Mutual_Annihilation.pdf'
S.repetitions = 10
S.run()

