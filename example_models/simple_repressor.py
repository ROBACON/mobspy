import sys, os
import matplotlib.pyplot as plt
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

# Model definition
Chemical, Promoter, Protein = BaseSpecies(3)

Rev[Chemical + Promoter.inactive >> Promoter.active][2, 1]
rate = lambda pa, pi: f'{pi}/({pa} + {pi})'
Promoter.active + Promoter.inactive >> Promoter.active + Promoter.inactive + Protein [rate]
Protein >> Zero [0.1]
Promoter(100), Chemical(100)

# Simulation
MySim = Simulation(Chemical | Promoter | Protein)
MySim.save_data = False
MySim.plot_data = False
MySim.duration = 100
MySim.run()
MySim.plot_deterministic(Protein)