import sys, os
import matplotlib.pyplot as plt
from mobspy import *

# Model definition
"""
    This model is a simple single repressor
"""

Chemical, Promoter, Protein = BaseSpecies()

Rev[Chemical + Promoter.inactive >> Promoter.active][2, 1]
Promoter >> Promoter + Protein [lambda promoter: 1 if promoter.inactive else 0]
Protein >> Zero [0.1]
Promoter.inactive(100), Chemical(1000)

# Simulation
MySim = Simulation(Chemical | Promoter | Protein)
MySim.save_data = False
MySim.plot_data = False
MySim.duration = 200
MySim.run()
MySim.plot_deterministic(Protein)