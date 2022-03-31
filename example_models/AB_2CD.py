import sys, os
from mobspy import *

A, B, C, D = BaseSpecies(4)
A(200) + B(100) >> 2*C + D [0.1]
My_Sim = Simulation(A | B | C | D)
My_Sim.save_data = False
My_Sim.duration = 5
My_Sim.volume = 10
My_Sim.simulation_method = 'stochastic'
My_Sim.run()
