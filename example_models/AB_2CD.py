import sys, os

# What does this part below do?
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

A, B, C, D = BaseSpecies(4)
A(200) + B(100) >> 2*C + D [0.1]
My_Sim = Simulation(A | B | C | D)
My_Sim.save_data = False
My_Sim.duration = 1
My_Sim.run()