import sys, os
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

A, B, C = BaseSpecies(3)
A(300) + B(200) >> C [0.1]
MySim = Simulation(A | B | C)
MySim.plot_data = False
MySim.save_data = False
MySim.run()

