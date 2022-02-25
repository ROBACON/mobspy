import sys, os
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

#TODO fix one characteristic error

A, S = BaseSpecies(2)
A.not_hi
A.hi(1000) >> A.hi + S[0.1]
Sim = Simulation(A | S)
Sim.simulation_method = 'stochastic'
Sim.repetitions = 1
Sim.duration = 10
Sim.save_data = False
Sim.run()