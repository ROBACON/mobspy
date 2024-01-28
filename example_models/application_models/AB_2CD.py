import sys, os
from mobspy import *

"""
    This is a basic MobsPy example
    Here we declare four species and the single reaction A + B >> 2*C + D
    It enunciates the basics of MobsPy, like declaring a species using the BaseSpecies constructor
    Species counts are assigned using the call method
    Rates are assigned using the getitem brackets []
"""
A, B, C, D = BaseSpecies()
A(200) + B(100) >> 2*C + D [0.1]
My_Sim = Simulation(A | B | C | D)
My_Sim.save_data = False
My_Sim.duration = 5
My_Sim.volume = 10
My_Sim.simulation_method = 'stochastic'
My_Sim.run()
