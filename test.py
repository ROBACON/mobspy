from mobspy import *
import os
from pint import UnitRegistry,Quantity

if __name__ == '__main__':

    # Check rate dimension
    A, B, C = BaseSpecies(3)
    A(100) >> Zero [0.1/u.minute]
    MySim = Simulation(A)
    MySim.simulation_method = 'stochastic'
    MySim.volume = 1
    MySim.save_data = False
    MySim.compile()