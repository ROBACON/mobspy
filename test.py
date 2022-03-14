from mobspy import *
import os
from pint import UnitRegistry,Quantity

if __name__ == '__main__':

    # Check rate dimension
    A, B, C = BaseSpecies(3)
    A(200) + B(100) >> C[1]
    MySim = Simulation(A | B | C)
    MySim.simulation_method = 'stochastic'
    MySim.volume = 1 * u.millilitre
    MySim.save_data = False
    MySim.run()