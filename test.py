from mobspy import *
import os
from pint import UnitRegistry,Quantity

if __name__ == '__main__':

    # Check rate dimension
    A, B, C = BaseSpecies(3)
    Rev[A(100) + B(100) >> Zero] [0.1*u.litre/u.minute, 0.1/u.minute]
    MySim = Simulation(A | B)
    MySim.simulation_method = 'stochastic'
    MySim.volume = 1
    MySim.save_data = False
    MySim.compile()