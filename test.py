from mobspy import *
import os
from pint import UnitRegistry,Quantity

if __name__ == '__main__':

    # Check rate dimension
    A, B, C = BaseSpecies(3)
    A.s1, A.s2, A.s3
    A.s2 + B >> A.s3 + A [1]
    MySim = Simulation(A | B | C)
    MySim.simulation_method = 'stochastic'
    MySim.volume = 1
    MySim.save_data = False
    MySim.default_order = All
    MySim.compile()