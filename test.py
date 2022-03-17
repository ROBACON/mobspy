from mobspy import *
import os
from pint import UnitRegistry,Quantity

if __name__ == '__main__':

    # Check rate dimension
    def local():
        A, B, C = BaseSpecies(3)
        A(100) + B(100) >> Zero [1.8e-3*u.liter/(u.nanomoles*u.second)]
        MySim = Simulation(A | B)
        MySim.simulation_method = 'stochastic'
        MySim.volume = 1
        MySim.save_data = False
        MySim.compile()

    local()
