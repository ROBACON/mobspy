from mobspy import *
import os

if __name__ == '__main__':

    # Check rate dimension
    def local():
        A, B, C = BaseSpecies(3)
        A(100*u.nanomolar*u.liter) + B(100) >> Zero [1]
        MySim = Simulation(A | B)
        MySim.simulation_method = 'stochastic'
        MySim.volume = 1*u.femtoliter
        MySim.save_data = False
        MySim.compile()

    local()
