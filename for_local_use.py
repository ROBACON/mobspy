import dateutil.rrule

from mobspy import *
import numpy as np

if __name__ == '__main__':

    A, B = BaseSpecies()
    A + B >> Zero [1]

    A(200), B(200)
    Sim = Simulation(A | B)
    Sim.level = -1
    descr = Sim.compile()
    Sim.plot_data = False
    Sim.duration = 30 * u.hour
    Sim.volume = 1*u.m**3
    Sim.run()

    assert Sim._parameters_for_sbml['volume'][0] > 100
    assert Sim.fres['Time'][-1] > 100





