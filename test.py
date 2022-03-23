from mobspy import *
import os

if __name__ == '__main__':

    Ecoli = BaseSpecies(1)
    Ecoli.correct >> Ecoli.incorrect[1 / u.day]
    Ecoli.correct(1e3)

    MySim = Simulation(Ecoli)
    MySim.save_data = False
    MySim.plot_data = False
    MySim.compile()
