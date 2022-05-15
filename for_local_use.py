from mobspy import *

if __name__ == '__main__':

    A, B = BaseSpecies(2)
    A >> B [1]
    MySim = Simulation(A | B)
    MySim.save_data = False
    MySim.plot_data = False
    MySim.run()
    pass
