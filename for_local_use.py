from mobspy import *

if __name__ == '__main__':

    C = BaseSpecies(1)
    C.b1(100), C.b2(100)
    C.b1 >> Zero [0.1]
    MySim = Simulation(C)
    MySim.simulation_method = 'stochastic'
    MySim.save_data = False
    MySim.plot_data = False
    MySim.run()
    MySim.plot_stochastic(C.b2)
    pass
