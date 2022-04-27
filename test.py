from mobspy import *

# TODO Plot has random order for species names

if __name__ == '__main__':

    A, B, AB, R = BaseSpecies(4)

    A(100) + B(10) >> AB + A[0.1]
    AB + R(1e4) >> 2*A[1]
    AB + R >> 2*B[1]

    MySim = Simulation(A | B | AB | R)
    MySim.save_data = False
    MySim.plot_data = True
    MySim.duration = 100
    MySim.simulation_method = 'stochastic'
    MySim.repetitions = 10
    MySim.plot.color = 'red'
    MySim.run()


