from mobspy import *

if __name__ == '__main__':

    A = BaseSpecies()

    A >> Zero [1]

    S = Simulation(A)
    S.duration = 10*u.m**2
    S.run()
