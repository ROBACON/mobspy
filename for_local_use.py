from mobspy import *

if __name__ == '__main__':

    A = BaseSpecies()

    A >> Zero [1]

    A(100)
    S = Simulation(A)
    S.duration = 10
    S.run()
