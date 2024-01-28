from mobspy import *

if __name__ == '__main__':

    A = BaseSpecies()

    A >> Zero [1]

    A(200)
    S = Simulation(A)
    S.duration = 100
    S.run()





