from mobspy import *
import os


if __name__ == '__main__':

    Yo, B, C =  BaseSpecies(3)
    Yo.burning, Yo.for_you
    A = New(Yo)
    A.to_b(20) >> B[5]
    A.to_c(20) >> C[5]
    MySim = Simulation(A | B | C).compile()
    exit()

    A, B, C, D = BaseSpecies(4)
    A(20) + B(10) >> 2*C + D [0.1]
    MySim = Simulation(A | B | C | D)
    MySim.duration = 100
    MySim.save_data = False
    MySim.plot_data = False
    MySim.run()
    MySim.plot_deterministic(A, C)
    exit()


