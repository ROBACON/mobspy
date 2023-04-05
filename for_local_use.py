from mobspy import *
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':

    Cell = BaseSpecies()
    Cell >> 2*Cell [1]
    A, B, C = New(Cell)

    S = Simulation(A | B | C)
    print(S.compile())

    def hi():
        Cell = BaseSpecies()
        Cell >> 2 * Cell[1]
        A, B, C = New(Cell)

        S = Simulation(A | B | C)
        print(S.compile())
    hi()

    def hi_inside_hi():
        hi()
    hi_inside_hi()






