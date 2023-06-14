from mobspy import *
import basico
import matplotlib.pyplot as plt
import os
import inspect
from mobspy.sbml_simulator.builder import build

if __name__ == '__main__':

    Age, Color, Mortal = BaseSpecies()

    Age.young >> Age.old [1]
    Color.red >> Color.blue [1]
    Color.blue >> Color.yellow [1]
    Mortal >> Zero [1]

    Cell = Age*Color*Mortal
    C1, C2  = New(Cell)
    S = Simulation(C1 | C2)
    S.compile()









