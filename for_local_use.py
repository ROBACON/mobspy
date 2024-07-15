from mobspy import *
import matplotlib.pyplot as plt
import os
import basico
import pint


if __name__ == '__main__':

    # Used for Model 56
    try:
        _S1 = BaseSpecies()
    except SystemExit:
        print('Nice')

    S0, S1, S2 = BaseSpecies()

    S0 >> Zero [1]
    S1 >> Zero [1]
    S2 >> Zero [1]

    S = Simulation(S0 | S1 | S2)
    print(S.compile())



