from mobspy import *
import matplotlib.pyplot as plt
import os
import basico
import pint


if __name__ == '__main__':

    _S0 = BaseSpecies()

    _S0 >> Zero [1]

    S = Simulation(_S0)
    print(S.compile())


