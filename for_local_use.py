from mobspy import *
import matplotlib.pyplot as plt
import os
import basico
import pint
from numpy import issubdtype, floating, integer, array


if __name__ == '__main__':

    idk = array([1, 2, 3])

    for a in idk:
        print(issubdtype(a, floating))
        print(issubdtype(a, integer))



