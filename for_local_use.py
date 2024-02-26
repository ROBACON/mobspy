from mobspy import *
import numpy as np

if __name__ == '__main__':

    Location = BaseSpecies()
    Location.l1, Location.l2
    A, B, C, D = New(Location)

    # Before
    A.l1 + B.l1 >> 2 * C.l1
    C.l1 >> D.l1

    # After
    with Location.l1:
        A + B >> 2 * C
        C >> D



