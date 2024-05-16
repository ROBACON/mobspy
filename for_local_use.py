from mobspy import *
import numpy as np

if __name__ == '__main__':

    # TODO Write texts for Rev operator
    def test_rev():

        A, B, C = BaseSpecies()

        Rev[A + 4*B >> C][1, 2]
        Rev[A + 4 * B >> C][lambda r1, r2: (100-r1)*(100-r2), lambda r: r**3]

        S = Simulation(A | B | C)
        print(S.compile())
    test_rev()

