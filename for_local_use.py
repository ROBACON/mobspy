from mobspy import *
import numpy as np

if __name__ == '__main__':

    A = BaseSpecies()
    B, C = New(A)

    A.a1, A.a2, A.a3

    B + C >> Zero [lambda r1, r2: 100 if A(r1) == A(r2) else 0]

    S = Simulation(B | C)
    print(S.compile())



