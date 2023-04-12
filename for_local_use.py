from mobspy import *
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':

    A, C = BaseSpecies()
    A.a1, A.a2
    B = New(A)
    B.b1, B.b2

    model = set_counts({All['B.a1']: 100, C: 200*u.mols, 'A.a1': 100, A.a2: 50})
    S = Simulation(model)
    print(S.compile())






