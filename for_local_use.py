from mobspy import *
import matplotlib.pyplot as plt


if __name__ == '__main__':
    A, B = BaseSpecies(2)

    A.a1, A.a2, B.b1, B.b2
    Combination = A*B

    Combination >> Zero [lambda r1: 0 if r1.b2 else 1]
    S = Simulation(Combination)
    print(S.compile())
