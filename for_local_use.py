from mobspy import *

if __name__ == "__main__":

    A, B = BaseSpecies()
    A.a1, A.a2, A.a3, B.b1, B.b2
    C = A*B

    ~C.a1 >> Zero [10]

    S = Simulation(A | C)
    print(S.compile())