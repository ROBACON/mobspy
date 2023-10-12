from mobspy import *

if __name__ == '__main__':
    A = BaseSpecies()

    A.a1, A.a2, A.a3

    S = Simulation(A)
    S.level = -1
    with S.event_time(5):
        All[A](f'{A} + 1')

    with S.event_time(10):
        All[A.a1](f'{A} + 1')

    with S.event_time(15):
        A.a1(f'{A} + 1')

    print(S.compile())

