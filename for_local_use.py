from mobspy import *


if __name__ == '__main__':
    A, B = BaseSpecies()

    A.a1.assign(B * 5)

    S = Simulation(A | B)
    print(S.compile())

    exit()
    A, B = BaseSpecies()

    A.assign(B*5)

    S = Simulation(A)
    print(S.compile())
    exit()


    A, B, C = BaseSpecies()

    A.assign(2*B*(u.l/u.min))

    B(100)
    S = Simulation(A | B | C)
    # S.plot_data = False
    S.duration = 10
    S.step_size = 5
    S.run()

    exit()

    A, B = BaseSpecies()

    A.assign(5 * B * (u.l / u.s) + 10 * B * (1 / u.s))

    exit()


