from mobspy import *

if __name__ == '__main__':
    # Assignment_Operator._asg_context = True

    A, B, C = BaseSpecies()

    A >> Zero [1]
    B.assign(2*A)

    A(10)
    S = Simulation(A | B)
    S.duration = 10
    S.run()
