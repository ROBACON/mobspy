from mobspy import *
from mobspy.modules.ode_operator import dt

if __name__ == "__main__":

    A, B = BaseSpecies()

    dt[A] >> -0.1*A

    A(100)
    S = Simulation(A)
    print(S.generate_sbml()[0])