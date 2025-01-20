from mobspy import *
from mobspy import modules

if __name__ == '__main__':

    k1 = ModelParameters(1)
    A, B, C, D = BaseSpecies()

    A >> 2 * B[modules.mobspy_parameters._Internal_Parameter_Constructor('k1', 3)]
    B >> C + D[modules.mobspy_parameters._Internal_Parameter_Constructor('k2', 1.4)]

    A(100)
    S = Simulation(A | B | C | D)
    S.duration = 10
    print(S.compile())
    # print(S.generate_sbml()[0])

    S.update_model()


