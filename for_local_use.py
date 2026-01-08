from mobspy import *
from mobspy.modules.ode_operator import dt
from testutils import compare_model, compare_model_ignore_order


if __name__ == "__main__":

    A, B = BaseSpecies()
    A.a1, A.a2
    B.b1, B.b2

    C = A * B

    All[C](100)
    C >> All[C] [1]

    S = Simulation(C)
    S.level = -1
    assert compare_model(S.compile(), "test_tools/model_21.txt")