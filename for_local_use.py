from mobspy import *
from mobspy.modules.ode_operator import dt
from testutils import compare_model, compare_model_ignore_order


if __name__ == "__main__":


    A, B, Hey = BaseSpecies()
    D = New(A)

    A >> 2 * A[lambda r: 1 / u.hour * (1 + 10 / r)]
    (
            A + B
            >> Zero [
                lambda r1, r2: (1 * u.millimolar / u.hour)
                               * (1 + 10 * u.millimolar / r1 + 20 * u.millimolar / r2)
            ]
    )
    Hey >> Zero[lambda r: 1 / u.hour * (20 * r + 30 * r + 40 * r)]
    D >> 2 * D[lambda r: 20 / u.hour * r]

    S = Simulation(A | B | Hey | D)
    S.level = -1
    assert compare_model(S.compile(), "test_tools/model_33.txt")
