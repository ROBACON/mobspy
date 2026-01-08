from mobspy import *
from mobspy.modules.ode_operator import dt
from testutils import compare_model, compare_model_ignore_order


if __name__ == "__main__":

    def ode_compartments():
        """ODE combined with regular CRN reactions"""
        A = BaseSpecies()
        A.c1, A.c2

        Zero >> A.c1[1]
        A.c1 >> A.c2[1]

        # ODE for A
        dt[A] >> -0.1 * A

        S = Simulation(A)
        assert compare_model(S.compile(), "test_tools/model_ode_compartments.txt")
    ode_compartments()
