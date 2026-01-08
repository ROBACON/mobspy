from mobspy.modules.mobspy_parameters import Internal_Parameter_Constructor
from mobspy.modules.ode_operator import dt
from mobspy import *
from testutils import compare_model, compare_model_ignore_order

def test_ode_syntax_basic():
    A = BaseSpecies()

    dt[A] >> -0.1 * A

    A(100)
    S = Simulation(A)
    assert compare_model(S.compile(), "test_tools/model_ode_syntax_basic.txt")


def test_ode_syntax_two_species():
    """ODE with two species: dA/dt = -0.1*A + 0.05*B"""
    A, B = BaseSpecies()

    dt[A] >> -0.1 * A + 0.05 * B

    A(100), B(50)
    S = Simulation(A | B)
    assert compare_model(S.compile(), "test_tools/model_ode_syntax_two_species.txt")

def test_ode_applied_to_species():

    A = BaseSpecies()
    B = BaseSpecies()
    B.b1

    dt[A] >> A
    dt[B.b1] >> B.b1

    S = Simulation(A | B)
    assert compare_model(S.compile(), "test_tools/model_ode_applied_to_species.txt")


def test_ode_neg_test():

    Neg, NegR = BaseSpecies()
    NegR.comp1

    dt[Neg] >> -Neg
    dt[NegR] >> -NegR.comp1

    S = Simulation(Neg | NegR)
    assert compare_model(S.compile(), "test_tools/model_ode_neg_test.txt")



