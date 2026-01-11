from mobspy.modules.mobspy_parameters import Internal_Parameter_Constructor
from mobspy.modules.ode_operator import dt
from mobspy import *
from testutils import compare_model, compare_model_ignore_order
from mobspy.modules.functions import ms_exp

def test_ode_syntax_basic():
    A = BaseSpecies()

    dt[A] += -0.1 * A

    A(100)
    S = Simulation(A)
    assert compare_model(S.compile(), "test_tools/model_ode_syntax_basic.txt")


def test_ode_syntax_two_species():
    """ODE with two species: dA/dt = -0.1*A + 0.05*B"""
    A, B = BaseSpecies()

    dt[A] += -0.1 * A + 0.05 * B

    A(100), B(50)
    S = Simulation(A | B)
    assert compare_model(S.compile(), "test_tools/model_ode_syntax_two_species.txt")

def test_ode_applied_to_species():

    A = BaseSpecies()
    B = BaseSpecies()
    B.b1

    dt[A] += A
    dt[B.b1] += B.b1

    S = Simulation(A | B)
    assert compare_model(S.compile(), "test_tools/model_ode_applied_to_species.txt")


def test_ode_neg_test():

    Neg, NegR = BaseSpecies()
    NegR.comp1

    dt[Neg] += -Neg
    dt[NegR] += -NegR.comp1

    S = Simulation(Neg | NegR)
    assert compare_model(S.compile(), "test_tools/model_ode_neg_test.txt")


def test_ode_compartments():

    """ODE combined with regular CRN reactions"""
    A = BaseSpecies()
    A.c1, A.c2

    Zero >> A.c1[1]
    A.c1 >> A.c2[1]

    # ODE for A
    dt[A] += -0.1 * A

    S = Simulation(A)
    assert compare_model(S.compile(), "test_tools/model_ode_compartments.txt")


def test_ode_complex_expressions():

    """ODE with various complex expressions: Hill functions, Michaelis-Menten, feedback loops"""
    A, B, C, D = BaseSpecies()

    # Hill-style repression: production inhibited by B
    dt[A] += 100 / (1 + B ** 2) - 0.1 * A

    # Michaelis-Menten with multiple substrates
    dt[B] += (A * C) / (10 + A + C) - B / (5 + B)

    # Nested fractions and mixed operations
    dt[C] += (A / (1 + A)) * (B / (1 + B)) - 0.05 * C * D

    # Complex feedback with powers and sums
    dt[D] += (A ** 2 + B ** 2) / (100 + A ** 2 + B ** 2) * (1 - D / 1000)

    A(10), B(10), C(10), D(10)
    S = Simulation(A | B | C | D)
    assert compare_model(S.compile(), "test_tools/model_ode_complex_expressions.txt")


def test_ode_inheritance():

    """ODE applied to parent affects all children"""
    Mortal = BaseSpecies()
    Human, Animal = New(Mortal)

    # Decay applied to Mortal affects both Human and Animal
    dt[Mortal] += -0.1 * Mortal

    Human(100), Animal(50)
    S = Simulation(Human | Animal)
    assert compare_model(S.compile(), "test_tools/model_ode_inheritance.txt")


def test_ode_with_functions():
    A = BaseSpecies()

    dt[A] += 1 / (1 + ms_exp(A / 1000))

    A(100)
    S = Simulation(A)
    assert compare_model(S.compile(), "test_tools/model_ode_with_functions.txt")


def test_ode_neg():

    A = BaseSpecies()

    dt[A] -= A

    S = Simulation(A)
    assert compare_model(S.compile(), "test_tools/model_ode_neg.txt")





