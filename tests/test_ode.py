"""Tests for ODE syntax."""

from __future__ import annotations

from mobspy import BaseSpecies, New, Simulation, Zero
from mobspy.modules.functions import ms_exp
from mobspy.modules.ode_operator import dt

from .conftest import compare_model


def test_ode_syntax_basic():
    A = BaseSpecies()
    dt[A] += -0.1 * A
    A(100)
    S = Simulation(A)
    assert compare_model(S.compile(), "model_ode_syntax_basic.txt")


def test_ode_syntax_two_species():
    A, B = BaseSpecies()
    dt[A] += -0.1 * A + 0.05 * B
    A(100), B(50)
    S = Simulation(A | B)
    assert compare_model(S.compile(), "model_ode_syntax_two_species.txt")


def test_ode_applied_to_species():
    A = BaseSpecies()
    B = BaseSpecies()
    B.b1
    dt[A] += A
    dt[B.b1] += B.b1
    S = Simulation(A | B)
    assert compare_model(S.compile(), "model_ode_applied_to_species.txt")


def test_ode_neg_test():
    Neg, NegR = BaseSpecies()
    NegR.comp1
    dt[Neg] += -Neg
    dt[NegR] += -NegR.comp1
    S = Simulation(Neg | NegR)
    assert compare_model(S.compile(), "model_ode_neg_test.txt")


def test_ode_compartments():
    A = BaseSpecies()
    A.c1, A.c2
    Zero >> A.c1[1]
    A.c1 >> A.c2[1]
    dt[A] += -0.1 * A
    S = Simulation(A)
    assert compare_model(S.compile(), "model_ode_compartments.txt")


def test_ode_complex_expressions():
    A, B, C, D = BaseSpecies()
    dt[A] += 100 / (1 + B**2) - 0.1 * A
    dt[B] += (A * C) / (10 + A + C) - B / (5 + B)
    dt[C] += (A / (1 + A)) * (B / (1 + B)) - 0.05 * C * D
    dt[D] += (A**2 + B**2) / (100 + A**2 + B**2) * (1 - D / 1000)
    A(10), B(10), C(10), D(10)
    S = Simulation(A | B | C | D)
    assert compare_model(S.compile(), "model_ode_complex_expressions.txt")


def test_ode_inheritance():
    Mortal = BaseSpecies()
    Human, Animal = New(Mortal)
    dt[Mortal] += -0.1 * Mortal
    Human(100), Animal(50)
    S = Simulation(Human | Animal)
    assert compare_model(S.compile(), "model_ode_inheritance.txt")


def test_ode_with_functions():
    A = BaseSpecies()
    dt[A] += 1 / (1 + ms_exp(A / 1000))
    A(100)
    S = Simulation(A)
    assert compare_model(S.compile(), "model_ode_with_functions.txt")


def test_ode_neg():
    A = BaseSpecies()
    dt[A] -= A
    S = Simulation(A)
    assert compare_model(S.compile(), "model_ode_neg.txt")
