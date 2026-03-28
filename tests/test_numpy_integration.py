"""Tests for numpy integration with MobsPy."""

from __future__ import annotations

import numpy as np

import mobspy
from mobspy import BaseSpecies, Simulation, set_counts, u

from .conftest import compare_model


def test_numpy_init_params():
    expected = """
Species
A,10.0
B,0

Mappings
A :
A
B :
B

Parameters
volume,1

Reactions
reaction_0,{'re': [(1, 'A')], 'pr': [(1, 'B')], 'kin': 'A * 1.0'}
"""
    A, B = BaseSpecies()
    params = np.array([10.0, 1.0])
    A(params[0])
    A >> B[params[1]]
    MySim = Simulation(A | B)
    MySim.level = -1
    result = MySim.compile()
    assert result == expected


def test_numpy_in_expression_function():
    def test_numpy_in_expression(r, op):
        np_array = np.array([3])

        for a in np_array:
            if op == 1:
                b = a + r
            elif op == 2:
                b = a - r
            elif op == 3:
                b = a * r
            elif op == 4:
                b = a / r
            else:
                b = 0
        return b

    A, B, C, D = BaseSpecies()
    A >> mobspy.Zero[lambda r: test_numpy_in_expression(r, 1)]
    B >> mobspy.Zero[lambda r: test_numpy_in_expression(r, 2)]
    C >> mobspy.Zero[lambda r: test_numpy_in_expression(r, 3)]
    D >> mobspy.Zero[lambda r: test_numpy_in_expression(r, 4)]
    A(100)
    S = Simulation(A | B | C | D)
    S.level = -1
    assert compare_model(S.compile(), "model_46.txt")


def test_numpy_with_units():
    np_array = np.array([3])

    def test_numpy_in_expression(r):

        for a in np_array:
            b = a / u.hour
        return b

    A, B, C, D = BaseSpecies()
    A >> mobspy.Zero[lambda r: test_numpy_in_expression(r)]
    for a in np_array:
        B >> mobspy.Zero[a / u.hour]
    S = Simulation(A | B | C | D)
    S.level = -1
    assert compare_model(S.compile(), "model_47.txt")


def test_numpy_in_rates():
    np_array = np.array([1])
    for a in np_array:
        b = a
    A = BaseSpecies()
    A >> mobspy.Zero[b]
    A(200)
    S = Simulation(A)
    S.level = -1
    S.plot_data = False
    S.step_size = 30
    S.run(plot_data=False)
    assert S.fres[A][-1] <= 10


def test_numpy_in_counts():
    np_array = np.array([200])
    for a in np_array:
        b = a
    A = BaseSpecies()
    A >> mobspy.Zero[b]
    A(b)
    S = Simulation(A)
    S.level = -1
    S.plot_data = False
    S.step_size = 30
    S.run(plot_data=False)
    assert S.fres[A][-1] <= 10


def test_numpy_in_set_counts():
    np_array = np.array([200])
    for a in np_array:
        b = a
    A = BaseSpecies()
    A >> mobspy.Zero[b]
    model = set_counts({A: b})
    S = Simulation(model)
    S.level = -1
    S.plot_data = False
    S.step_size = 30
    S.run(plot_data=False)
    assert S.fres[A][-1] <= 10
