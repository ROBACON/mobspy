"""Tests for numpy integration with MobsPy."""

from __future__ import annotations

import numpy as np

import mobspy
from mobspy import BaseSpecies, Simulation, set_counts, u


def test_numpy_init_params():
    A, B = BaseSpecies()
    params = np.array([10.0, 1.0])
    A(params[0])
    A >> B @ params[1]
    MySim = Simulation(A | B)
    MySim.level = -1
    MySim.compile(verbose=False)
    cm = MySim._concrete_model
    assert cm is not None
    assert len(cm.species) == 2
    assert cm.species["A"] == 10.0
    assert cm.species["B"] == 0
    assert len(cm.reactions) == 1
    rxn = next(iter(cm.reactions.values()))
    assert "1.0" in rxn.kinetics


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
    A >> mobspy.Zero @ (lambda r: test_numpy_in_expression(r, 1))
    B >> mobspy.Zero @ (lambda r: test_numpy_in_expression(r, 2))
    C >> mobspy.Zero @ (lambda r: test_numpy_in_expression(r, 3))
    D >> mobspy.Zero @ (lambda r: test_numpy_in_expression(r, 4))
    A(100)
    S = Simulation(A | B | C | D)
    S.level = -1
    S.compile(verbose=False)
    cm = S._concrete_model
    assert cm is not None
    assert len(cm.species) == 4
    assert cm.species["A"] == 100
    assert len(cm.reactions) == 4
    kinetics = sorted(r.kinetics for r in cm.reactions.values())
    assert any("3+" in k or "(3+A)" in k or "+A)" in k for k in kinetics)
    assert any("3*" in k or "(3*" in k for k in kinetics)


def test_numpy_with_units():
    np_array = np.array([3])

    def test_numpy_in_expression(r):
        for a in np_array:
            b = a / u.hour
        return b

    A, B, C, D = BaseSpecies()
    A >> mobspy.Zero @ test_numpy_in_expression
    for a in np_array:
        B >> mobspy.Zero @ (a / u.hour)
    S = Simulation(A | B | C | D)
    S.level = -1
    S.compile(verbose=False)
    cm = S._concrete_model
    assert cm is not None
    assert len(cm.species) == 4
    assert len(cm.reactions) == 2
    # Both reactions should have numerical rate constants (3/hour -> ~0.000833/s)
    for rxn in cm.reactions.values():
        assert "0.000" in rxn.kinetics


def test_numpy_in_rates():
    np_array = np.array([1])
    for a in np_array:
        b = a
    A = BaseSpecies()
    A >> mobspy.Zero @ b
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
    A >> mobspy.Zero @ b
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
    A >> mobspy.Zero @ b
    model = set_counts({A: b})
    S = Simulation(model)
    S.level = -1
    S.plot_data = False
    S.step_size = 30
    S.run(plot_data=False)
    assert S.fres[A][-1] <= 10
