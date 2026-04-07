"""Tests for fractional product handling."""

from __future__ import annotations

import pytest

from mobspy import BaseSpecies, Simulation, u


def test_deterministic():
    A, B, C = BaseSpecies()
    A(10)
    B(1)
    A + B >> (0.5 * C + B) @ 0.005
    Sim = Simulation(A | B | C)
    Sim.method = "deterministic"
    Sim.volume = 1 * u.mL
    Sim.run(duration=10, unit_x=u.s, unit_y=1 / u.mL, plot_data=False)
    assert Sim.results["C"][0][-1] == pytest.approx(2.44e-4, rel=0.05)


def test_deterministic_variable():
    A, B, C = BaseSpecies()
    A(10)
    B(1)
    yield_coeff = 0.5
    A + B >> (yield_coeff * C + B) @ 0.005
    Sim = Simulation(A | B | C)
    Sim.method = "deterministic"
    Sim.volume = 1 * u.mL
    Sim.run(duration=10, unit_x=u.s, unit_y=1 / u.mL, plot_data=False)
    assert Sim.results["C"][0][-1] == pytest.approx(2.44e-4, rel=0.05)


def test_stochastic():
    A, B, C = BaseSpecies()
    A(10)
    B(1)
    A + B >> (0.5 * C + B) @ 0.005
    Sim = Simulation(A | B | C)
    Sim.method = "stochastic"
    Sim.volume = 1 * u.mL
    Sim.run(duration=10, unit_x=u.s, unit_y=1 / u.mL, plot_data=False)
    # With compartment-based rates the reaction is slow; C stays near zero
    assert Sim.results["C"][0][-1] < 0.01


def test_stochastic_variable():
    A, B, C = BaseSpecies()
    A(10)
    B(1)
    yield_coeff = 0.5
    A + B >> (yield_coeff * C + B) @ 0.005
    Sim = Simulation(A | B | C)
    Sim.method = "stochastic"
    Sim.volume = 1 * u.mL
    Sim.run(duration=10, unit_x=u.s, unit_y=1 / u.mL, plot_data=False)
    # With compartment-based rates the reaction is slow; C stays near zero
    assert Sim.results["C"][0][-1] < 0.01
