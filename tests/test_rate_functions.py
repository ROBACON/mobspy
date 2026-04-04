"""Tests for rate function variations (multiplication order, units)."""

from __future__ import annotations

from mobspy import BaseSpecies, Simulation, u


def test_rate_function_right_multiplication():
    default_rate = 3e-11 / u.min

    def antibio_uptake_senders(a):
        return default_rate * a

    A = BaseSpecies()
    A(1000 * u.counts / u.mL)
    A >> 2 * A @ antibio_uptake_senders
    S = Simulation(A)
    S.run(duration=1222 * u.min, unit_x=u.min, unit_y=u.counts / u.mL, plot_data=False)


def test_rate_function_inside_multiplication():
    def antibio_uptake_senders(a):
        return (3e-11 / u.min) * a

    A = BaseSpecies()
    A(1000 * u.counts / u.mL)
    A >> 2 * A @ antibio_uptake_senders
    S = Simulation(A)
    S.run(duration=1222 * u.min, unit_x=u.min, unit_y=u.counts / u.mL, plot_data=False)


def test_rate_function_left_multiplication():
    default_rate = 3e-11 / u.min

    def antibio_uptake_senders(a):
        return a * default_rate

    A = BaseSpecies()
    A(1000 * u.counts / u.mL)
    A >> 2 * A @ antibio_uptake_senders
    S = Simulation(A)
    S.run(duration=1222 * u.min, unit_x=u.min, unit_y=u.counts / u.mL, plot_data=False)
