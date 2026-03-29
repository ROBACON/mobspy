"""Tests for exponent handling in rate functions."""

from __future__ import annotations

from mobspy import BaseSpecies, Simulation, u


def test_exponent_basic():
    default_rate = 3e-11 / u.min
    k_coeff = 1208585 * u.counts / u.mL
    n_hill = 1.0
    init_a_qty = 1000 * u.counts / u.mL

    def antibio_uptake_senders(a):
        mu_monod = (default_rate) * (
            k_coeff**n_hill / ((k_coeff**n_hill) + init_a_qty**n_hill)
        )
        return mu_monod * a

    A = BaseSpecies()
    A(init_a_qty)
    A >> 2 * A[antibio_uptake_senders]
    S = Simulation(A)
    S.duration = 50 * u.min
    S.unit_x = u.min
    S.unit_y = u.counts / u.mL
    S.plot_data = False
    S.save_data = False
    S.run()


def test_rate_function_with_exp():
    default_rate_ = 3e-11 / u.min
    qty_A = 1000 * u.counts / u.mL
    K = 1e10 * u.counts / u.mL
    n = 1.0

    def antibio_uptake_senders(a):
        rate_cst = default_rate_ * (K**n / (K**n + qty_A**n))
        return rate_cst * a

    A = BaseSpecies()
    A(qty_A)
    A >> 2 * A[antibio_uptake_senders]
    S = Simulation(A)
    S.run(duration=50 * u.min, unit_x=u.min, unit_y=u.counts / u.mL, plot_data=False)
