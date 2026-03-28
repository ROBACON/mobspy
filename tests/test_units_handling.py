"""Tests for unit handling with Pint."""

from __future__ import annotations

from pathlib import Path

from mobspy import BaseSpecies, Simulation, u


def test_def_by_str():
    import yaml

    config_path = Path(__file__).resolve().parent / "testconfig_units.yaml"
    with open(config_path, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    A0 = u(config["A0"])
    rate_constant = u(config["rate_constant"])

    assert str(A0) == "1.0 / milliliter"
    assert str(rate_constant) == "1.0 / minute"

    def rate_fn(a):
        assert str(A0) == "1.0 / milliliter"
        assert str(rate_constant) == "1.0 / minute"
        return rate_constant

    A = BaseSpecies()
    A(A0)
    A >> 2 * A[rate_fn]
    MySim = Simulation(A)
    MySim.volume = 1 * u.mL
    MySim.compile()


def test_def_by_quantity():
    A0 = 1 / u.mL
    rate_constant = 1 / u.min

    assert str(A0) == "1.0 / milliliter"
    assert str(rate_constant) == "1.0 / minute"

    def rate_fn(a):
        assert str(A0) == "1.0 / milliliter"
        assert str(rate_constant) == "1.0 / minute"
        return rate_constant

    A = BaseSpecies()
    A(A0)
    A >> 2 * A[rate_fn]
    MySim = Simulation(A)
    MySim.volume = 1 * u.mL
    MySim.compile()


def test_def_by_quantity_with_run():
    A0 = 1 / u.mL
    rate_constant = 1 / u.min
    duration = 5 * u.min

    def rate_fn(a):
        assert str(A0) == "1.0 / milliliter"
        assert str(rate_constant) == "1.0 / minute"
        return rate_constant

    A = BaseSpecies()
    A(A0)
    A >> 2 * A[rate_fn]
    MySim = Simulation(A)
    MySim.volume = 1 * u.mL
    MySim.compile()
    MySim.run(duration=duration, unit_x=u.min, unit_y=1 / u.mL, plot_data=False)


def test_print_vs_str():
    A0 = 1 / u.mL
    rate_constant = 1 / u.min
    duration = 5 * u.min

    def rate_fn(a):
        def fn_a0():
            return 1 / u.mL

        fn_a0()
        assert str(A0) == "1.0 / milliliter"
        assert str(rate_constant) == "1.0 / minute"
        return rate_constant

    A = BaseSpecies()
    A(A0)
    A >> 2 * A[rate_fn]
    MySim = Simulation(A)
    MySim.volume = 1 * u.mL
    MySim.compile()
    MySim.run(duration=duration, unit_x=u.min, unit_y=1 / u.mL, plot_data=False)
