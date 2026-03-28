"""Tests for parametric sweeps with ModelParameters."""

from __future__ import annotations

import mobspy
from mobspy import BaseSpecies, ModelParameters, Simulation, set_counts, simlog, u

simlog.global_simlog_level = -1


class TestParameterInReactions:
    def test_single_parameter_in_rate(self):
        A = BaseSpecies()
        k = ModelParameters(0.5)
        A >> mobspy.Zero[k]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert "k" in result
        assert result is not None

    def test_parameter_expression_in_rate(self):
        A = BaseSpecies()
        k = ModelParameters(2.0)
        A >> mobspy.Zero[3 * k]
        A(50)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert result is not None

    def test_two_parameters_in_rate(self):
        A = BaseSpecies()
        k1, k2 = ModelParameters(1.0, 0.5)
        A >> mobspy.Zero[k1 * k2]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert result is not None

    def test_parameter_as_initial_count(self):
        A = BaseSpecies()
        n0 = ModelParameters(50)
        A(n0)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert result is not None


class TestParameterSweeps:
    def test_single_parameter_sweep(self):
        A = BaseSpecies()
        k = ModelParameters([1, 2, 3])
        A >> mobspy.Zero[k]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert result is not None

    def test_parameter_sweep_with_expression(self):
        A = BaseSpecies()
        k = ModelParameters([0.5, 1.0, 1.5])
        A >> mobspy.Zero[2 * k]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert result is not None

    def test_parameter_sweep_as_count(self):
        A = BaseSpecies()
        n = ModelParameters([10, 50, 100])
        A >> mobspy.Zero[1]
        set_counts({"A": n})
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert result is not None

    def test_two_sweep_parameters(self):
        A = BaseSpecies()
        k1, k2 = ModelParameters([1, 2], [0.1, 0.2, 0.3])
        A >> mobspy.Zero[k1 + k2]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert result is not None


class TestParameterRenaming:
    def test_rename_parameter(self):
        A = BaseSpecies()
        k = ModelParameters(1.0)
        k.rename("rate_constant")
        A >> mobspy.Zero[k]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert "rate_constant" in result

    def test_rename_sweep_parameter(self):
        A = BaseSpecies()
        k = ModelParameters([1, 2, 3])
        k.rename("decay_rate")
        A >> mobspy.Zero[k]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert "decay_rate" in result

    def test_get_name_after_rename(self):
        k = ModelParameters(5.0)
        assert k.get_name() == "k"
        k.rename("my_param")
        assert k.get_name() == "my_param"


class TestParameterWithUnits:
    def test_parameter_with_rate_unit(self):
        A = BaseSpecies()
        k = ModelParameters(1 / u.hour)
        A >> mobspy.Zero[k]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert result is not None

    def test_parameter_sweep_with_units(self):
        A = BaseSpecies()
        k = ModelParameters([1 / u.hour, 2 / u.hour, 3 / u.hour])
        A >> mobspy.Zero[k]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert result is not None

    def test_parameter_has_units_flag(self):
        k_unit = ModelParameters(1 / u.hour)
        k_plain = ModelParameters(1.0)
        assert k_unit.has_units() is True
        assert k_plain.has_units() is False

    def test_convert_to_original_unit(self):
        k = ModelParameters(2 * u.mol / u.l)
        k.convert_to_original_unit()
        assert k.value.magnitude == (2 * u.mol / u.l).magnitude
        assert k.value.units == (2 * u.mol / u.l).units

    def test_two_parameters_with_different_units(self):
        A = BaseSpecies()
        k1, k2 = ModelParameters(1, [1 / u.hour, 2 / u.hour])
        A >> mobspy.Zero[k1 * k2]
        A(100)
        S = Simulation(A)
        S.level = -1
        result = S.compile()
        assert result is not None
