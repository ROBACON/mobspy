"""Tests for unit handling bugs in rate functions.

Covers:
- Bug 1: OverrideQuantity dual unit tracking (substance normalization in execute_quantity_op)
- Bug 2: QuantityConverter mol_power parsing (dimensionality dict instead of string index)
- Bug 4: _has_units boolean sentinel (True instead of "T")
"""

from __future__ import annotations

import pytest
from scipy.constants import N_A

from mobspy.modules.mobspy_expressions import (
    MobsPyExpression,
    OverrideQuantity,
    QuantityConverter,
    u,
)

ur = u.unit_registry_object


class TestQuantityConverterMolPower:
    """Bug 2: mol_power extraction uses dimensionality dict."""

    def test_mol_per_second(self) -> None:
        q = OverrideQuantity(1.0 * ur.moles / ur.seconds)
        result = QuantityConverter.convert_received_unit(q)
        assert abs(result.magnitude - N_A) / N_A < 1e-6
        assert "second" in str(result.units)

    def test_liters_per_mol_second(self) -> None:
        q = OverrideQuantity(1.0 * ur.liters / (ur.moles * ur.seconds))
        result = QuantityConverter.convert_received_unit(q)
        expected = 1.0 / N_A
        assert abs(result.magnitude - expected) / expected < 1e-6

    def test_mol_squared_per_second(self) -> None:
        q = OverrideQuantity(1.0 * ur.moles**2 / ur.seconds)
        result = QuantityConverter.convert_received_unit(q)
        expected = N_A**2
        assert abs(result.magnitude - expected) / expected < 1e-6

    def test_liters_squared_per_mol_squared_second(self) -> None:
        q = OverrideQuantity(1.0 * ur.liters**2 / (ur.moles**2 * ur.seconds))
        result = QuantityConverter.convert_received_unit(q)
        expected = 1.0 / N_A**2
        assert abs(result.magnitude - expected) / expected < 1e-6

    def test_no_substance(self) -> None:
        q = OverrideQuantity(1.0 / ur.seconds)
        result = QuantityConverter.convert_received_unit(q)
        assert abs(result.magnitude - 1.0) < 1e-10

    def test_concentration_mol_per_liter(self) -> None:
        q = OverrideQuantity(100.0 * ur.moles / ur.liters)
        result = QuantityConverter.convert_received_unit(q)
        expected = 100.0 * N_A
        assert abs(result.magnitude - expected) / expected < 1e-6

    def test_time_conversion(self) -> None:
        q = OverrideQuantity(60.0 / ur.minutes)
        result = QuantityConverter.convert_received_unit(q)
        assert abs(result.magnitude - 1.0) < 1e-10

    def test_per_hour(self) -> None:
        q = OverrideQuantity(3600.0 / ur.hour)
        result = QuantityConverter.convert_received_unit(q)
        assert abs(result.magnitude - 1.0) < 1e-10


class TestHasUnitsBoolean:
    """Bug 4: _has_units uses proper boolean."""

    def test_override_quantity_has_units(self) -> None:
        q = OverrideQuantity(1.0 * ur.seconds)
        assert q._has_units is True
        assert isinstance(q._has_units, bool)

    def test_expression_definer_default(self) -> None:
        expr = MobsPyExpression("A", None)
        assert expr._has_units is False
        assert isinstance(expr._has_units, bool)

    def test_has_units_propagates_in_operations(self) -> None:
        q = OverrideQuantity(5.0 / ur.seconds)
        expr = MobsPyExpression(
            "A",
            None,
            count_in_model=True,
            concentration_in_model=False,
            count_in_expression=False,
            concentration_in_expression=False,
        )
        # Trigger expression mode for the operation
        from mobspy.modules.mobspy_expressions import _ms_active_ctx

        token = _ms_active_ctx.set(True)
        try:
            result = q * expr
            assert result._has_units is True
            assert isinstance(result._has_units, bool)
        finally:
            _ms_active_ctx.reset(token)


class TestSubstanceNormalization:
    """Bug 1: execute_quantity_op normalizes [substance] on retry."""

    def test_mol_per_liter_plus_dimensionless_conc_path(self) -> None:
        """mol/L + 1/L should succeed on the concentration path."""
        from mobspy.modules.mobspy_expressions import _ms_active_ctx

        K_m = OverrideQuantity(100.0 * ur.moles / ur.liters)
        expr = MobsPyExpression(
            "A",
            None,
            count_in_model=True,
            concentration_in_model=False,
            count_in_expression=False,
            concentration_in_expression=False,
        )

        token = _ms_active_ctx.set(True)
        try:
            # This should not raise - the conc path should work after
            # normalizing [substance] -> counts
            result = K_m + expr
            assert isinstance(result, MobsPyExpression)
            # The conc_op should not be an Exception
            assert not isinstance(result._unit_conc_op, Exception), (
                f"conc_op should succeed but got: {result._unit_conc_op}"
            )
        finally:
            _ms_active_ctx.reset(token)

    def test_michaelis_menten_compiles(self) -> None:
        """Full Michaelis-Menten with concentration K_m should compile."""
        from mobspy.modules.mobspy_expressions import _ms_active_ctx

        k = OverrideQuantity(1.0 / ur.seconds)
        K_m = OverrideQuantity(100.0 * ur.moles / ur.liters)
        expr = MobsPyExpression(
            "A",
            None,
            count_in_model=True,
            concentration_in_model=False,
            count_in_expression=False,
            concentration_in_expression=False,
        )

        token = _ms_active_ctx.set(True)
        try:
            # k * r1 / (K_m + r1)
            numerator = k * expr
            denominator = K_m + expr
            result = numerator / denominator
            assert isinstance(result, MobsPyExpression)
        finally:
            _ms_active_ctx.reset(token)

    def test_hill_function_compiles(self) -> None:
        """Hill function with concentration K should compile."""
        from mobspy.modules.mobspy_expressions import _ms_active_ctx

        k = OverrideQuantity(1.0 / ur.seconds)
        K = OverrideQuantity(50.0 * ur.moles / ur.liters)
        n = 2
        expr = MobsPyExpression(
            "A",
            None,
            count_in_model=True,
            concentration_in_model=False,
            count_in_expression=False,
            concentration_in_expression=False,
        )

        token = _ms_active_ctx.set(True)
        try:
            result = k * expr**n / (K**n + expr**n)
            assert isinstance(result, MobsPyExpression)
        finally:
            _ms_active_ctx.reset(token)

    def test_incompatible_units_still_fail(self) -> None:
        """Adding meters to seconds should still fail both paths."""
        from mobspy.modules.mobspy_expressions import _ms_active_ctx

        q_meters = OverrideQuantity(1.0 * ur.meters)
        q_seconds = OverrideQuantity(1.0 * ur.seconds)

        token = _ms_active_ctx.set(True)
        try:
            with pytest.raises(Exception):  # noqa: B017
                _ = q_meters + q_seconds
        finally:
            _ms_active_ctx.reset(token)
