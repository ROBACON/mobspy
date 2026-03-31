"""Tests for the unit_handler module."""

from __future__ import annotations

import warnings

import pytest

from mobspy.exceptions import UnitError
from mobspy.modules.mobspy_expressions import u
from mobspy.modules.unit_handler import (
    check_dimension,
    convert_counts,
    convert_rate,
    convert_time,
    convert_volume,
    deep_copy_quantities,
    extract_length_dimension,
    time_convert_to_other_unit,
)


class TestDeepCopyQuantities:
    def test_copy_quantity(self):
        q = 5.0 * u.mL
        copied = deep_copy_quantities(q)
        assert copied.magnitude == q.magnitude
        assert str(copied.units) == str(q.units)

    def test_passthrough_non_quantity(self):
        assert deep_copy_quantities(42) == 42
        assert deep_copy_quantities("hello") == "hello"


class TestCheckDimension:
    def test_first_dimension_stored(self):
        assert check_dimension(None, 3) == 3

    def test_consistent_dimension(self):
        assert check_dimension(3, 3) == 3

    def test_inconsistent_dimension(self):
        with pytest.raises(UnitError, match="dimensions are not consistent"):
            check_dimension(3, 2)

    def test_inconsistent_with_context(self):
        with pytest.raises(UnitError, match="some context"):
            check_dimension(3, 2, error_context="some context")


class TestExtractLengthDimension:
    def test_no_length(self):
        result = extract_length_dimension("[time] ** -1", None)
        assert result is False

    def test_with_length_no_order(self):
        result = extract_length_dimension("[length] ** 3 / [time]", None)
        assert result == 3

    def test_with_order_2(self):
        result = extract_length_dimension(
            "[length] ** 3 / [time]", None, reaction_order=2
        )
        assert result == 3

    def test_order_1_no_length(self):
        result = extract_length_dimension("[time] ** -1", None, reaction_order=1)
        assert result == 0

    def test_order_1_with_length_raises(self):
        with pytest.raises(UnitError, match="Unimolecular reaction"):
            extract_length_dimension("[length] ** 3 / [time]", None, reaction_order=1)


class TestConvertRate:
    def test_plain_number_passthrough(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result, dimension, _ = convert_rate(5.0, 1, None)
            assert result == 5.0
            assert dimension is None

    def test_time_only_rate(self):
        rate = 1.0 / u.s
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result, dim, converted = convert_rate(rate, 1, None)
        assert converted is True
        assert isinstance(result, float)
        assert abs(result - 1.0) < 1e-10

    def test_time_rate_minutes(self):
        rate = 60.0 / u.min
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result, _, converted = convert_rate(rate, 1, None)
        assert converted is True
        assert abs(result - 1.0) < 1e-10


class TestConvertCounts:
    def test_plain_number_passthrough(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = convert_counts(100, 1.0, 3)
        assert result == 100

    def test_dimensionless_quantity(self):
        q = 50.0 * u.dimensionless
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = convert_counts(q, 1.0, 3)
        assert result == 50.0

    def test_invalid_unit_raises(self):
        q = 5.0 * u.s
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            with pytest.raises(UnitError, match="neither a count nor a concentration"):
                convert_counts(q, 1.0, 3)


class TestConvertVolume:
    def test_plain_number_passthrough(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = convert_volume(5.0, 3)
        assert result == 5.0

    def test_quantity_conversion(self):
        vol = 1.0 * u.L
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = convert_volume(vol, 3)
        assert isinstance(result, float)
        assert result == pytest.approx(1.0, rel=1e-6)


class TestConvertTime:
    def test_plain_number_passthrough(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = convert_time(30.0)
        assert result == 30.0

    def test_seconds_quantity(self):
        t = 60.0 * u.s
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = convert_time(t)
        assert result == pytest.approx(60.0)

    def test_minutes_to_seconds(self):
        t = 2.0 * u.min
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = convert_time(t)
        assert result == pytest.approx(120.0)

    def test_non_time_quantity_returns_none(self):
        q = 5.0 * u.meter
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = convert_time(q)
        assert result is None


class TestTimeConvertToOtherUnit:
    def test_convert_seconds_to_minutes(self):
        t = 120.0 * u.s
        result = time_convert_to_other_unit(t, "minute")
        assert result == pytest.approx(2.0)

    def test_plain_number_passthrough(self):
        result = time_convert_to_other_unit(30, "minute")
        assert result == 30

    def test_non_time_returns_none(self):
        q = 5.0 * u.meter
        result = time_convert_to_other_unit(q, "second")
        assert result is None
