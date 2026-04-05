"""Tests for parameter estimation helper functions.

Tests the pure validation/helper functions in
``mobspy.parameter_estimation_data_loader.parameter_estimation_scripts``
without requiring BasiCO/COPASI.
"""

from __future__ import annotations

from typing import Any
from unittest.mock import MagicMock

import pandas as pd
import pytest

from mobspy.exceptions import ParameterError
from mobspy.parameter_estimation_data_loader.parameter_estimation_scripts import (
    _extract_results,
    _resolve_experimental_data,
    _validate_bound,
    _validate_estimation_inputs,
)

# ------------------------------------------------------------------
# Helpers for building mock parameter objects
# ------------------------------------------------------------------


def _make_param(name: str, value: float) -> MagicMock:
    """Create a mock parameter with ``__str__`` returning *name* and a numeric *value*."""
    param = MagicMock()
    param.__str__ = lambda self: name
    param.value = value
    return param


# ==================================================================
# _validate_estimation_inputs
# ==================================================================


class TestValidateEstimationInputs:
    """Tests for ``_validate_estimation_inputs``."""

    def test_rejects_non_collection(self) -> None:
        with pytest.raises(ParameterError, match="must be inside a list"):
            _validate_estimation_inputs("not_a_list", None)

    def test_rejects_single_param_object(self) -> None:
        param = _make_param("k1", 5.0)
        with pytest.raises(ParameterError, match="must be inside a list"):
            _validate_estimation_inputs(param, None)

    @pytest.mark.parametrize("container", [list, tuple, set])
    def test_accepts_list_tuple_set(
        self, container: type[list[Any] | tuple[Any, ...] | set[Any]]
    ) -> None:
        param = _make_param("k1", 10.0)
        original, converted, bound = _validate_estimation_inputs(
            container([param]), None
        )
        assert "k1" in converted
        assert isinstance(bound, dict)

    def test_auto_sets_bounds_when_none(self) -> None:
        p1 = _make_param("k1", 100.0)
        p2 = _make_param("k2", 0.5)
        _, converted, bound = _validate_estimation_inputs([p1, p2], None)

        assert converted == ["k1", "k2"]
        assert isinstance(bound, dict)
        assert bound["k1"] == pytest.approx([0.1, 100_000.0])
        assert bound["k2"] == pytest.approx([0.0005, 500.0])

    def test_preserves_explicit_bound_dict(self) -> None:
        p = _make_param("k1", 1.0)
        explicit_bound: dict[str, list[float]] = {"k1": [0.01, 100.0]}
        original, converted, bound = _validate_estimation_inputs([p], explicit_bound)

        assert bound == {"k1": [0.01, 100.0]}
        assert original == [p]

    def test_preserves_explicit_bound_list(self) -> None:
        p = _make_param("k1", 1.0)
        explicit_bound: list[float] = [0.01, 100.0]
        _, _, bound = _validate_estimation_inputs([p], explicit_bound)

        assert bound == [0.01, 100.0]

    def test_converted_parameters_are_strings(self) -> None:
        p1 = _make_param("alpha", 1.0)
        p2 = _make_param("beta", 2.0)
        _, converted, _ = _validate_estimation_inputs([p1, p2], None)

        assert all(isinstance(c, str) for c in converted)
        assert converted == ["alpha", "beta"]


# ==================================================================
# _resolve_experimental_data
# ==================================================================


class TestResolveExperimentalData:
    """Tests for ``_resolve_experimental_data``."""

    def test_uses_sim_experimental_data_when_present(self) -> None:
        df = pd.DataFrame({"time": [0, 1], "A": [10, 5]})
        sim = MagicMock()
        sim.experimental_data.return_pandas.return_value = [df]

        result = _resolve_experimental_data(sim, None)
        assert isinstance(result, pd.DataFrame)
        pd.testing.assert_frame_equal(result, df)

    def test_uses_argument_when_sim_data_is_none(self) -> None:
        sim = MagicMock()
        sim.experimental_data = None
        df = pd.DataFrame({"time": [0, 1], "B": [0, 10]})

        result = _resolve_experimental_data(sim, df)
        pd.testing.assert_frame_equal(result, df)

    def test_raises_when_no_data_at_all(self) -> None:
        sim = MagicMock()
        sim.experimental_data = None

        with pytest.raises(ParameterError, match="No experimental data found"):
            _resolve_experimental_data(sim, None)

    def test_raises_for_invalid_type(self) -> None:
        sim = MagicMock()
        sim.experimental_data = None

        with pytest.raises(ParameterError, match="list of pandas dataframes"):
            _resolve_experimental_data(sim, "not_a_dataframe")

    def test_accepts_list_of_dataframes(self) -> None:
        sim = MagicMock()
        sim.experimental_data = None
        dfs = [
            pd.DataFrame({"time": [0], "A": [1]}),
            pd.DataFrame({"time": [0], "A": [2]}),
        ]

        result = _resolve_experimental_data(sim, dfs)
        assert result == dfs

    def test_accepts_tuple_of_dataframes(self) -> None:
        sim = MagicMock()
        sim.experimental_data = None
        dfs = (pd.DataFrame({"time": [0], "A": [1]}),)

        result = _resolve_experimental_data(sim, dfs)
        assert result == dfs

    def test_sim_data_takes_precedence_over_argument(self) -> None:
        sim_df = pd.DataFrame({"time": [0], "X": [42]})
        arg_df = pd.DataFrame({"time": [0], "Y": [99]})

        sim = MagicMock()
        sim.experimental_data.return_pandas.return_value = [sim_df]

        result = _resolve_experimental_data(sim, arg_df)
        pd.testing.assert_frame_equal(result, sim_df)


# ==================================================================
# _validate_bound
# ==================================================================


class TestValidateBound:
    """Tests for ``_validate_bound``."""

    def test_dict_bound_passes_when_all_params_present(self) -> None:
        bound: dict[str, list[float]] = {"k1": [0.1, 10.0], "k2": [1.0, 100.0]}
        result = _validate_bound(bound, ["k1", "k2"])
        assert result == {"k1": [0.1, 10.0], "k2": [1.0, 100.0]}

    def test_dict_bound_raises_when_param_missing(self) -> None:
        bound: dict[str, list[float]] = {"k1": [0.1, 10.0]}
        with pytest.raises(ParameterError, match="all parameters range"):
            _validate_bound(bound, ["k1", "k2"])

    def test_list_bound_of_length_two_passes(self) -> None:
        result = _validate_bound([0.01, 100.0], ["k1"])
        assert result == [0.01, 100.0]

    def test_tuple_bound_of_length_two_passes(self) -> None:
        result = _validate_bound((0.01, 100.0), ["k1"])
        assert result == (0.01, 100.0)

    @pytest.mark.parametrize(
        "bad_bound",
        [
            [1.0],
            [1.0, 2.0, 3.0],
            [],
        ],
        ids=["too_short", "too_long", "empty"],
    )
    def test_list_bound_wrong_length_raises(self, bad_bound: list[float]) -> None:
        with pytest.raises(ParameterError, match="lower and upper bound"):
            _validate_bound(bad_bound, ["k1"])

    def test_non_iterable_bound_raises(self) -> None:
        with pytest.raises(ParameterError, match="lower and upper bound"):
            _validate_bound(42, ["k1"])

    def test_dict_keys_are_stringified(self) -> None:
        """Non-string keys in bound dict are converted to strings."""
        param = _make_param("k1", 1.0)
        bound = {param: [0.1, 10.0]}
        result = _validate_bound(bound, ["k1"])
        assert "k1" in result

    def test_dict_bound_empty_params_list(self) -> None:
        """An empty parameter list with an empty dict passes."""
        result = _validate_bound({}, [])
        assert result == {}


# ==================================================================
# _extract_results
# ==================================================================


class TestExtractResults:
    """Tests for ``_extract_results``."""

    def test_basic_extraction(self) -> None:
        basico_results = MagicMock()
        basico_results.to_dict.return_value = {
            "sol": {"Values[k1]": 1.5, "Values[k2]": 0.3}
        }

        results = _extract_results(basico_results)
        assert results == {"k1": 1.5, "k2": 0.3}

    def test_single_parameter(self) -> None:
        basico_results = MagicMock()
        basico_results.to_dict.return_value = {"sol": {"Values[alpha]": 42.0}}

        results = _extract_results(basico_results)
        assert results == {"alpha": 42.0}

    def test_empty_results(self) -> None:
        basico_results = MagicMock()
        basico_results.to_dict.return_value = {"sol": {}}

        results = _extract_results(basico_results)
        assert results == {}

    def test_preserves_numeric_types(self) -> None:
        basico_results = MagicMock()
        basico_results.to_dict.return_value = {
            "sol": {
                "Values[int_param]": 5,
                "Values[float_param]": 3.14,
            }
        }

        results = _extract_results(basico_results)
        assert results["int_param"] == 5
        assert results["float_param"] == pytest.approx(3.14)

    def test_parameter_name_with_underscores(self) -> None:
        basico_results = MagicMock()
        basico_results.to_dict.return_value = {
            "sol": {"Values[my_long_param_name]": 0.001}
        }

        results = _extract_results(basico_results)
        assert "my_long_param_name" in results
        assert results["my_long_param_name"] == pytest.approx(0.001)
