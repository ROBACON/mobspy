"""Regression tests for unit conversion with non-default simulation volumes.

BasiCO/COPASI always reports species trajectories as a concentration (amount /
declared compartment volume), never as raw counts. ``convert_data_to_desired_unit``
in ``mobspy/results/process_data.py`` must either leave that concentration as-is
(converted to the requested ``unit_y``) or multiply it back by the volume to
recover counts. For any volume other than the model's native unit, both paths
used to return values off by a volume-dependent factor.
"""

from __future__ import annotations

import pytest

from mobspy import BaseSpecies, Simulation, Zero, u

N = 100  # items placed in the compartment


def _run_at_t0(n_items, volume, *, unit_y=None, output_concentration=None):
    """Place n_items in volume, run a near-zero decay, return the t=0 value."""
    A = BaseSpecies(["A"])
    A >> Zero @ (1e-6 / u.min)  # near-zero decay so the initial value is stable
    A(n_items)
    S = Simulation(A)
    S.run(
        volume=volume,
        duration=0.01 * u.min,
        plot_data=False,
        level=-1,
        unit_y=unit_y,
        output_concentration=output_concentration,
    )
    return S.results["A"][0][0]


@pytest.mark.slow
class TestConcentrationOutputUnits:
    """N items placed in volume V must report as N / V in the requested unit_y."""

    @pytest.mark.parametrize(
        ("volume", "unit_y", "expected"),
        [
            (1.0 * u.L, 1 / u.L, 100),  # 100 items / 1 L     = 100 item/L
            (0.001 * u.L, 1 / u.L, 100_000),  # 100 items / 0.001 L = 100 000 item/L
            (1.0 * u.mL, 1 / u.mL, 100),  # 100 items / 1 mL    = 100 item/mL
            (0.2 * u.mL, 1 / u.mL, 500),  # 100 items / 0.2 mL  = 500 item/mL
        ],
    )
    def test_concentration_output_matches_expected(self, volume, unit_y, expected):
        got = _run_at_t0(N, volume, unit_y=unit_y)
        assert abs(got - expected) / expected < 0.01, (
            f"volume={volume}, unit_y={unit_y}: got {got:.2f}, expected {expected}"
        )

    def test_default_output_concentration_matches_native_volume_unit(self):
        # No explicit unit_y/output_concentration: falls back to the model's
        # native (declared) volume unit, still N / V.
        got = _run_at_t0(N, 0.2 * u.mL)
        expected = 500  # 100 items / 0.2 mL
        assert abs(got - expected) / expected < 0.01, f"got {got:.2f}, expected {expected}"


@pytest.mark.slow
class TestCountOutputUnits:
    """Requesting raw counts must return the true item count, independent of volume."""

    @pytest.mark.parametrize(
        "volume",
        [1.0 * u.L, 1.0 * u.mL, 0.2 * u.mL, 0.001 * u.L],
    )
    def test_count_output_matches_initial_count(self, volume):
        got = _run_at_t0(N, volume, output_concentration=False)
        assert abs(got - N) / N < 0.01, f"volume={volume}: got {got:.2f}, expected {N}"
