"""Test fractional product handling in Mobspy."""

import pytest
import mobspy
from mobspy import (
    BaseSpecies,
    Simulation,
    New,
    u,
    ModelParameters,
    simlog,
    Assign,
    Rev,
    All,
    set_counts,
)


def test_deterministic():
    """Test fractional products in deterministic simulation."""
    A, B, C = BaseSpecies()
    A(10)
    B(1)
    A + B >> 0.5 * C + B [0.005]
    Sim = Simulation(A | B | C)
    Sim.method = "deterministic"
    Sim.volume = 1 * u.mL
    Sim.run(duration=10, unit_x=u.s, unit_y=1 / u.mL, plot_data=False)
    
    # approximate check
    assert Sim.results["C"][0][-1] == pytest.approx(5)


def test_stochastic():
    """Test fractional products in stochastic simulation."""
    A, B, C = BaseSpecies()
    A(10)
    B(1)
    A + B >> 0.5 * C + B [0.005]
    Sim = Simulation(A | B | C)
    Sim.method = "stochastic"
    Sim.volume = 1 * u.mL
    Sim.run(duration=10, unit_x=u.s, unit_y=1 / u.mL, plot_data=False)
    
    # approximate check
    assert Sim.results["C"][0][-1] == pytest.approx(5)


if __name__ == "__main__":
    test_deterministic()
    test_stochastic()