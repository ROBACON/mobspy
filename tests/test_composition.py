"""Tests for SimulationComposition (concatenating simulations with +)."""

from __future__ import annotations

import pytest

import mobspy
from mobspy import BaseSpecies, Simulation, SimulationComposition, simlog
from mobspy.exceptions import SimulationError

simlog.global_simlog_level = -1


class TestBasicComposition:
    def test_two_simulations(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        A(100)
        S1 = Simulation(A)
        S1.duration = 3

        S2 = Simulation(A)
        S2.duration = 5

        S = S1 + S2
        assert isinstance(S, SimulationComposition)
        assert len(S) == 2

    def test_three_simulations(self):
        A, B, C = BaseSpecies(3)
        A >> mobspy.Zero[1]
        A(50)
        S1 = Simulation(A)
        S1.duration = 2

        B >> mobspy.Zero[1]
        B(50)
        S2 = Simulation(B)
        S2.duration = 3

        C >> mobspy.Zero[1]
        C(50)
        S3 = Simulation(C)
        S3.duration = 4

        S = S1 + S2 + S3
        assert len(S) == 3

    def test_composition_is_iterable(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        sims = list(S)
        assert sims[0] is S1
        assert sims[1] is S2

    def test_base_sim_is_first(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        assert S.base_sim is S1

    def test_is_simulation_classmethod(self):
        assert SimulationComposition.is_simulation() is True


class TestBroadcastParameters:
    """Broadcast parameters are set identically on all child simulations."""

    def test_level_broadcast(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        S.level = -1
        assert S1.__dict__["parameters"]["level"] == -1
        assert S2.__dict__["parameters"]["level"] == -1

    def test_repetitions_broadcast(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        S.repetitions = 5
        assert S1.__dict__["parameters"]["repetitions"] == 5
        assert S2.__dict__["parameters"]["repetitions"] == 5


class TestMulticastParameters:
    """Multicast parameters require a list matching the number of simulations."""

    def test_duration_multicast(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S2 = Simulation(A)

        S = S1 + S2
        S.duration = [3, 7]
        assert S1.__dict__["parameters"]["duration"] == 3
        assert S2.__dict__["parameters"]["duration"] == 7

    def test_duration_scalar_raises(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S2 = Simulation(A)

        S = S1 + S2
        with pytest.raises(SimulationError):
            S.duration = 5

    def test_duration_wrong_length_raises(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S2 = Simulation(A)

        S = S1 + S2
        with pytest.raises(SimulationError):
            S.duration = [1, 2, 3]


class TestDoublecastParameters:
    """Double-cast parameters broadcast on scalar, multicast on list."""

    def test_volume_scalar_broadcasts(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        S.volume = 10
        assert S1.__dict__["parameters"]["volume"] == 10
        assert S2.__dict__["parameters"]["volume"] == 10

    def test_volume_list_multicasts(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        S.volume = [5, 20]
        assert S1.__dict__["parameters"]["volume"] == 5
        assert S2.__dict__["parameters"]["volume"] == 20

    def test_volume_wrong_length_raises(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        with pytest.raises(SimulationError):
            S.volume = [1, 2, 3]

    def test_method_scalar_broadcasts(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        S.method = "stochastic"
        assert S1.__dict__["parameters"]["method"] == "stochastic"
        assert S2.__dict__["parameters"]["method"] == "stochastic"

    def test_simulation_method_list_multicasts(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        S.simulation_method = ["deterministic", "stochastic"]
        assert S1.__dict__["parameters"]["simulation_method"] == "deterministic"
        assert S2.__dict__["parameters"]["simulation_method"] == "stochastic"


class TestPlotConfig:
    def test_plot_config_delegates_to_base_sim(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        S.plot_config.title = "Test Title"
        assert S1.plot_config.title == "Test Title"
        assert S.plot_config.title == "Test Title"

    def test_plot_config_shares_state_with_base_sim(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        S.plot_config.xlabel = "Time (s)"
        assert S1.plot_config.xlabel == "Time (s)"


class TestCompilation:
    def test_compile_returns_string(self):
        A = BaseSpecies()
        A >> mobspy.Zero[1]
        A(100)
        S1 = Simulation(A)
        S1.duration = 3

        S2 = Simulation(A)
        S2.duration = 5

        S = S1 + S2
        result = S.compile()
        assert isinstance(result, str)
        assert len(result) > 0

    def test_compile_with_different_species(self):
        A, B = BaseSpecies(2)
        A >> mobspy.Zero[1]
        A(50)
        S1 = Simulation(A)
        S1.duration = 3

        B >> mobspy.Zero[0.5]
        B(100)
        S2 = Simulation(A | B)
        S2.duration = 5

        S = S1 + S2
        result = S.compile()
        assert result is not None

    def test_compile_detects_modified_species(self):
        """Species modified between simulations should raise an error."""
        A = BaseSpecies()
        S1 = Simulation(A)
        S1.duration = 1

        # Re-create A with same name but different characteristics
        A = BaseSpecies()
        A.a1, A.a2
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        with pytest.raises((SystemExit, SimulationError)):
            S.compile()


class TestErrorCases:
    def test_add_non_simulation_raises(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        with pytest.raises(SimulationError):
            SimulationComposition(S1, 42)  # type: ignore[arg-type]

    def test_multicast_duration_length_mismatch(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S2 = Simulation(A)
        S3 = Simulation(A)

        S = S1 + S2 + S3
        with pytest.raises(SimulationError):
            S.duration = [1, 2]

    def test_simulation_method_list_length_mismatch(self):
        A = BaseSpecies()
        A(10)
        S1 = Simulation(A)
        S1.duration = 1
        S2 = Simulation(A)
        S2.duration = 2

        S = S1 + S2
        with pytest.raises(SimulationError):
            S.simulation_method = ["deterministic", "stochastic", "deterministic"]
