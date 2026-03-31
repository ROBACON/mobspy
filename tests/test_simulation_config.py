"""Tests for SimulationConfig dataclass."""

from __future__ import annotations

import pytest

from mobspy.simulation_config import SimulationConfig


class TestDefaults:
    def test_default_values(self):
        config = SimulationConfig()
        assert config.duration == 60
        assert config.volume == 1
        assert config.repetitions == 1
        assert config.simulation_method == "deterministic"
        assert config.level == 3
        assert config.jobs == -1
        assert config.save_data is False
        assert config.plot_data is True

    def test_custom_values(self):
        config = SimulationConfig(duration=200, volume=5, repetitions=3)
        assert config.duration == 200
        assert config.volume == 5
        assert config.repetitions == 3


class TestValidation:
    def test_invalid_simulation_method(self):
        with pytest.raises(ValueError, match="Invalid simulation_method"):
            SimulationConfig(simulation_method="invalid")

    def test_valid_methods(self):
        for method in (
            "deterministic",
            "stochastic",
            "hybrid",
            "hybridode45",
            "hybridlsoda",
            "tauleap",
            "directmethod",
            "sde",
            "lsoda",
        ):
            config = SimulationConfig(simulation_method=method)
            assert config.simulation_method == method

    def test_repetitions_less_than_one(self):
        with pytest.raises(ValueError, match="repetitions must be >= 1"):
            SimulationConfig(repetitions=0)

    def test_invalid_level_too_high(self):
        with pytest.raises(ValueError, match="level must be 0-3"):
            SimulationConfig(level=4)

    def test_invalid_level_negative(self):
        with pytest.raises(ValueError, match="level must be 0-3"):
            SimulationConfig(level=-1)

    def test_volume_zero(self):
        with pytest.raises(ValueError, match="volume must be > 0"):
            SimulationConfig(volume=0)

    def test_volume_negative(self):
        with pytest.raises(ValueError, match="volume must be > 0"):
            SimulationConfig(volume=-1)

    def test_negative_duration(self):
        with pytest.raises(ValueError, match="duration must be >= 0"):
            SimulationConfig(duration=-5)

    def test_zero_duration_allowed(self):
        config = SimulationConfig(duration=0)
        assert config.duration == 0

    def test_invalid_r_tol(self):
        with pytest.raises(ValueError, match="r_tol must be > 0"):
            SimulationConfig(r_tol=0)

    def test_invalid_a_tol(self):
        with pytest.raises(ValueError, match="a_tol must be > 0"):
            SimulationConfig(a_tol=-1e-10)


class TestDictCompat:
    def test_getitem(self):
        config = SimulationConfig(duration=99)
        assert config["duration"] == 99

    def test_getitem_missing(self):
        config = SimulationConfig()
        with pytest.raises(KeyError):
            config["nonexistent"]

    def test_setitem(self):
        config = SimulationConfig()
        config["duration"] = 200
        assert config.duration == 200

    def test_setitem_unknown_key(self):
        config = SimulationConfig()
        with pytest.raises(KeyError, match="Unknown parameter"):
            config["nonexistent"] = 42

    def test_contains(self):
        config = SimulationConfig()
        assert "duration" in config
        assert "nonexistent" not in config
        assert 42 not in config

    def test_get_existing(self):
        config = SimulationConfig(volume=10)
        assert config.get("volume") == 10

    def test_get_missing_default(self):
        config = SimulationConfig()
        assert config.get("nonexistent", "fallback") == "fallback"

    def test_keys(self):
        config = SimulationConfig()
        k = config.keys()
        assert "duration" in k
        assert "volume" in k

    def test_items(self):
        config = SimulationConfig(duration=42)
        items = dict(config.items())
        assert items["duration"] == 42

    def test_values(self):
        config = SimulationConfig()
        vals = config.values()
        assert 60 in vals  # default duration

    def test_to_dict(self):
        config = SimulationConfig(duration=77, volume=3)
        d = config.to_dict()
        assert isinstance(d, dict)
        assert d["duration"] == 77
        assert d["volume"] == 3

    def test_iter(self):
        config = SimulationConfig()
        keys = list(config)
        assert "duration" in keys
        assert "volume" in keys

    def test_update(self):
        config = SimulationConfig()
        config.update({"duration": 500, "volume": 10})
        assert config.duration == 500
        assert config.volume == 10

    def test_update_skips_comments(self):
        config = SimulationConfig()
        config.update({"__comment": "ignored", "duration": 123})
        assert config.duration == 123

    def test_from_dict(self):
        config = SimulationConfig.from_dict(
            {"duration": 300, "unknown_key": "ignored", "volume": 7}
        )
        assert config.duration == 300
        assert config.volume == 7
