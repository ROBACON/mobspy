"""Typed configuration dataclass for MobsPy simulations."""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Any


@dataclass
class SimulationConfig:
    """Typed configuration for MobsPy simulations.

    Provides IDE autocompletion, validation, and type safety for all
    simulation parameters. Supports dict-style access for backward
    compatibility.

    Examples:
        >>> config = SimulationConfig()
        >>> config.duration
        60
        >>> config["volume"]
        1
        >>> config.simulation_method
        'deterministic'
        >>> config.update({"duration": 100, "volume": 2})
        >>> config.duration
        100
        >>> "duration" in config
        True
        >>> config.to_dict()["volume"]
        2
    """

    # Model parameters
    volume: int | float = 1
    repetitions: int = 1
    level: int = 3
    rate_type: str | None = None

    # Simulation engine parameters
    simulation_method: str = "deterministic"
    method: str | None = None
    start_time: int | float = 0
    duration: int | float = 60
    r_tol: float = 1e-8
    a_tol: float = 1e-10
    step_size: int | float | None = None
    seeds: list[int] | None = None

    # Parallelism
    jobs: int = -1

    # Output configuration
    output_dir: str = "outputs/"
    output_file: str | None = None
    output_event: bool = False
    unit_x: Any = None
    unit_y: Any = None
    skip_expression_check: bool = False
    output_concentration: bool = True
    save_data: bool = False
    plot_data: bool = True
    plot_type: str | None = None
    absolute_output_file: str = ""

    # Internal
    _continuous_simulation: bool = False
    _end_condition: Any = None
    _with_event: bool = False
    initial_conditional_duration: int | float = 0

    def __post_init__(self) -> None:
        """Validate configuration values."""
        _valid_methods = {
            "deterministic",
            "stochastic",
            "hybrid",
            "hybridode45",
            "hybridlsoda",
            "tauleap",
            "directmethod",
            "sde",
            "lsoda",
        }
        if self.simulation_method not in _valid_methods:
            msg = (
                f"Invalid simulation_method: {self.simulation_method!r}. "
                "Must be 'deterministic', 'stochastic', or a valid BasiCO method."
            )
            raise ValueError(msg)
        if self.repetitions < 1:
            msg = f"repetitions must be >= 1, got {self.repetitions}"
            raise ValueError(msg)
        if not (0 <= self.level <= 3):
            msg = f"level must be 0-3, got {self.level}"
            raise ValueError(msg)
        if self.volume <= 0:
            msg = f"volume must be > 0, got {self.volume}"
            raise ValueError(msg)
        if self.duration < 0:
            msg = f"duration must be >= 0, got {self.duration}"
            raise ValueError(msg)
        if self.r_tol <= 0:
            msg = f"r_tol must be > 0, got {self.r_tol}"
            raise ValueError(msg)
        if self.a_tol <= 0:
            msg = f"a_tol must be > 0, got {self.a_tol}"
            raise ValueError(msg)

    # --- Dict-style access for backward compatibility ---

    def __getitem__(self, key: str) -> Any:
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(key) from None

    def __setitem__(self, key: str, value: Any) -> None:
        if hasattr(self, key):
            object.__setattr__(self, key, value)
        else:
            raise KeyError(f"Unknown parameter: {key!r}")

    def __contains__(self, key: object) -> bool:
        if not isinstance(key, str):
            return False
        return key in {f.name for f in fields(self)}

    def get(self, key: str, default: Any = None) -> Any:
        """Dict-compatible get with default."""
        try:
            return getattr(self, key)
        except AttributeError:
            return default

    def keys(self) -> list[str]:
        """Return parameter names."""
        return [f.name for f in fields(self)]

    def items(self) -> list[tuple[str, Any]]:
        """Dict-like items() for backward compat."""
        return [(f.name, getattr(self, f.name)) for f in fields(self)]

    def values(self) -> list[Any]:
        """Dict-like values() for backward compat."""
        return [getattr(self, f.name) for f in fields(self)]

    def to_dict(self) -> dict[str, Any]:
        """Convert to plain dict."""
        return {f.name: getattr(self, f.name) for f in fields(self)}

    def update(self, other: dict[str, Any]) -> None:
        """Merge values from a dict, skipping comment keys."""
        for key, value in other.items():
            if key.startswith("__comment"):
                continue
            if hasattr(self, key):
                object.__setattr__(self, key, value)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> SimulationConfig:
        """Create a SimulationConfig from a dict, ignoring unknown keys."""
        valid_keys = {f.name for f in fields(cls)}
        filtered = {k: v for k, v in data.items() if k in valid_keys}
        return cls(**filtered)

    def __iter__(self) -> Any:
        """Allow iteration over keys (for `for key in config` patterns)."""
        return iter(self.keys())
