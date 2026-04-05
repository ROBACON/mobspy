"""Simulation object utilities."""

from __future__ import annotations

from typing import Any


def sim_remove_reaction(sim: Any, reaction: Any, Simulation_Constructor: Any) -> Any:  # noqa: N803
    """Create a new simulation with the given reaction removed."""
    new_sim = Simulation_Constructor(sim.model)
    new_sim._reactions_set.remove(reaction)
    return new_sim
