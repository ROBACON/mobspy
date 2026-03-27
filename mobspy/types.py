"""Type definitions for MobsPy.

Provides TypedDicts and type aliases used throughout the
MobsPy package for structured, type-safe dictionary access.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, TypedDict

# --- SBML data structures (compiler -> builder -> SBMLWriter) ---


class ReactionData(TypedDict):
    """Single reaction in SBML format."""

    re: list[tuple[float | int, str]]
    pr: list[tuple[float | int, str]]
    kin: str


class EventData(TypedDict):
    """Single event in SBML format."""

    trigger: str
    delay: str | float | int
    assignments: list[tuple[str, str | int | float]]


class AssignmentData(TypedDict):
    """Single assignment in SBML format."""

    species: str
    expression: str


class ParameterUsedInfo(TypedDict):
    """Metadata for a parameter tracked during compilation."""

    name: str
    values: float | int | list[float | int]
    used_in: set[str]
    object: Any


# --- Composite dict structures ---


class SBMLModelDict(TypedDict):
    """Dict passed to builder.build() / SBMLWriter."""

    species_for_sbml: dict[str, int | float]
    parameters_for_sbml: dict[str, tuple[float | int, str]]
    reactions_for_sbml: dict[str, ReactionData]
    events_for_sbml: dict[str, EventData]
    assignments_for_sbml: dict[str, AssignmentData]


class CompiledModelDict(TypedDict):
    """Full compiled model stored in _list_of_models."""

    species_for_sbml: dict[str, int | float]
    parameters_for_sbml: dict[str, tuple[float | int, str]]
    reactions_for_sbml: dict[str, ReactionData]
    events_for_sbml: dict[str, EventData]
    assignments_for_sbml: dict[str, AssignmentData]
    species_not_mapped: dict[str, int | float]
    mappings: dict[str, list[str]]
    assigned_species: list[str]


class SimulationEventData(TypedDict):
    """Internal event data collected during event context."""

    event_time: float
    event_counts: list[Any]
    trigger: str


class TimeSeriesDataDict(TypedDict):
    """Data structure passed to MobsPyTimeSeries constructor."""

    data: dict[str, list[float]]
    params: SimulationParameters
    models: list[CompiledModelDict]


# --- Simulation parameters ---


class SimulationParameters(TypedDict, total=False):
    """Parameters dict used throughout Simulation (from default_reader).

    Uses total=False because comment keys and optional keys coexist.
    """

    volume: float | int
    repetitions: int
    level: int
    rate_type: str | None
    simulation_method: str
    method: str | None
    start_time: float | int
    duration: float | int
    r_tol: float
    a_tol: float
    jobs: int
    output_dir: str
    output_file: str | None
    output_event: bool
    unit_x: str | None
    unit_y: str | None
    skip_expression_check: bool
    output_concentration: bool
    save_data: bool
    plot_data: bool
    plot_type: str | None
    _continuous_simulation: bool
    _end_condition: str | None
    _with_event: bool


# --- Type aliases ---

SpeciesForSbml = dict[str, int | float]
ParameterValue = tuple[float | int, str]
ParametersForSbml = dict[str, ParameterValue]
ReactionsForSbml = dict[str, ReactionData]
EventsForSbml = dict[str, EventData]
AssignmentsForSbml = dict[str, AssignmentData]
MappingsForSbml = dict[str, list[str]]
ParametersUsed = dict[str, ParameterUsedInfo]
ParameterSweepList = list[list[CompiledModelDict]]
SimParams = SimulationParameters


# --- Compiler result ---


@dataclass
class CompilerResult:
    """Structured result from Compiler.compile().

    Replaces the 11-element tuple previously returned.
    """

    species_for_sbml: dict[str, int | float] = field(default_factory=dict)
    reactions_for_sbml: dict[str, ReactionData] = field(default_factory=dict)
    parameters_for_sbml: dict[str, tuple[float | int, str]] = field(
        default_factory=dict
    )
    mappings_for_sbml: dict[str, list[str]] = field(default_factory=dict)
    model_str: str = ""
    events_for_sbml: dict[str, EventData] = field(default_factory=dict)
    assigned_species: list[str] = field(default_factory=list)
    parameters_used: dict[str, ParameterUsedInfo] = field(default_factory=dict)
    parameter_object_dict: dict[str, Any] = field(default_factory=dict)
    assignments_for_sbml: dict[str, AssignmentData] = field(default_factory=dict)
    has_mole: bool = False

    def to_compiled_model_dict(
        self,
        species_not_mapped: dict[str, int | float],
        mappings: dict[str, list[str]],
    ) -> CompiledModelDict:
        """Create a CompiledModelDict from compiler output."""
        return {
            "species_for_sbml": self.species_for_sbml,
            "parameters_for_sbml": self.parameters_for_sbml,
            "reactions_for_sbml": self.reactions_for_sbml,
            "events_for_sbml": self.events_for_sbml,
            "assignments_for_sbml": self.assignments_for_sbml,
            "species_not_mapped": species_not_mapped,
            "mappings": mappings,
            "assigned_species": self.assigned_species,
        }
