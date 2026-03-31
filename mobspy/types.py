"""Type definitions for MobsPy.

Provides dataclasses and type aliases used throughout the
MobsPy package for structured, type-safe data access.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# --- SBML data structures (compiler -> builder -> sbml_writer) ---


@dataclass
class ReactionData:
    """Single reaction in SBML format.

    Examples:
        >>> r = ReactionData()
        >>> r.reactants
        []
        >>> r.kinetics
        ''
    """

    reactants: list[tuple[float | int, str]] = field(default_factory=list)
    products: list[tuple[float | int, str]] = field(default_factory=list)
    kinetics: str = ""

    def __str__(self) -> str:
        """Backward-compatible dict-style string representation."""
        return str({"re": self.reactants, "pr": self.products, "kin": self.kinetics})


@dataclass
class EventData:
    """Single event in SBML format."""

    trigger: str = ""
    delay: str | float | int = "0"
    assignments: list[tuple[str, str | int | float]] = field(default_factory=list)

    def __str__(self) -> str:
        """Backward-compatible dict-style string representation."""
        return str(
            {
                "trigger": self.trigger,
                "delay": self.delay,
                "assignments": self.assignments,
            }
        )


@dataclass
class AssignmentData:
    """Single assignment in SBML format."""

    species: str = ""
    expression: str = ""

    def __str__(self) -> str:
        """Backward-compatible dict-style string representation."""
        return str({"species": self.species, "expression": self.expression})


@dataclass
class ParameterUsedInfo:
    """Metadata for a parameter tracked during compilation."""

    name: str = ""
    values: float | int | list[float | int] = 0
    used_in: set[str] = field(default_factory=set)
    object: Any = None


# --- Composite structures ---


@dataclass
class SBMLModelData:
    """Data passed to builder.build() / sbml_writer."""

    species_for_sbml: dict[str, int | float] = field(default_factory=dict)
    parameters_for_sbml: dict[str, tuple[float | int, str]] = field(
        default_factory=dict
    )
    reactions_for_sbml: dict[str, ReactionData] = field(default_factory=dict)
    events_for_sbml: dict[str, EventData] = field(default_factory=dict)
    assignments_for_sbml: dict[str, AssignmentData] = field(default_factory=dict)

    # Backward compatibility: dict-like access for gradual migration
    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        setattr(self, key, value)


@dataclass
class CompiledModel:
    """Full compiled model stored in _list_of_models.

    Examples:
        >>> m = CompiledModel()
        >>> m.species_for_sbml
        {}
        >>> m["species_for_sbml"]
        {}
    """

    species_for_sbml: dict[str, int | float] = field(default_factory=dict)
    parameters_for_sbml: dict[str, tuple[float | int, str]] = field(
        default_factory=dict
    )
    reactions_for_sbml: dict[str, ReactionData] = field(default_factory=dict)
    events_for_sbml: dict[str, EventData] = field(default_factory=dict)
    assignments_for_sbml: dict[str, AssignmentData] = field(default_factory=dict)
    species_not_mapped: dict[str, int | float] = field(default_factory=dict)
    mappings: dict[str, list[str]] = field(default_factory=dict)
    assigned_species: list[str] = field(default_factory=list)
    model_context: Any = None  # ModelUnitContext | None (avoid circular import)

    # Backward compatibility: dict-like access for gradual migration
    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        setattr(self, key, value)

    def __contains__(self, key: str) -> bool:
        return hasattr(self, key)

    def items(self) -> list[tuple[str, Any]]:
        """Dict-like items() for backward compat."""
        return [
            (f.name, getattr(self, f.name)) for f in self.__dataclass_fields__.values()
        ]


@dataclass
class SimulationEventData:
    """Internal event data collected during event context."""

    event_time: Any = 0.0
    event_counts: list[Any] = field(default_factory=list)
    trigger: str = ""


class TimeSeriesDataDict:
    """Data structure passed to MobsPyTimeSeries constructor."""

    def __init__(
        self,
        data: dict[str, list[float]] | None = None,
        params: SimulationParameters | Any | None = None,
        models: list[CompiledModel] | None = None,
    ) -> None:
        self.data: dict[str, list[float]] = data or {}
        self.params: SimulationParameters | Any | None = params
        self.models: list[CompiledModel] = models or []

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        setattr(self, key, value)


# --- Simulation parameters ---
# Kept as TypedDict because it's used as a loose config bag with
# optional keys, comment keys, and dict-spread patterns.

from typing import TypedDict  # noqa: E402


class SimulationParameters(TypedDict):
    """Parameters dict used throughout Simulation (from default_reader)."""

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
    absolute_output_file: str
    initial_conditional_duration: float | int
    step_size: float | None
    seeds: list[int] | None


# --- Type aliases ---

SpeciesForSbml = dict[str, int | float]
ParameterValue = tuple[float | int, str]
ParametersForSbml = dict[str, ParameterValue]
ReactionsForSbml = dict[str, ReactionData]
EventsForSbml = dict[str, EventData]
AssignmentsForSbml = dict[str, AssignmentData]
MappingsForSbml = dict[str, list[str]]
ParametersUsed = dict[str, ParameterUsedInfo]
ParameterSweepList = list[list[CompiledModel]]
SimParams = SimulationParameters

# Backward compat aliases for the old TypedDict names
SBMLModelDict = SBMLModelData
CompiledModelDict = CompiledModel


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
    model_context: Any = None  # ModelUnitContext | None (avoid circular import)

    def to_compiled_model(
        self,
        species_not_mapped: dict[str, int | float],
        mappings: dict[str, list[str]],
    ) -> CompiledModel:
        """Create a CompiledModel from compiler output."""
        return CompiledModel(
            species_for_sbml=self.species_for_sbml,
            parameters_for_sbml=self.parameters_for_sbml,
            reactions_for_sbml=self.reactions_for_sbml,
            events_for_sbml=self.events_for_sbml,
            assignments_for_sbml=self.assignments_for_sbml,
            species_not_mapped=species_not_mapped,
            mappings=mappings,
            assigned_species=self.assigned_species,
            model_context=self.model_context,
        )

    # Keep old name as alias
    def to_compiled_model_dict(
        self,
        species_not_mapped: dict[str, int | float],
        mappings: dict[str, list[str]],
    ) -> CompiledModel:
        """Backward compat alias for to_compiled_model."""
        return self.to_compiled_model(species_not_mapped, mappings)


@dataclass
class CompilationContext:
    """Shared context threaded through compiler helper methods.

    Groups parameters that are always passed together during
    compilation: volume, dimension, model type, unit context,
    and the parameter registry.
    """

    volume: int | float = 1
    dimension: int | None = 3
    type_of_model: str = "deterministic"
    model_context: Any = None  # ModelUnitContext | None
    parameter_exist: dict[str, Any] = field(default_factory=dict)
    skip_expression_check: bool = False
    continuous_sim: bool = False
    ending_condition: Any = None
    event_dictionary: list[Any] | None = None


@dataclass
class CountAccumulator:
    """Mutable state accumulated while assigning initial species counts."""

    species_for_sbml: dict[str, int | float] = field(default_factory=dict)
    assigned_species: list[str] = field(default_factory=list)
    parameters_in_counts: set[Any] = field(default_factory=set)
    parameters_used: dict[str, Any] = field(default_factory=dict)


@dataclass
class RenderContext:
    """Flags controlling how species references are resolved to SBML strings."""

    count_in_model: bool = True
    concentration_in_model: bool = False
    count_in_expression: bool = True
    concentration_in_expression: bool = False
