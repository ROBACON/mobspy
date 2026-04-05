"""Resolve meta-species inheritance and expand meta-reactions.

The primary entry point is :func:`compile_model`, a standalone pure
function that orchestrates a typed phased pipeline:

1. **Species setup** -- validate names, enumerate concrete species
2. **Volume resolution** -- resolve volume/dimension from user input
3. **Count assignment** -- assign initial counts to species
4. **Reaction expansion** -- expand meta-reactions to concrete
5. **Duplicate detection** -- warn on duplicate reactions (O(n))
6. **Event building** -- compile events and phantom reactions
7. **Parameter validation** -- check uniqueness and collisions
8. **Assignment building** -- compile ODE assignment rules
9. **Model string** -- generate human-readable output

Each phase is a standalone function returning a typed result.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from pint import Quantity

from mobspy.constants import (
    ALL_CHAR,
    DOT_SEPARATOR,
    END_FLAG_SPECIES_NAME,
    SBML_LOCATION,
    STD_CHAR,
)
from mobspy.exceptions import CompilationError
from mobspy.mobspy_logging import get_logger
from mobspy.modules.assignments_implementation import (
    Assign as asgi_Assign,
)
from mobspy.modules.compiler_operators import (
    create_all_not_reactions as cof_create_all_not_reactions,
)
from mobspy.modules.event_functions import (
    format_event_dictionary_for_sbml as eh_format_event_dictionary_for_sbml,
)
from mobspy.modules.mobspy_parameters import (
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)
from mobspy.modules.reaction_expansion import (
    create_all_reactions as rc_create_all_reactions,
)
from mobspy.modules.species_constructors import EndFlagSpecies
from mobspy.modules.species_string_generator import (
    construct_all_combinations as ssg_construct_all_combinations,
)
from mobspy.modules.species_string_generator import (
    construct_all_species_ids as ssg_construct_all_species_ids,
)
from mobspy.modules.species_string_generator import (
    construct_species_id as ssg_construct_species_id,
)
from mobspy.modules.unit_handler import (
    convert_counts as uh_convert_counts,
)
from mobspy.modules.unit_handler import (
    convert_volume as uh_convert_volume,
)
from mobspy.modules.unit_handler import (
    extract_length_dimension as uh_extract_length_dimension,
)
from mobspy.types import (
    CompilationContext,
    CompilerResult,
    ConcreteModel,
    ConcreteSpeciesId,
    CountAccumulator,
    CountAssignmentResult,
    EventBuildResult,
    EventData,
    ParameterUsedInfo,
    ReactionData,
    ReactionExpansionResult,
    SpeciesSetupResult,
    VolumeResolutionResult,
)

if TYPE_CHECKING:
    from mobspy.modules.list_species import List_Species
    from mobspy.modules.model_unit_context import ModelUnitContext
    from mobspy.types import (
        MappingsForSbml,
        ParametersForSbml,
        ParametersUsed,
        ReactionsForSbml,
        SpeciesForSbml,
    )

_logger = get_logger(__name__)


# ------------------------------------------------------------------
# Phase 1: Species setup -- validate and enumerate
# ------------------------------------------------------------------


def phase_species_setup(
    meta_species_to_simulate: List_Species,
    orthogonal_vector_structure: dict[str, Any],
) -> SpeciesSetupResult:
    """Validate species names and enumerate all concrete species.

    Args:
        meta_species_to_simulate: Deduplicated list of meta-species.
        orthogonal_vector_structure: Characteristic-to-species mapping.

    Returns:
        SpeciesSetupResult with species dict, mappings, and validated names.

    Raises:
        CompilationError: If species names are invalid or duplicated.
    """
    names_used = _validate_species_names(meta_species_to_simulate)
    species, mappings = _build_species_and_mappings(
        meta_species_to_simulate, orthogonal_vector_structure
    )
    return SpeciesSetupResult(
        species=species,
        mappings=mappings,
        names_used=frozenset(names_used),
    )


# ------------------------------------------------------------------
# Phase 2: Volume resolution
# ------------------------------------------------------------------


def phase_volume_resolution(
    volume: int | float | Quantity,
    dimension: int | None,
    species_counts: list[dict[str, Any]],
    model_context: ModelUnitContext | None = None,
) -> VolumeResolutionResult:
    """Resolve volume and spatial dimension from user input.

    Args:
        volume: System volume (possibly with units).
        dimension: Spatial dimension (None for auto-detect).
        species_counts: Count assignments (may contain Pint quantities).
        model_context: Unit context for conversion.

    Returns:
        VolumeResolutionResult with resolved volume, dimension, and
        the volume parameter entry.
    """
    resolved_volume, resolved_dimension, params = _resolve_volume_and_dimension(
        volume, dimension, species_counts, model_context
    )
    vol_unit_id = params["volume"][1]
    return VolumeResolutionResult(
        volume=resolved_volume,
        dimension=resolved_dimension,
        volume_parameter=(resolved_volume, vol_unit_id),
    )


# ------------------------------------------------------------------
# Phase 3: Count assignment
# ------------------------------------------------------------------


def phase_count_assignment(
    species_counts: list[dict[str, Any]],
    species: dict[str, int | float],
    orthogonal_vector_structure: dict[str, Any],
    parameters_used: ParametersUsed,
    ctx: CompilationContext,
) -> CountAssignmentResult:
    """Assign initial counts to concrete species.

    Mutates ``species`` and ``parameters_used`` in place.

    Args:
        species_counts: Raw count declarations.
        species: Species dict to update with counts.
        orthogonal_vector_structure: Characteristic mapping.
        parameters_used: Parameter tracking dict (mutated).
        ctx: Compilation context.

    Returns:
        CountAssignmentResult with assigned species and parameter ids.
    """
    assigned, params_in_counts = _assign_initial_counts(
        species_counts, species, orthogonal_vector_structure, parameters_used, ctx
    )
    return CountAssignmentResult(
        assigned_species=tuple(assigned),
        parameters_in_counts=frozenset(params_in_counts),
    )


# ------------------------------------------------------------------
# Phase 4: Reaction expansion
# ------------------------------------------------------------------


def phase_reaction_expansion(
    reactions_set: set[Any],
    meta_species_to_simulate: List_Species,
    orthogonal_vector_structure: dict[str, Any],
    ctx: CompilationContext,
) -> ReactionExpansionResult:
    """Expand meta-reactions into concrete SBML-ready reactions.

    Args:
        reactions_set: Set of meta-reactions.
        meta_species_to_simulate: Species in the model.
        orthogonal_vector_structure: Characteristic mapping.
        ctx: Compilation context.

    Returns:
        ReactionExpansionResult with concrete reactions and parameter ids.
    """
    reactions, params = _build_reactions(
        reactions_set, meta_species_to_simulate, orthogonal_vector_structure, ctx
    )
    return ReactionExpansionResult(
        reactions=reactions,
        parameters_in_reactions=frozenset(params),
    )


# ------------------------------------------------------------------
# Phase 5: Duplicate detection (O(n) via hashing)
# ------------------------------------------------------------------


def phase_duplicate_detection(
    reactions: ReactionsForSbml,
) -> None:
    """Warn about duplicate reactions.

    Uses hash-based comparison for O(n) performance.

    Args:
        reactions: Concrete reactions dict.
    """
    seen: dict[tuple[tuple[Any, ...], tuple[Any, ...], str], str] = {}
    for name, rxn in reactions.items():
        key = (tuple(rxn.reactants), tuple(rxn.products), rxn.kinetics)
        if key in seen:
            _logger.warning(
                "The following reaction: \n"
                + f"{rxn} \n"
                + "Is doubled. Was that intentional? \n"
            )
        else:
            seen[key] = name


# ------------------------------------------------------------------
# Phase 6: Event building
# ------------------------------------------------------------------


def phase_event_building(
    species: SpeciesForSbml,
    reactions: ReactionsForSbml,
    orthogonal_vector_structure: dict[str, Any],
    meta_species_to_simulate: List_Species,
    ctx: CompilationContext,
) -> EventBuildResult:
    """Compile events and add phantom reactions for event-only species.

    Mutates ``reactions`` by adding phantom reactions.

    Args:
        species: Species dict.
        reactions: Reactions dict (mutated with phantom reactions).
        orthogonal_vector_structure: Characteristic mapping.
        meta_species_to_simulate: Species in the model.
        ctx: Compilation context.

    Returns:
        EventBuildResult with events and parameter ids.
    """
    events, params = _build_events(
        species, reactions, orthogonal_vector_structure, meta_species_to_simulate, ctx
    )
    return EventBuildResult(
        events=events,
        parameters_in_events=frozenset(params),
    )


# ------------------------------------------------------------------
# Phase 7: Parameter validation
# ------------------------------------------------------------------


def phase_parameter_validation(
    parameters: ParametersForSbml,
    names_used: frozenset[str],
    all_param_objects: set[mp_Mobspy_Parameter],
) -> None:
    """Validate parameter names for uniqueness and species collision.

    Args:
        parameters: Parameters dict.
        names_used: Set of species names already in use.
        all_param_objects: All parameter objects from counts, reactions, events.

    Raises:
        CompilationError: On name collision or duplicate parameter names.
    """
    mutable_names = set(names_used)
    for p in parameters:
        if p in mutable_names:
            raise CompilationError(
                "Parameters names must be unique"
                " and they must not share a"
                " name with a species"
            )
        mutable_names.add(p)

    for p1 in all_param_objects:
        for p2 in all_param_objects:
            if p1 is not p2 and p1.name == p2.name:
                raise CompilationError(
                    "There are two different Parameter Objects with the same name"
                )


# ------------------------------------------------------------------
# Phase 8: Assignment building
# ------------------------------------------------------------------


def phase_assignments(
    meta_species_to_simulate: List_Species,
    orthogonal_vector_structure: dict[str, Any],
) -> dict[str, Any]:
    """Compile species assignment rules for SBML.

    Args:
        meta_species_to_simulate: Species in the model.
        orthogonal_vector_structure: Characteristic mapping.

    Returns:
        Dict of compiled assignments.
    """
    return _build_assignments(meta_species_to_simulate, orthogonal_vector_structure)


# ------------------------------------------------------------------
# Private helpers (unchanged logic, used by phases)
# ------------------------------------------------------------------


def _add_to_parameters_to_sbml(
    parameters_used: ParametersUsed,
    parameters_for_sbml: ParametersForSbml,
    parameters_to_add: set[mp_Mobspy_Parameter],
    model_context: ModelUnitContext | None = None,
) -> None:
    """Register model parameters into the SBML dictionaries."""
    param_unit = "dimensionless"
    if model_context is not None:
        rate_unit_id = f"per_{model_context.get_sbml_time_units_id()}"
        param_unit = rate_unit_id

    for parameter in parameters_to_add:
        if parameter.name in parameters_used:
            parameters_used[parameter.name].used_in.add(SBML_LOCATION)
        else:
            parameters_used[parameter.name] = ParameterUsedInfo(
                name=parameter.name,
                values=parameter.value,
                used_in={SBML_LOCATION},
                object=parameter,
            )

        try:
            parameters_for_sbml[parameter.name] = (
                parameter.value[0],
                param_unit,
            )
        except (IndexError, TypeError):
            parameters_for_sbml[parameter.name] = (  # pyright: ignore[reportArgumentType]
                parameter.value,
                param_unit,
            )


def _add_phantom_reactions(
    reactions_for_sbml: ReactionsForSbml,
    species_to_add_phantom_reaction: set[str],
) -> None:
    """Add phantom reactions for event-only species.

    Basico cannot compile species assigned to events
    but not present in reactions. This error does not
    happen in python 3.10. To be compatible with
    earlier versions we add an extremely slow phantom
    reaction for each such species with a coefficient
    of 1e-100.
    """
    for i, spe in enumerate(species_to_add_phantom_reaction):
        reactions_for_sbml["phantom_reaction_" + str(i)] = ReactionData(
            reactants=[(1, spe)],
            products=[],
            kinetics=str(spe) + " * 1e-100",
        )


def _validate_species_names(
    meta_species_to_simulate: List_Species,
) -> set[str]:
    """Check all species have valid, unique names.

    Returns the set of names used.
    """
    names_used: set[str] = set()
    black_listed_names = {"Time", "Rev", "All"}

    for i, species in enumerate(meta_species_to_simulate):
        name = species.get_name()
        if DOT_SEPARATOR in name:
            raise CompilationError(
                f"In species: {name} \n _dot_ cannot be used in meta-species names"
            )
        if name in black_listed_names:
            raise CompilationError(
                f"The name {name} is not allowed for meta-species. Please change it."
            )
        if "$" in name:
            raise CompilationError(
                f"In species: {name} \n"
                "An error has occurred and one of"
                " the species was either not named"
                " or named with the"
                " restricted $ symbol"
            )
        if name in names_used:
            raise CompilationError(
                "Names must be unique for all species\n"
                + f"The repeated name is {name} "
                + f"in position {i}\n"
                + "Another possibility could be a"
                " repeated meta-species in the model"
            )
        names_used.add(name)

    return names_used


def _build_species_and_mappings(
    meta_species_to_simulate: List_Species,
    orthogonal_vector_structure: dict[str, Any],
) -> tuple[SpeciesForSbml, MappingsForSbml]:
    """Build species_for_sbml and mappings_for_sbml dicts."""
    species_for_sbml: SpeciesForSbml = {}
    mappings_for_sbml: MappingsForSbml = {}

    for spe_object in meta_species_to_simulate:
        species_string_list = ssg_construct_all_combinations(
            spe_object, STD_CHAR, orthogonal_vector_structure
        )
        for x in species_string_list:
            x[0] = x[0].get_name()
        mappings_for_sbml[spe_object.get_name()] = [
            ".".join(x) for x in species_string_list
        ]
        for species_string in species_string_list:
            sid = ConcreteSpeciesId(
                base=species_string[0],
                characteristics=tuple(species_string[1:]),
            )
            species_for_sbml[sid.to_sbml_id()] = 0

    return species_for_sbml, mappings_for_sbml


def _resolve_volume_and_dimension(
    volume: int | float | Quantity,
    dimension: int | None,
    species_counts: list[dict[str, Any]],
    model_context: ModelUnitContext | None = None,
) -> tuple[int | float, int, ParametersForSbml]:
    """Resolve volume and dimension, returning converted values."""
    if isinstance(volume, Quantity):
        dimension = uh_extract_length_dimension(str(volume.dimensionality), dimension)
    elif dimension is None:
        dimension = 3

    for count in species_counts:
        if isinstance(count["quantity"], Quantity) and uh_extract_length_dimension(
            str(count["quantity"].dimensionality), dimension
        ):
            dimension = uh_extract_length_dimension(
                str(count["quantity"].dimensionality), dimension
            )

    volume = uh_convert_volume(volume, dimension, model_context=model_context)

    if model_context is not None:
        vol_unit_id = model_context.get_sbml_volume_units_id()
    else:
        vol_unit_id = "dimensionless"
    parameters_for_sbml: ParametersForSbml = {"volume": (volume, vol_unit_id)}
    return volume, dimension, parameters_for_sbml


def _apply_count_to_species(
    quantity: Any,
    species_strings: set[str] | list[str],
    acc: CountAccumulator,
    ctx: CompilationContext,
) -> None:
    """Register a parameter if applicable, convert count, and assign."""
    if isinstance(quantity, mp_Mobspy_Parameter):
        acc.parameters_in_counts.add(quantity)
        if quantity.name in acc.parameters_used:
            acc.parameters_used[quantity.name].used_in = acc.parameters_used[
                quantity.name
            ].used_in.union(set(species_strings))
        else:
            acc.parameters_used[quantity.name] = ParameterUsedInfo(
                name=quantity.name,
                values=quantity.value,
                used_in=set(species_strings),
                object=quantity,
            )

    temp_count = uh_convert_counts(
        quantity, ctx.volume, ctx.dimension, model_context=ctx.model_context
    )
    for spe_str in species_strings:
        if isinstance(temp_count, float) and ctx.type_of_model != "deterministic":
            _logger.warning("The stochastic simulation rounds floats to integers")
            acc.species_for_sbml[spe_str] = int(temp_count)
        else:
            acc.species_for_sbml[spe_str] = temp_count
        acc.assigned_species.append(spe_str)


def _assign_initial_counts(
    species_counts: list[dict[str, Any]],
    species_for_sbml: SpeciesForSbml,
    orthogonal_vector_structure: dict[str, Any],
    parameters_used: ParametersUsed,
    ctx: CompilationContext,
) -> tuple[list[str], set[mp_Mobspy_Parameter]]:
    """Process species counts and assign initial values.

    Returns (assigned_species, parameters_in_counts).
    """
    acc = CountAccumulator(
        species_for_sbml=species_for_sbml,
        parameters_used=parameters_used,
    )

    # All-assignments first (lower priority, can be overridden)
    for count in species_counts:
        if ALL_CHAR not in count["characteristics"]:
            continue

        temp_set = set(count["characteristics"])
        temp_set.remove(ALL_CHAR)
        species_strings = [
            sid.to_sbml_id()
            for sid in ssg_construct_all_species_ids(
                count["object"],
                temp_set,
                orthogonal_vector_structure,
            )
        ]

        _apply_count_to_species(
            count["quantity"],
            species_strings,
            acc,
            ctx,
        )

    # Specific assignments (higher priority)
    for count in species_counts:
        if ALL_CHAR in count["characteristics"]:
            continue

        species_string = ssg_construct_species_id(
            count["object"],
            count["characteristics"],
            orthogonal_vector_structure,
        ).to_sbml_id()

        _apply_count_to_species(
            count["quantity"],
            [species_string],
            acc,
            ctx,
        )

    return acc.assigned_species, acc.parameters_in_counts


def _build_reactions(
    reactions_set: set[Any],
    meta_species_to_simulate: List_Species,
    orthogonal_vector_structure: dict[str, Any],
    ctx: CompilationContext,
) -> tuple[ReactionsForSbml, set[mp_Mobspy_Parameter]]:
    """Expand meta-reactions into concrete SBML reactions."""
    reactions_set = cof_create_all_not_reactions(reactions_set)

    parameters_in_reaction: set[mp_Mobspy_Parameter] = set()
    reactions_for_sbml, parameters_in_reaction = rc_create_all_reactions(
        reactions_set,
        meta_species_to_simulate,
        orthogonal_vector_structure,
        parameters_in_reaction,
        ctx,
    )
    return reactions_for_sbml, parameters_in_reaction


def _build_events(
    species_for_sbml: SpeciesForSbml,
    reactions_for_sbml: ReactionsForSbml,
    orthogonal_vector_structure: dict[str, Any],
    meta_species_to_simulate: List_Species,
    ctx: CompilationContext,
) -> tuple[dict[str, Any], set[mp_Mobspy_Parameter]]:
    """Build events and add phantom reactions for event-only species.

    Returns (events_for_sbml, parameters_in_events).
    """
    parameters_in_events: set[mp_Mobspy_Parameter] = set()
    events_for_sbml, species_in_events = eh_format_event_dictionary_for_sbml(
        species_for_sbml,
        ctx.event_dictionary or [],
        orthogonal_vector_structure,
        meta_species_to_simulate,
        parameters_in_events,
        ctx,
    )

    # Phantom reactions for species in events but not in reactions
    species_in_reactions: set[str] = set()
    for reaction in reactions_for_sbml.values():
        for reactant in reaction.reactants:
            species_in_reactions.add(reactant[1])
        for product in reaction.products:
            species_in_reactions.add(product[1])
    _add_phantom_reactions(
        reactions_for_sbml, species_in_events.difference(species_in_reactions)
    )

    # End condition event for continuous simulations
    if ctx.continuous_sim:
        end_event = EventData(
            trigger=ctx.ending_condition.generate_string(orthogonal_vector_structure),
            delay="0",
            assignments=[(END_FLAG_SPECIES_NAME, "1")],
        )
        reactions_for_sbml["phantom_reaction_end"] = ReactionData(
            reactants=[(10, END_FLAG_SPECIES_NAME)],
            products=[],
            kinetics=END_FLAG_SPECIES_NAME + " * 1e-100",
        )
        events_for_sbml["end_event"] = end_event

    return events_for_sbml, parameters_in_events


def _build_assignments(
    meta_species_to_simulate: List_Species,
    orthogonal_vector_structure: dict[str, Any],
) -> dict[str, Any]:
    """Compile species assignments for SBML."""
    non_processed_assignments: dict[str, Any] = {
        asgn: expression
        for spe in meta_species_to_simulate
        for asgn, expression in spe._assignments.items()
    }
    return asgi_Assign.compile_assignments_for_sbml(
        non_processed_assignments,
        orthogonal_vector_structure,
        meta_species_to_simulate,
    )


def _display_species(sbml_id: str) -> str:
    """Convert an SBML species ID to human-readable dotted form."""
    return ConcreteSpeciesId.from_sbml_id(sbml_id).to_display()


def _generate_model_string(  # noqa: PLR0913
    species_for_sbml: SpeciesForSbml,
    mappings_for_sbml: MappingsForSbml,
    parameters_for_sbml: ParametersForSbml,
    reactions_for_sbml: ReactionsForSbml,
    events_for_sbml: dict[str, Any],
    assignments_for_sbml: dict[str, Any],
    parameters_used: ParametersUsed | None = None,
) -> str:
    """Generate a human-readable model string for verbose output."""
    model_str = "\n"

    model_str += "Species\n"
    for spe in sorted(species_for_sbml):
        model_str += _display_species(spe) + "," + str(species_for_sbml[spe]) + "\n"

    model_str += "\nMappings\n"
    for map_key in sorted(mappings_for_sbml):
        model_str += map_key + " :\n"
        for element in sorted(mappings_for_sbml[map_key]):
            model_str += element + "\n"

    model_str += "\nParameters\n"
    for par in sorted(parameters_for_sbml):
        value = parameters_for_sbml[par][0]
        model_str += par + "," + str(value)
        if parameters_used and par in parameters_used:
            vals = parameters_used[par].values
            if isinstance(vals, list) and len(vals) > 1:
                model_str += f" (sweep: {vals})"
        model_str += "\n"

    model_str += "\nReactions\n"
    visible_reactions = {
        k: v for k, v in reactions_for_sbml.items() if "phantom" not in k
    }
    reaction_alpha = [
        str(x[1]).replace(DOT_SEPARATOR, ".")
        for x in sorted(visible_reactions.items(), key=lambda x: str(x[1]))
    ]
    for i, reac in enumerate(reaction_alpha):
        model_str += "reaction_" + str(i) + "," + reac + "\n"

    if events_for_sbml:
        model_str += "\nEvents\n"
        list_to_sort = sorted(str(events_for_sbml[key]) for key in events_for_sbml)
        for i, event_str in enumerate(list_to_sort):
            model_str += ("event_" + str(i) + "," + event_str + "\n").replace(
                DOT_SEPARATOR, "."
            )

    if assignments_for_sbml:
        model_str += "\nAssignments\n"
        list_to_sort = sorted(
            str(assignments_for_sbml[key]) for key in assignments_for_sbml
        )
        for i, asgn_str in enumerate(list_to_sort):
            model_str += ("assignment_" + str(i) + "," + asgn_str + "\n").replace(
                DOT_SEPARATOR, "."
            )

    return model_str


# ------------------------------------------------------------------
# Public entry point
# ------------------------------------------------------------------


def compile_model(  # noqa: PLR0913
    meta_species_to_simulate: List_Species,
    reactions_set: set[Any],
    species_counts: list[dict[str, Any]],
    orthogonal_vector_structure: dict[str, Any],
    *,
    volume: int | float | Quantity = 1,
    dimension: int | None = None,
    type_of_model: str = "deterministic",
    verbose: bool = True,
    event_dictionary: list[Any] | None = None,
    continuous_sim: bool = False,
    ending_condition: Any = None,
    skip_expression_check: bool = False,
    parameter_context: dict[str, mp_Mobspy_Parameter] | None = None,
    model_context: ModelUnitContext | None = None,
) -> CompilerResult:
    """Compile a MobsPy model into a backend-agnostic IR.

    Orchestrates a typed phased pipeline, building a
    :class:`~mobspy.types.ConcreteModel` internally and returning
    a :class:`~mobspy.types.CompilerResult` for backward compatibility.

    The ``ConcreteModel`` can be obtained from the result via
    ``result.to_concrete_model()``.

    Args:
        meta_species_to_simulate: Species in the model.
        reactions_set: Set of meta-reactions.
        species_counts: Initial count assignments.
        orthogonal_vector_structure: Characteristic-to-species mapping.
        volume: System volume.
        dimension: Spatial dimension (0-3).
        type_of_model: ``"deterministic"`` or ``"stochastic"``.
        verbose: Generate human-readable model string.
        event_dictionary: Packed event data from Simulation.
        continuous_sim: Whether this is a continuous simulation.
        ending_condition: End condition for continuous sims.
        skip_expression_check: Skip rate expression validation.
        parameter_context: Pre-built parameter registry.
        model_context: Unit context for SBML generation.

    Returns:
        CompilerResult with all compiled data. Use
        ``result.to_concrete_model()`` for the backend-agnostic IR.
    """
    # --- Phase 0: Freeze inheritance ---
    # All characteristics have been added by this point.
    # Eagerly compute ordered references for every species so
    # the lazy cache is never hit during compilation.
    meta_species_to_simulate = meta_species_to_simulate.remove_repeated_elements()
    for spe in meta_species_to_simulate:
        spe.freeze_references()

    # --- Phase 1: Species setup ---
    setup = phase_species_setup(meta_species_to_simulate, orthogonal_vector_structure)
    species_for_sbml = setup.species
    mappings_for_sbml = setup.mappings

    # --- Phase 2: Volume resolution ---
    vol_result = phase_volume_resolution(
        volume, dimension, species_counts, model_context
    )
    parameters_for_sbml: ParametersForSbml = {
        "volume": vol_result.volume_parameter,
    }

    # Build shared compilation context
    parameter_exist: dict[str, mp_Mobspy_Parameter]
    if parameter_context is None:
        parameter_exist = dict(mp_Mobspy_Parameter.parameter_stack)
    else:
        parameter_exist = parameter_context

    ctx = CompilationContext(
        volume=vol_result.volume,
        dimension=vol_result.dimension,
        type_of_model=type_of_model,
        model_context=model_context,
        parameter_exist=parameter_exist,
        skip_expression_check=skip_expression_check,
        continuous_sim=continuous_sim,
        ending_condition=ending_condition,
        event_dictionary=event_dictionary,
    )

    # Add end flag species for continuous simulations
    if continuous_sim:
        species_for_sbml[EndFlagSpecies.get_name()] = 0

    # Check for mole units in counts
    has_mole = any(
        isinstance(count["quantity"], Quantity)
        and "substance" in str(count["quantity"].dimensionality)
        for count in species_counts
    )

    # --- Phase 3: Count assignment ---
    parameters_used: ParametersUsed = {}
    count_result = phase_count_assignment(
        species_counts,
        species_for_sbml,
        orthogonal_vector_structure,
        parameters_used,
        ctx,
    )
    _add_to_parameters_to_sbml(
        parameters_used,
        parameters_for_sbml,
        count_result.parameters_in_counts,  # type: ignore[arg-type]
        model_context=model_context,
    )

    # --- Phase 4: Reaction expansion ---
    rxn_result = phase_reaction_expansion(
        reactions_set,
        meta_species_to_simulate,
        orthogonal_vector_structure,
        ctx,
    )
    reactions_for_sbml = rxn_result.reactions
    _add_to_parameters_to_sbml(
        parameters_used,
        parameters_for_sbml,
        rxn_result.parameters_in_reactions,  # type: ignore[arg-type]
        model_context=model_context,
    )

    # --- Phase 5: Duplicate detection (O(n)) ---
    phase_duplicate_detection(reactions_for_sbml)

    # --- Phase 6: Event building ---
    evt_result = phase_event_building(
        species_for_sbml,
        reactions_for_sbml,
        orthogonal_vector_structure,
        meta_species_to_simulate,
        ctx,
    )
    events_for_sbml = evt_result.events
    _add_to_parameters_to_sbml(
        parameters_used,
        parameters_for_sbml,
        evt_result.parameters_in_events,  # type: ignore[arg-type]
        model_context=model_context,
    )

    # --- Phase 7: Parameter validation ---
    all_params = (
        count_result.parameters_in_counts
        | rxn_result.parameters_in_reactions
        | evt_result.parameters_in_events
    )
    phase_parameter_validation(
        parameters_for_sbml,
        setup.names_used,
        all_params,  # type: ignore[arg-type]
    )

    # --- Phase 8: Assignments ---
    assignments_for_sbml = phase_assignments(
        meta_species_to_simulate,
        orthogonal_vector_structure,
    )

    # --- Build parameter object dict ---
    parameter_object_dict = {key: parameter_exist[key] for key in parameters_used}

    # --- Phase 9: Model string ---
    model_str = ""
    if verbose:
        model_str = _generate_model_string(
            species_for_sbml,
            mappings_for_sbml,
            parameters_for_sbml,
            reactions_for_sbml,
            events_for_sbml,
            assignments_for_sbml,
            parameters_used,
        )

    # --- Assemble ConcreteModel (the backend-agnostic IR) ---
    concrete = ConcreteModel(
        species=species_for_sbml,
        reactions=reactions_for_sbml,
        parameters=parameters_for_sbml,
        events=events_for_sbml,
        assignments=assignments_for_sbml,
        mappings=mappings_for_sbml,
        assigned_species=tuple(count_result.assigned_species),
        parameters_used=parameters_used,
        parameter_objects=parameter_object_dict,
        model_string=model_str,
        has_mole=has_mole,
        unit_context=model_context,
    )

    # Return CompilerResult for backward compatibility
    return concrete.to_compiler_result()
