from __future__ import annotations

from typing import TYPE_CHECKING, Any

from pint import Quantity

from mobspy.mobspy_logging import get_logger
from mobspy.modules.assignments_implementation import (
    Assign as asgi_Assign,
)
from mobspy.modules.compiler_operator_functions import (
    create_all_not_reactions as cof_create_all_not_reactions,
)
from mobspy.modules.event_functions import (
    format_event_dictionary_for_sbml as eh_format_event_dictionary_for_sbml,
)
from mobspy.modules.meta_class import EndFlagSpecies
from mobspy.modules.mobspy_parameters import (
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)
from mobspy.modules.reaction_construction_nb import (
    create_all_reactions as rc_create_all_reactions,
)
from mobspy.modules.species_string_generator import (
    construct_all_combinations as ssg_construct_all_combinations,
)
from mobspy.modules.species_string_generator import (
    construct_species_char_list as ssg_construct_species_char_list,
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
from mobspy.types import CompilerResult

if TYPE_CHECKING:
    from mobspy.types import (
        MappingsForSbml,
        ParametersForSbml,
        ParametersUsed,
        ReactionsForSbml,
        SpeciesForSbml,
    )

_logger = get_logger(__name__)


class Compiler:
    """Compiler for constructing and validating MobsPy models.

    Responsible for constructing the model from the
    meta-species and meta-reactions given and checking
    if it is a valid MobsPy model.
    """

    # Basic storage of defined variables for Compiler
    entity_counter: int = 0
    ref_characteristics_to_object: dict[str, Any] = {}
    last_rate: Any = None

    @classmethod
    def add_to_parameters_to_sbml(
        cls,
        parameters_used: ParametersUsed,
        parameters_for_sbml: ParametersForSbml,
        parameters_to_add: set[Any],
    ) -> None:
        for parameter in parameters_to_add:
            if parameter.name in parameters_used:
                parameters_used[parameter.name]["used_in"].add("$sbml")
            else:
                temp = {
                    "name": parameter.name,
                    "values": parameter.value,
                    "used_in": {"$sbml"},
                    "object": parameter,
                }
                parameters_used[parameter.name] = temp

            try:
                parameters_for_sbml[parameter.name] = (
                    parameter.value[0],
                    "dimensionless",
                )
            except Exception:
                parameters_for_sbml[parameter.name] = (
                    parameter.value,
                    "dimensionless",
                )

    @classmethod
    def override_get_item(cls, object_to_return: Any, item: Any) -> Any:
        """Store a rate before the reaction is defined.

        Due to priority in Python the item is stored
        before the reaction. So it is stored in the
        Compiler level and passed to the reaction object
        in the end. Important: the rate is used before
        the compilation level, it's used when the
        reaction has been completely defined.

        :param object_to_return: (Species or Reacting)
            returns the object __getitem__ was called on
        :param item: (int, float, callable, Quantity)
            stored reaction rate
        """
        cls.last_rate = item
        return object_to_return

    @classmethod
    def add_phantom_reactions(
        cls,
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
            reactions_for_sbml["phantom_reaction_" + str(i)] = {
                "re": [(1, spe)],
                "pr": [],
                "kin": str(spe) + " * 1e-100",
            }

    # ------------------------------------------------------------------
    # Private helpers: each handles one phase of compilation
    # ------------------------------------------------------------------

    @classmethod
    def _validate_species_names(
        cls,
        meta_species_to_simulate: Any,
    ) -> set[str]:
        """Check all species have valid, unique names.

        Returns the set of names used.
        """
        names_used: set[str] = set()
        black_listed_names = {"Time", "Rev", "All"}

        for i, species in enumerate(meta_species_to_simulate):
            name = species.get_name()
            if "_dot_" in name:
                _logger.error(
                    f"In species: {name} \n _dot_ cannot be used in meta-species names"
                )
            if name in black_listed_names:
                _logger.error(
                    f"The name {name} is not allowed for meta-species please change it"
                )
            if "$" in name:
                _logger.error(
                    f"In species: {name} \n"
                    "An error has occurred and one of"
                    " the species was either not named"
                    " or named with the"
                    " restricted $ symbol"
                )
            if name in names_used:
                _logger.error(
                    "Names must be unique for all species\n"
                    + f"The repeated name is {name} "
                    + f"in position {i}\n"
                    + "Another possibility could be a"
                    " repeated meta-species in the model"
                )
            names_used.add(name)

        return names_used

    @classmethod
    def _build_species_and_mappings(
        cls,
        meta_species_to_simulate: Any,
        orthogonal_vector_structure: dict[str, Any],
    ) -> tuple[SpeciesForSbml, MappingsForSbml]:
        """Build species_for_sbml and mappings_for_sbml dicts."""
        species_for_sbml: SpeciesForSbml = {}
        mappings_for_sbml: MappingsForSbml = {}

        for spe_object in meta_species_to_simulate:
            species_string_list = ssg_construct_all_combinations(
                spe_object, "std$", orthogonal_vector_structure
            )
            for x in species_string_list:
                x[0] = x[0].get_name()
            mappings_for_sbml[spe_object.get_name()] = [
                ".".join(x) for x in species_string_list
            ]
            for species_string in species_string_list:
                species_for_sbml["_dot_".join(species_string)] = 0

        return species_for_sbml, mappings_for_sbml

    @classmethod
    def _resolve_volume_and_dimension(
        cls,
        volume: int | float | Quantity,
        dimension: int | None,
        species_counts: list[dict[str, Any]],
    ) -> tuple[int | float, int, ParametersForSbml]:
        """Resolve volume and dimension, returning converted values."""
        if isinstance(volume, Quantity):
            dimension = uh_extract_length_dimension(
                str(volume.dimensionality), dimension
            )
        elif dimension is None:
            dimension = 3

        for count in species_counts:
            if isinstance(count["quantity"], Quantity):  # noqa: SIM102
                if uh_extract_length_dimension(
                    str(count["quantity"].dimensionality), dimension
                ):
                    dimension = uh_extract_length_dimension(
                        str(count["quantity"].dimensionality), dimension
                    )

        volume = uh_convert_volume(volume, dimension)
        parameters_for_sbml: ParametersForSbml = {"volume": (volume, "dimensionless")}
        return volume, dimension, parameters_for_sbml

    @classmethod
    def _assign_initial_counts(
        cls,
        species_counts: list[dict[str, Any]],
        species_for_sbml: SpeciesForSbml,
        orthogonal_vector_structure: dict[str, Any],
        volume: int | float,
        dimension: int,
        type_of_model: str,
        parameters_used: ParametersUsed,
    ) -> tuple[list[str], set[Any]]:
        """Process species counts and assign initial values.

        Returns (assigned_species, parameters_in_counts).
        """
        assigned_species: list[str] = []
        parameters_in_counts: set[Any] = set()

        # All-assignments first (lower priority, can be overridden)
        for count in species_counts:
            if "all$" not in count["characteristics"]:
                continue

            temp_set = set(count["characteristics"])
            temp_set.remove("all$")
            species_strings = ssg_construct_all_combinations(
                count["object"], temp_set, orthogonal_vector_structure, symbol="_dot_"
            )

            if isinstance(count["quantity"], mp_Mobspy_Parameter):
                parameters_in_counts.add(count["quantity"])
                if count["quantity"].name in parameters_used:
                    parameters_used[count["quantity"].name]["used_in"] = (
                        parameters_used[count["quantity"].name]["used_in"].union(
                            set(species_strings)
                        )
                    )
                else:
                    temp = {
                        "name": count["quantity"].name,
                        "values": count["quantity"].value,
                        "used_in": set(species_strings),
                        "object": count["quantity"],
                    }
                    parameters_used[count["quantity"].name] = temp

            temp_count = uh_convert_counts(count["quantity"], volume, dimension)
            for spe_str in species_strings:
                if type(temp_count) == float and not type_of_model == "deterministic":  # noqa: SIM201, E721
                    _logger.warning(
                        "The stochastic simulation rounds floats to integers"
                    )
                    species_for_sbml[spe_str] = int(temp_count)
                else:
                    species_for_sbml[spe_str] = temp_count
                assigned_species.append(spe_str)

        # Specific assignments (higher priority)
        for count in species_counts:
            if "all$" in count["characteristics"]:
                continue

            species_string = ssg_construct_species_char_list(
                count["object"],
                count["characteristics"],
                orthogonal_vector_structure,
                symbol="_dot_",
            )

            if isinstance(count["quantity"], mp_Mobspy_Parameter):
                parameters_in_counts.add(count["quantity"])
                if count["quantity"].name in parameters_used:
                    parameters_used[count["quantity"].name]["used_in"].add(
                        species_string
                    )
                else:
                    temp = {
                        "name": count["quantity"].name,
                        "values": count["quantity"].value,
                        "used_in": {species_string},
                        "object": count["quantity"],
                    }
                    parameters_used[count["quantity"].name] = temp

            temp_count = uh_convert_counts(count["quantity"], volume, dimension)
            if type(temp_count) == float and not type_of_model == "deterministic":  # noqa: SIM201, E721
                _logger.warning("The stochastic simulation rounds floats to integers")
                species_for_sbml[species_string] = int(temp_count)
            else:
                species_for_sbml[species_string] = temp_count
            assigned_species.append(species_string)

        return assigned_species, parameters_in_counts

    @classmethod
    def _build_reactions(
        cls,
        reactions_set: set[Any],
        meta_species_to_simulate: Any,
        orthogonal_vector_structure: dict[str, Any],
        type_of_model: str,
        dimension: int,
        parameter_exist: dict[str, Any],
        skip_expression_check: bool,
    ) -> tuple[ReactionsForSbml, set[Any]]:
        """Expand meta-reactions into concrete SBML reactions."""
        reactions_set = cof_create_all_not_reactions(reactions_set)

        parameters_in_reaction: set[Any] = set()
        reactions_for_sbml, parameters_in_reaction = rc_create_all_reactions(
            reactions_set,
            meta_species_to_simulate,
            orthogonal_vector_structure,
            type_of_model,
            dimension,
            parameter_exist,
            parameters_in_reaction,
            skip_expression_check,
        )
        return reactions_for_sbml, parameters_in_reaction

    @classmethod
    def _check_duplicate_reactions(
        cls,
        reactions_for_sbml: ReactionsForSbml,
    ) -> None:
        """Warn about duplicate reactions (O(n^2) check)."""
        keys = list(reactions_for_sbml)
        for i, r1 in enumerate(keys):
            for j, r2 in enumerate(keys):
                if i >= j:
                    continue
                if (
                    reactions_for_sbml[r1]["re"] == reactions_for_sbml[r2]["re"]
                    and reactions_for_sbml[r1]["pr"] == reactions_for_sbml[r2]["pr"]
                    and reactions_for_sbml[r1]["kin"] == reactions_for_sbml[r2]["kin"]
                ):
                    _logger.warning(
                        "The following reaction: \n"
                        + f"{reactions_for_sbml[r1]} \n"
                        + "Is doubled. Was that intentional? \n"
                    )

    @classmethod
    def _build_events(
        cls,
        species_for_sbml: SpeciesForSbml,
        reactions_for_sbml: ReactionsForSbml,
        event_dictionary: list[Any] | None,
        orthogonal_vector_structure: dict[str, Any],
        volume: int | float,
        dimension: int,
        meta_species_to_simulate: Any,
        parameter_exist: dict[str, Any],
        continuous_sim: bool,
        ending_condition: Any,
    ) -> tuple[dict[str, Any], set[Any]]:
        """Build events and add phantom reactions for event-only species.

        Returns (events_for_sbml, parameters_in_events).
        """
        parameters_in_events: set[Any] = set()
        events_for_sbml, species_in_events = eh_format_event_dictionary_for_sbml(
            species_for_sbml,
            event_dictionary,
            orthogonal_vector_structure,
            volume,
            dimension,
            meta_species_to_simulate,
            parameter_exist,
            parameters_in_events,
        )

        # Phantom reactions for species in events but not in reactions
        species_in_reactions: set[str] = set()
        for reaction in reactions_for_sbml.values():
            for reactant in reaction["re"]:
                species_in_reactions.add(reactant[1])
            for product in reaction["pr"]:
                species_in_reactions.add(product[1])
        cls.add_phantom_reactions(
            reactions_for_sbml, species_in_events.difference(species_in_reactions)
        )

        # End condition event for continuous simulations
        if continuous_sim:
            end_event = {
                "trigger": ending_condition.generate_string(
                    orthogonal_vector_structure
                ),
                "delay": "0",
                "assignments": [("_End_Flag_MetaSpecies", "1")],
            }
            reactions_for_sbml["phantom_reaction_end"] = {
                "re": [(10, "_End_Flag_MetaSpecies")],
                "pr": [],
                "kin": "_End_Flag_MetaSpecies * 1e-100",
            }
            events_for_sbml["end_event"] = end_event

        return events_for_sbml, parameters_in_events

    @classmethod
    def _validate_parameters(
        cls,
        parameters_for_sbml: ParametersForSbml,
        names_used: set[str],
        parameters_in_counts: set[Any],
        parameters_in_reaction: set[Any],
        parameters_in_events: set[Any],
    ) -> None:
        """Check parameter names are unique and don't collide with species."""
        for p in parameters_for_sbml:
            if p in names_used:
                _logger.error(
                    "Parameters names must be unique"
                    " and they must not share a"
                    " name with a species"
                )
            names_used.add(p)

        all_params = (
            set()
            .union(parameters_in_counts)
            .union(parameters_in_reaction)
            .union(parameters_in_events)
        )
        for p1 in all_params:
            for p2 in all_params:
                if p1 is not p2 and p1.name == p2.name:
                    _logger.error(
                        "There are two different Parameter Objects with the same name"
                    )

    @classmethod
    def _build_assignments(
        cls,
        meta_species_to_simulate: Any,
        orthogonal_vector_structure: dict[str, Any],
    ) -> dict[str, Any]:
        """Compile species assignments for SBML."""
        non_processed_assignments: dict[str, Any] = {}
        for spe in meta_species_to_simulate:
            for asgn, expression in spe._assignments.items():
                non_processed_assignments[asgn] = expression
        return asgi_Assign.compile_assignments_for_sbml(
            non_processed_assignments,
            orthogonal_vector_structure,
            meta_species_to_simulate,
        )

    @classmethod
    def _generate_model_string(
        cls,
        species_for_sbml: SpeciesForSbml,
        mappings_for_sbml: MappingsForSbml,
        parameters_for_sbml: ParametersForSbml,
        reactions_for_sbml: ReactionsForSbml,
        events_for_sbml: dict[str, Any],
        assignments_for_sbml: dict[str, Any],
    ) -> str:
        """Generate a human-readable model string for verbose output."""
        model_str = "\n"

        model_str += "Species\n"
        for spe in sorted(species_for_sbml):
            model_str += (
                spe.replace("_dot_", ".") + "," + str(species_for_sbml[spe]) + "\n"
            )

        model_str += "\nMappings\n"
        for map_key in sorted(mappings_for_sbml):
            model_str += map_key + " :\n"
            for element in sorted(mappings_for_sbml[map_key]):
                model_str += element + "\n"

        model_str += "\nParameters\n"
        for par in sorted(parameters_for_sbml):
            model_str += par + "," + str(parameters_for_sbml[par][0]) + "\n"

        model_str += "\nReactions\n"
        visible_reactions = {
            k: v for k, v in reactions_for_sbml.items() if "phantom" not in k
        }
        reaction_alpha = [
            str(x[1]).replace("_dot_", ".")
            for x in sorted(visible_reactions.items(), key=lambda x: str(x[1]))
        ]
        for i, reac in enumerate(reaction_alpha):
            model_str += "reaction_" + str(i) + "," + reac + "\n"

        if events_for_sbml:
            model_str += "\nEvents\n"
            list_to_sort = sorted(str(events_for_sbml[key]) for key in events_for_sbml)
            for i, event_str in enumerate(list_to_sort):
                model_str += ("event_" + str(i) + "," + event_str + "\n").replace(
                    "_dot_", "."
                )

        if assignments_for_sbml:
            model_str += "\nAssignments\n"
            list_to_sort = sorted(
                str(assignments_for_sbml[key]) for key in assignments_for_sbml
            )
            for i, asgn_str in enumerate(list_to_sort):
                model_str += ("assignment_" + str(i) + "," + asgn_str + "\n").replace(
                    "_dot_", "."
                )

        return model_str

    # ------------------------------------------------------------------
    # Main compile entry point
    # ------------------------------------------------------------------

    @classmethod
    def compile(
        cls,
        meta_species_to_simulate: Any,
        reactions_set: set[Any],
        species_counts: list[dict[str, Any]],
        orthogonal_vector_structure: dict[str, Any],
        volume: int | float | Quantity = 1,
        dimension: int | None = None,
        type_of_model: str = "deterministic",
        verbose: bool = True,
        event_dictionary: list[Any] | None = None,
        continuous_sim: bool = False,
        ending_condition: Any = None,
        skip_expression_check: bool = False,
    ) -> CompilerResult:
        """Compile a MobsPy model into SBML-ready data structures.

        Orchestrates validation, species expansion, count assignment,
        reaction construction, event building, and model string generation.
        """
        # Phase 1: Validate species names
        meta_species_to_simulate = meta_species_to_simulate.remove_repeated_elements()
        names_used = cls._validate_species_names(meta_species_to_simulate)

        # Phase 2: Initialize parameters
        parameters_used: ParametersUsed = {}
        parameter_exist = {}
        if mp_Mobspy_Parameter.parameter_stack != {}:
            parameter_exist = mp_Mobspy_Parameter.parameter_stack

        for species in meta_species_to_simulate:
            species.order_references()

        # Phase 3: Build species and mappings
        species_for_sbml, mappings_for_sbml = cls._build_species_and_mappings(
            meta_species_to_simulate, orthogonal_vector_structure
        )

        # Phase 4: Resolve volume and dimension
        volume, dimension, parameters_for_sbml = cls._resolve_volume_and_dimension(
            volume, dimension, species_counts
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

        # Phase 5: Assign initial counts
        assigned_species, parameters_in_counts = cls._assign_initial_counts(
            species_counts,
            species_for_sbml,
            orthogonal_vector_structure,
            volume,
            dimension,
            type_of_model,
            parameters_used,
        )
        cls.add_to_parameters_to_sbml(
            parameters_used, parameters_for_sbml, parameters_in_counts
        )

        # Phase 6: Build reactions
        reactions_for_sbml, parameters_in_reaction = cls._build_reactions(
            reactions_set,
            meta_species_to_simulate,
            orthogonal_vector_structure,
            type_of_model,
            dimension,
            parameter_exist,
            skip_expression_check,
        )
        cls.add_to_parameters_to_sbml(
            parameters_used, parameters_for_sbml, parameters_in_reaction
        )

        # Phase 7: Check for duplicate reactions
        cls._check_duplicate_reactions(reactions_for_sbml)

        # Phase 8: Build events
        events_for_sbml, parameters_in_events = cls._build_events(
            species_for_sbml,
            reactions_for_sbml,
            event_dictionary,
            orthogonal_vector_structure,
            volume,
            dimension,
            meta_species_to_simulate,
            parameter_exist,
            continuous_sim,
            ending_condition,
        )
        cls.add_to_parameters_to_sbml(
            parameters_used, parameters_for_sbml, parameters_in_events
        )

        # Phase 9: Validate parameters
        cls._validate_parameters(
            parameters_for_sbml,
            names_used,
            parameters_in_counts,
            parameters_in_reaction,
            parameters_in_events,
        )

        # Phase 10: Store parameter objects
        parameter_object_dict = {
            key: mp_Mobspy_Parameter.parameter_stack[key] for key in parameters_used
        }

        # Phase 11: Build assignments
        assignments_for_sbml = cls._build_assignments(
            meta_species_to_simulate, orthogonal_vector_structure
        )

        # Phase 12: Generate model string
        model_str = ""
        if verbose:
            model_str = cls._generate_model_string(
                species_for_sbml,
                mappings_for_sbml,
                parameters_for_sbml,
                reactions_for_sbml,
                events_for_sbml,
                assignments_for_sbml,
            )

        return CompilerResult(
            species_for_sbml=species_for_sbml,
            reactions_for_sbml=reactions_for_sbml,
            parameters_for_sbml=parameters_for_sbml,
            mappings_for_sbml=mappings_for_sbml,
            model_str=model_str,
            events_for_sbml=events_for_sbml,
            assigned_species=assigned_species,
            parameters_used=parameters_used,
            parameter_object_dict=parameter_object_dict,
            assignments_for_sbml=assignments_for_sbml,
            has_mole=has_mole,
        )
