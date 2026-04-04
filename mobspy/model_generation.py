"""
Model generation module: SBML and Antimony export from compiled MobsPy models.

Provides the ModelGenerationMixin class that adds generate_sbml(), generate_antimony(),
and compose_sbml() capabilities to the Simulation class.
"""

from __future__ import annotations

from random import randint as rd_randint
from typing import TYPE_CHECKING, Any

from mobspy.exceptions import SBMLError
from mobspy.mobspy_logging import get_logger
from mobspy.sbml_simulator.builder import build as sbml_build
from mobspy.types import EventData, ReactionData, SBMLModelData

if TYPE_CHECKING:
    from mobspy.types import (
        CompiledModelDict,
        ParameterSweepList,
        SBMLModelDict,
        SimulationParameters,
    )

_logger = get_logger(__name__)


def _antimony_species_block(sbml_data: SBMLModelDict) -> str:
    """Format species declarations for Antimony."""
    result = ""
    if sbml_data.species_for_sbml:
        for species_name, species_count in sbml_data.species_for_sbml.items():
            result += f"    {species_name} = {species_count} dimensionless\n"
    return result


def _antimony_parameters_block(sbml_data: SBMLModelDict) -> str:
    """Format parameter declarations for Antimony."""
    result = ""
    if sbml_data.parameters_for_sbml:
        for parameter_name, parameter_value in sbml_data.parameters_for_sbml.items():
            if parameter_name == "volume":
                continue
            result += f"    {parameter_name} = {parameter_value[0]} dimensionless\n"
    return result


def _antimony_assignments_block(sbml_data: SBMLModelDict) -> str:
    """Format assignment rules for Antimony."""
    result = ""
    if sbml_data.assignments_for_sbml:
        for assign_data in sbml_data.assignments_for_sbml.values():
            result += f"    {assign_data.species} := {assign_data.expression}\n"
    return result


def _antimony_reactions_block(sbml_data: SBMLModelDict) -> str:
    """Format reactions for Antimony."""
    result = ""
    for reaction_name, reaction_data in sbml_data.reactions_for_sbml.items():
        if "phantom" in reaction_name:
            continue

        result += f"    {reaction_name}: "

        for i, r in enumerate(reaction_data.reactants):
            if i == 0 and r[0] > 1:
                result += f"{r[0]}*{r[1]}"
            elif i == 0:
                result += f"{r[1]}"
            elif r[0] > 1:
                result += f" + {r[0]} {r[1]}"
            else:
                result += f" + {r[1]}"

        result += " -> "

        for i, p in enumerate(reaction_data.products):
            if i == 0 and p[0] > 1:
                result += f" {p[0]} {p[1]}"
            elif i == 0:
                result += f" {p[1]}"
            elif p[0] > 1:
                result += f" + {p[0]}*{p[1]}"
            else:
                result += f" + {p[1]}"

        result += f"; {reaction_data.kinetics}\n"
    return result


def _antimony_events_block(sbml_data: SBMLModelDict) -> str:
    """Format events for Antimony."""
    result = ""
    if sbml_data.events_for_sbml:
        for event_name, event_data in sbml_data.events_for_sbml.items():
            if event_data.trigger == "true":
                result += f"    {event_name}: at(time > {event_data.delay}): "
            else:
                result += f"    {event_name}: at({event_data.trigger}): "

            for asg in event_data.assignments:
                result += f" {asg[0]}={asg[1]},"
            result = result[:-1] + "\n"
    return result


def _compose_reactions(
    new_sbml_file: SBMLModelDict,
    i: int,
    flag_species_name: str,
    current_sbml_reaction: dict[str, ReactionData],
) -> None:
    """Gate each reaction's kinetics by the flag species for simulation *i*."""
    for reaction_key, reaction in current_sbml_reaction.items():
        if "phantom" in reaction_key:
            continue
        new_reaction = ReactionData(
            reactants=reaction.reactants,
            products=reaction.products,
            kinetics="("
            + reaction.kinetics.replace("volume", f"_vol{i}")
            + ") * "
            + str(flag_species_name),
        )
        reaction_number = len(new_sbml_file.reactions_for_sbml)
        new_sbml_file.reactions_for_sbml["reaction_" + str(reaction_number)] = (
            new_reaction
        )


def _compose_parameters(
    new_sbml_file: SBMLModelDict,
    sim_index: int,
    sim_sbml: CompiledModelDict,
) -> None:
    """Copy parameters into the composite model.

    Renames volume per simulation index.
    """
    for par in sim_sbml.parameters_for_sbml:
        if par == "volume":
            new_sbml_file.parameters_for_sbml["_vol" + str(sim_index)] = (
                sim_sbml.parameters_for_sbml[par]
            )
        else:
            new_sbml_file.parameters_for_sbml[par] = sim_sbml.parameters_for_sbml[par]


def generate_sbml_strings(
    sbml_data_list: list[list[Any]],
) -> list[str]:
    """Generate SBML XML strings from compiled model data.

    Standalone function that does not require a Simulation instance.

    Args:
        sbml_data_list: Nested list of SBMLModelData (parameter sweeps).

    Returns:
        List of SBML XML strings.
    """
    results: list[str] = []
    for parameter_sweep in sbml_data_list:
        for sbml_data in parameter_sweep:
            model_ctx = getattr(sbml_data, "model_context", None)
            results.append(
                sbml_build(
                    SBMLModelData(
                        species_for_sbml=sbml_data.species_for_sbml,
                        parameters_for_sbml=sbml_data.parameters_for_sbml,
                        reactions_for_sbml=sbml_data.reactions_for_sbml,
                        events_for_sbml=sbml_data.events_for_sbml,
                        assignments_for_sbml=sbml_data.assignments_for_sbml,
                    ),
                    model_context=model_ctx,
                )
            )
    return results


def generate_antimony_strings(
    sbml_data_list: list[list[Any]],
    model_name: str | None = None,
) -> list[str]:
    """Generate Antimony model strings from compiled model data.

    Standalone function that does not require a Simulation instance.

    Args:
        sbml_data_list: Nested list of SBMLModelData (parameter sweeps).
        model_name: Optional model name (random if None).

    Returns:
        List of Antimony model strings.
    """
    results: list[str] = []
    for parameter_sweep in sbml_data_list:
        for sbml_data in parameter_sweep:
            if model_name is None:
                antimony_model = f"model mobspy_{rd_randint(0, 100000)} \n"  # noqa: S311
            else:
                antimony_model = f"model {model_name} \n"

            antimony_model += _antimony_species_block(sbml_data)
            antimony_model += _antimony_parameters_block(sbml_data)
            antimony_model += _antimony_assignments_block(sbml_data)
            antimony_model += _antimony_reactions_block(sbml_data)
            antimony_model += _antimony_events_block(sbml_data)
            antimony_model += "end"
            results.append(antimony_model)
    return results


class ModelGenerationMixin:
    """Mixin providing SBML and Antimony model generation for Simulation."""

    # Attributes provided by Simulation
    _list_of_parameters: list[SimulationParameters]
    sbml_data_list: ParameterSweepList
    _species_for_sbml: dict[str, Any] | None

    def compile(self, verbose: bool = True) -> Any:
        """Provided by Simulation at runtime."""
        ...

    def _assemble_multi_simulation_structure(self) -> None:
        """Provided by Simulation at runtime."""
        ...

    def compose_sbml(self) -> list[list[SBMLModelDict]]:
        """Merge concatenated simulations into single SBML models using flag species."""
        self._check_compose_convertible()

        return [
            [self._compose_single_chain(multi_sims)]
            for multi_sims in self.sbml_data_list
        ]

    def _check_compose_convertible(self) -> None:
        """Raise if the simulation chain cannot be composed into a single SBML."""
        if len(self._list_of_parameters) == 1:
            raise SBMLError(
                "Single simulations cannot generate a composed sbml or antimony string"
            )
        for i in range(len(self._list_of_parameters)):
            if self._list_of_parameters[i]["_end_condition"] is not None:
                raise SBMLError(
                    "Composite Simulations with conditional "
                    "duration cannot be converted to "
                    "sbml or antimony"
                )

    def _compose_single_chain(
        self, multi_sims: list[CompiledModelDict]
    ) -> SBMLModelDict:
        """Compose one chain of simulations into a single SBML model."""
        new_sbml_file: SBMLModelDict = SBMLModelData()
        initial_sim = multi_sims[0]

        new_sbml_file.species_for_sbml = initial_sim.species_for_sbml
        _compose_parameters(new_sbml_file, 0, initial_sim)
        new_sbml_file.assignments_for_sbml = initial_sim.assignments_for_sbml
        _compose_reactions(new_sbml_file, 0, "_SFS_0", initial_sim.reactions_for_sbml)
        self._compose_event(new_sbml_file, "_SFS_1", 0, 0, initial_sim)

        self._compose_subsequent_sims(new_sbml_file, multi_sims)

        for i in range(len(multi_sims)):
            new_sbml_file.species_for_sbml["_SFS_" + str(i)] = 0
        new_sbml_file.species_for_sbml["_SFS_0"] = 1

        return new_sbml_file

    def _compose_event(
        self,
        new_sbml_file: SBMLModelDict,
        next_spe: str,
        simulation_index: int,
        cul_duration: float | int,
        sim_sbml: CompiledModelDict,
    ) -> None:
        """Create a time- or condition-triggered event for the next simulation phase."""
        if self._list_of_parameters[simulation_index]["_end_condition"] is None:
            event_name = "e" + str(len(new_sbml_file.events_for_sbml))
            event = EventData(
                trigger="true",
                delay=cul_duration,
                assignments=[(next_spe, 1)],
            )
            new_sbml_file.events_for_sbml[event_name] = event
        else:
            end_trigger = sim_sbml.events_for_sbml["end_event"].trigger
            event_name = "e" + str(len(new_sbml_file.events_for_sbml))
            event = EventData(
                trigger=end_trigger,
                delay=0,
                assignments=[(next_spe, 1)],
            )
            new_sbml_file.events_for_sbml[event_name] = event

    def _compose_a_sim(  # noqa: PLR0913
        self,
        new_sbml_file: SBMLModelDict,
        i: int,
        sim_sbml: CompiledModelDict,
        pre_spe: str,
        next_spe: str,
        cul_duration: float | int,
        skip_end_event: bool,
    ) -> None:
        """Integrate one simulation's components into the composite model."""
        _compose_parameters(new_sbml_file, i, sim_sbml)
        _compose_reactions(new_sbml_file, i, pre_spe, sim_sbml.reactions_for_sbml)
        if not skip_end_event:
            self._compose_event(new_sbml_file, next_spe, i, cul_duration, sim_sbml)

        for spe in sim_sbml.species_for_sbml:
            if spe not in new_sbml_file.species_for_sbml and spe[0] != "_":
                event_number = len(new_sbml_file.events_for_sbml)
                spe_event = EventData(
                    trigger=f"_SFS_{i!s} > 0",
                    delay=0,
                    assignments=[(spe, sim_sbml.species_for_sbml[spe])],
                )
                new_sbml_file.events_for_sbml["e" + str(event_number)] = spe_event
                new_sbml_file.species_for_sbml[spe] = 0

    def _compose_subsequent_sims(
        self,
        new_sbml_file: SBMLModelDict,
        multi_sims: list[CompiledModelDict],
    ) -> None:
        """Iterate over all simulations after the first and compose them."""
        cul_duration: float | int = 0
        skip_end_event = False
        for i, sim_sbml in enumerate(multi_sims):
            if i == 0:
                cul_duration = self._list_of_parameters[i]["duration"]
                continue
            if i == len(multi_sims) - 1:
                skip_end_event = True
            else:
                cul_duration = cul_duration + self._list_of_parameters[i]["duration"]

            if sim_sbml.assignments_for_sbml != {}:
                _logger.warning("Assignments beyond the initial simulation are ignored")

            pre_spe = "_SFS_" + str(i)
            next_spe = "_SFS_" + str(i + 1)
            self._compose_a_sim(
                new_sbml_file,
                i,
                sim_sbml,
                pre_spe,
                next_spe,
                cul_duration,
                skip_end_event,
            )

    def parse_volume_name_for_antimony(self) -> list[list[SBMLModelDict]]:
        """Rename 'volume' to '_vol' in reaction kinetics for Antimony compatibility."""
        new_sims = []
        for multi_sims in self.sbml_data_list:
            sim_sbml = multi_sims[0]
            new_sbml_file: SBMLModelDict = SBMLModelData(
                species_for_sbml=sim_sbml.species_for_sbml,
                parameters_for_sbml=sim_sbml.parameters_for_sbml,
                events_for_sbml=sim_sbml.events_for_sbml,
                assignments_for_sbml=sim_sbml.assignments_for_sbml,
            )

            new_sbml_file.parameters_for_sbml["_vol"] = sim_sbml.parameters_for_sbml[
                "volume"
            ]

            for re_name, reaction in sim_sbml.reactions_for_sbml.items():
                new_reaction = ReactionData(
                    reactants=reaction.reactants,
                    products=reaction.products,
                    kinetics=reaction.kinetics.replace("volume", "_vol"),
                )
                new_sbml_file.reactions_for_sbml[re_name] = new_reaction

            new_sims.append([new_sbml_file])
        return new_sims

    def generate_sbml(self, compose: bool = False) -> list[str]:
        """Generate SBML strings from the current stored models.

        Args:
            compose: Join composite simulations into a single SBML.

        Returns:
            List of SBML XML strings.
        """
        if self._species_for_sbml is None:
            self.compile(verbose=False)
        self._assemble_multi_simulation_structure()

        data = self.compose_sbml() if compose else self.sbml_data_list
        return generate_sbml_strings(data)

    def generate_antimony(
        self, compose: bool = False, model_name: str | None = None
    ) -> list[str]:
        """Generate Antimony model strings from the current stored models.

        Args:
            compose: Join composite simulations into a single SBML.
            model_name: Desired name of the model.
        """
        if self._species_for_sbml is None:
            self.compile(verbose=False)
        self._assemble_multi_simulation_structure()

        data = self.compose_sbml() if compose else self.parse_volume_name_for_antimony()
        return generate_antimony_strings(data, model_name=model_name)
