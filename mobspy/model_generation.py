"""
Model generation module: SBML and Antimony export from compiled MobsPy models.

Provides the ModelGenerationMixin class that adds generate_sbml(), generate_antimony(),
and compose_sbml() capabilities to the Simulation class.
"""

from __future__ import annotations

from random import randint as rd_randint
from typing import TYPE_CHECKING, Any

from mobspy.mobspy_logging import get_logger
from mobspy.sbml_simulator.builder import build as sbml_build

if TYPE_CHECKING:
    from mobspy.types import (
        CompiledModelDict,
        EventData,
        ParameterSweepList,
        ReactionData,
        SBMLModelDict,
        SimulationParameters,
    )

logger = get_logger(__name__)


class ModelGenerationMixin:
    """Mixin providing SBML and Antimony model generation for Simulation."""

    # Attributes provided by Simulation
    _list_of_parameters: list[SimulationParameters]
    sbml_data_list: ParameterSweepList
    _species_for_sbml: dict[str, Any] | None

    def compose_sbml(self) -> list[list[SBMLModelDict]]:
        list_of_composite_dicts_for_sbml = []

        def check_convertible() -> None:
            if len(self._list_of_parameters) == 1:
                logger.error(
                    "Single simulations cannot generate a composed "
                    "sbml or antimony string"
                )

            for i in range(len(self._list_of_parameters)):
                if self._list_of_parameters[i]["_end_condition"] is not None:
                    logger.error(
                        "Composite Simulations with conditional "
                        "duration cannot be converted to "
                        "sbml or antimony"
                    )

        check_convertible()

        def reaction_process(
            i: int,
            flag_species_name: str,
            current_sbml_reaction: dict[str, ReactionData],
        ) -> None:
            for reaction_key, reaction in current_sbml_reaction.items():
                if "phantom" in reaction_key:
                    continue
                new_reaction: ReactionData = {
                    "re": reaction["re"],
                    "pr": reaction["pr"],
                    "kin": "("
                    + reaction["kin"].replace("volume", f"_vol{i}")
                    + ") * "
                    + str(flag_species_name),
                }

                reaction_number = len(new_sbml_file["reactions_for_sbml"])
                new_sbml_file["reactions_for_sbml"][
                    "reaction_" + str(reaction_number)
                ] = new_reaction

        def event_process(
            next_spe: str,
            simulation_index: int,
            cul_duration: float | int,
            sim_sbml: CompiledModelDict,
        ) -> None:
            if self._list_of_parameters[simulation_index]["_end_condition"] is None:
                event_name = "e" + str(len(new_sbml_file["events_for_sbml"]))
                event: EventData = {
                    "trigger": "true",
                    "delay": cul_duration,
                    "assignments": [(next_spe, 1)],
                }
                new_sbml_file["events_for_sbml"][event_name] = event
            else:
                end_trigger = sim_sbml["events_for_sbml"]["end_event"]["trigger"]
                event_name = "e" + str(len(new_sbml_file["events_for_sbml"]))
                event = {
                    "trigger": end_trigger,
                    "delay": 0,
                    "assignments": [(next_spe, 1)],
                }
                new_sbml_file["events_for_sbml"][event_name] = event

        def parameter_process(sim_index: int, sim_sbml: CompiledModelDict) -> None:
            for par in sim_sbml["parameters_for_sbml"]:
                if par == "volume":
                    new_sbml_file["parameters_for_sbml"]["_vol" + str(sim_index)] = (
                        sim_sbml["parameters_for_sbml"][par]
                    )
                else:
                    new_sbml_file["parameters_for_sbml"][par] = sim_sbml[
                        "parameters_for_sbml"
                    ][par]

        def process_a_sim(
            i: int,
            sim_sbml: CompiledModelDict,
            pre_spe: str,
            next_spe: str,
            cul_duration: float | int,
            skip_end_event: bool,
        ) -> None:
            parameter_process(i, sim_sbml)
            reaction_process(i, pre_spe, sim_sbml["reactions_for_sbml"])
            if not skip_end_event:
                event_process(next_spe, i, cul_duration, sim_sbml)

            for spe in sim_sbml["species_for_sbml"]:
                if spe not in new_sbml_file["species_for_sbml"] and spe[0] != "_":
                    event_number = len(new_sbml_file["events_for_sbml"])
                    spe_event: EventData = {
                        "trigger": f"_SFS_{str(i)} > 0",
                        "delay": 0,
                        "assignments": [(spe, sim_sbml["species_for_sbml"][spe])],
                    }
                    new_sbml_file["events_for_sbml"]["e" + str(event_number)] = (
                        spe_event
                    )
                    new_sbml_file["species_for_sbml"][spe] = 0

        def process_simulations(
            multi_sims: list[CompiledModelDict],
        ) -> None:
            cul_duration: float | int = 0
            skip_end_event = False
            for i, sim_sbml in enumerate(multi_sims):
                if i == 0:
                    cul_duration = self._list_of_parameters[i]["duration"]
                    continue
                elif i == len(multi_sims) - 1:
                    skip_end_event = True
                else:
                    cul_duration = (
                        cul_duration + self._list_of_parameters[i]["duration"]
                    )

                if sim_sbml["assignments_for_sbml"] != {}:
                    logger.warning(
                        "Assignments beyond the initial simulation are ignored"
                    )

                pre_spe = "_SFS_" + str(i)
                next_spe = "_SFS_" + str(i + 1)

                process_a_sim(
                    i, sim_sbml, pre_spe, next_spe, cul_duration, skip_end_event
                )

        for multi_sims in self.sbml_data_list:
            new_sbml_file: SBMLModelDict = {
                "species_for_sbml": {},
                "parameters_for_sbml": {},
                "reactions_for_sbml": {},
                "events_for_sbml": {},
                "assignments_for_sbml": {},
            }

            initial_sim = multi_sims[0]

            new_sbml_file["species_for_sbml"] = initial_sim["species_for_sbml"]
            parameter_process(0, initial_sim)
            new_sbml_file["assignments_for_sbml"] = initial_sim["assignments_for_sbml"]

            reaction_process(0, "_SFS_0", initial_sim["reactions_for_sbml"])

            event_process("_SFS_1", 0, 0, initial_sim)

            process_simulations(multi_sims)

            for i in range(len(multi_sims)):
                new_sbml_file["species_for_sbml"]["_SFS_" + str(i)] = 0
            new_sbml_file["species_for_sbml"]["_SFS_0"] = 1

            list_of_composite_dicts_for_sbml.append([new_sbml_file])
        return list_of_composite_dicts_for_sbml

    def parse_volume_name_for_antimony(self) -> list[list[SBMLModelDict]]:
        new_sims = []
        for multi_sims in self.sbml_data_list:
            sim_sbml = multi_sims[0]
            new_sbml_file: SBMLModelDict = {
                "species_for_sbml": sim_sbml["species_for_sbml"],
                "parameters_for_sbml": sim_sbml["parameters_for_sbml"],
                "reactions_for_sbml": {},
                "events_for_sbml": sim_sbml["events_for_sbml"],
                "assignments_for_sbml": sim_sbml["assignments_for_sbml"],
            }

            new_sbml_file["parameters_for_sbml"]["_vol"] = sim_sbml[
                "parameters_for_sbml"
            ]["volume"]

            for re_name, reaction in sim_sbml["reactions_for_sbml"].items():
                new_reaction: ReactionData = {
                    "re": reaction["re"],
                    "pr": reaction["pr"],
                    "kin": reaction["kin"].replace("volume", "_vol"),
                }
                new_sbml_file["reactions_for_sbml"][re_name] = new_reaction

            new_sims.append([new_sbml_file])
        return new_sims

    def generate_sbml(self, compose: bool = False) -> list[str]:
        """
        Generates sbml strings from the current stored models in the simulation.

        :param compose: (bool) Join composite simulations into a single sbml
        :return: list of sbml strings from all stored simulations
        """
        to_return = []
        if self._species_for_sbml is None:
            self.compile(verbose=False)
        self._assemble_multi_simulation_structure()

        sbml_dict_list = self.compose_sbml() if compose else self.sbml_data_list

        for parameter_sweep in sbml_dict_list:
            for sbml_data in parameter_sweep:
                to_return.append(
                    sbml_build(
                        sbml_data["species_for_sbml"],
                        sbml_data["parameters_for_sbml"],
                        sbml_data["reactions_for_sbml"],
                        sbml_data["events_for_sbml"],
                        sbml_data["assignments_for_sbml"],
                    )
                )
        return to_return

    def generate_antimony(
        self, compose: bool = False, model_name: str | None = None
    ) -> list[str]:
        """
        Generates a string with an Antimony model from a respective MobsPy model.

        :param compose: (bool) Join composite simulations into a single sbml
        :param model_name: (str) desired name of the model
        """
        if self._species_for_sbml is None:
            self.compile(verbose=False)
        self._assemble_multi_simulation_structure()

        sbml_dict_list = (
            self.compose_sbml() if compose else self.parse_volume_name_for_antimony()
        )

        model_list = []

        for parameter_sweep in sbml_dict_list:
            for sbml_data in parameter_sweep:
                if model_name is None:
                    antimony_model = f"model mobspy_{rd_randint(0, 100000)} \n"
                else:
                    antimony_model = f"model {model_name} \n"
                if sbml_data["species_for_sbml"]:
                    for species_name, species_count in sbml_data[
                        "species_for_sbml"
                    ].items():
                        antimony_model = (
                            antimony_model
                            + f"    {species_name} = {species_count} dimensionless"
                        )
                        antimony_model = antimony_model + "\n"

                if sbml_data["parameters_for_sbml"]:
                    for parameter_name, parameter_value in sbml_data[
                        "parameters_for_sbml"
                    ].items():
                        if parameter_name == "volume":
                            continue
                        antimony_model = (
                            antimony_model
                            + f"    {parameter_name} = {parameter_value[0]} "
                            f"dimensionless\n"
                        )

                if sbml_data["assignments_for_sbml"]:
                    for _assign_name, assign_data in sbml_data[
                        "assignments_for_sbml"
                    ].items():
                        antimony_model = (
                            antimony_model + f"    {assign_data['species']}"
                            f" := {assign_data['expression']}\n"
                        )

                for reaction_name, reaction_data in sbml_data[
                    "reactions_for_sbml"
                ].items():
                    if "phantom" in reaction_name:
                        continue

                    antimony_model = antimony_model + f"    {reaction_name}: "
                    for i, r in enumerate(reaction_data["re"]):
                        if i == 0 and r[0] > 1:
                            antimony_model = antimony_model + f"{r[0]}*{r[1]}"
                            continue
                        if i == 0:
                            antimony_model = antimony_model + f"{r[1]}"
                            continue

                        if r[0] > 1:
                            antimony_model = antimony_model + f" + {r[0]} {r[1]}"
                        else:
                            antimony_model = antimony_model + f" + {r[1]}"

                    antimony_model = antimony_model + " -> "
                    for i, p in enumerate(reaction_data["pr"]):
                        if i == 0 and p[0] > 1:
                            antimony_model = antimony_model + f" {p[0]} {p[1]}"
                            continue
                        if i == 0:
                            antimony_model = antimony_model + f" {p[1]}"
                            continue

                        if p[0] > 1:
                            antimony_model = antimony_model + f" + {p[0]}*{p[1]}"
                        else:
                            antimony_model = antimony_model + f" + {p[1]}"

                    antimony_model = antimony_model + f"; {reaction_data['kin']}"

                    antimony_model = antimony_model + "\n"

                if sbml_data["events_for_sbml"]:
                    for event_name, event_data in sbml_data["events_for_sbml"].items():
                        if event_data["trigger"] == "true":
                            antimony_model = (
                                antimony_model + f"    {event_name}: at("
                                f"time > {event_data['delay']}): "
                            )
                        else:
                            antimony_model = (
                                antimony_model
                                + f"    {event_name}: at({event_data['trigger']}): "
                            )

                        for asg in event_data["assignments"]:
                            antimony_model = antimony_model + f" {asg[0]}={asg[1]},"
                        antimony_model = antimony_model[:-1] + "\n"

                antimony_model = antimony_model + "end"
                model_list.append(antimony_model)

        return model_list
