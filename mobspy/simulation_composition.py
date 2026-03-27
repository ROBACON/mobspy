"""
Simulation composition module: concatenating multiple MobsPy simulations.

Provides the SimulationComposition class that allows combining multiple
Simulation objects using the + operator for sequential simulation execution.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from pint import Quantity

from mobspy.mobspy_logging import get_logger
from mobspy.parameter_scripts.parameter_reader import (
    manually_process_each_parameter as pr_manually_process_each_parameter,
)
from mobspy.parameter_scripts.parametric_sweeps import (
    unite_parameter_dictionaries as ps_unite_parameter_dictionaries,
)

logger = get_logger(__name__)
simlog = logger


class SimulationComposition:
    def update_model(self, *args: Any) -> None:
        self.base_sim.update_model(*args)

    def _compile_multi_simulation(self) -> None:
        for sim1 in self.list_of_simulations:
            for sim2 in self.list_of_simulations:
                if sim1 == sim2:
                    continue

                for spe1 in sim1.model:
                    for spe2 in sim2.model:
                        if (
                            spe1.get_name() == spe2.get_name()
                            and spe1.get_all_characteristics()
                            != spe2.get_all_characteristics()
                        ):
                            simlog.error(
                                f"Species {spe1.get_name()} was "
                                "modified through simulations.\n"
                                "Although reactions can be removed,"
                                " the characteristics inherited "
                                "must remain the same"
                            )

    def __len__(self) -> int:
        return len(self.list_of_simulations)

    def __iter__(self) -> Any:
        yield from self.list_of_simulations

    def __init__(self, S1: Any, S2: Any) -> None:
        # Import here to avoid circular imports
        from mobspy.simulation import Simulation

        if isinstance(S1, Simulation) and isinstance(S2, Simulation):
            self.list_of_simulations = [S1] + [S2]
        elif isinstance(S1, SimulationComposition) and isinstance(S2, Simulation):
            self.list_of_simulations = S1.list_of_simulations + [S2]
        elif isinstance(S1, Simulation) and isinstance(S2, SimulationComposition):
            self.list_of_simulations = [S1] + S2.list_of_simulations
        elif isinstance(S1, SimulationComposition) and isinstance(
            S2, SimulationComposition
        ):
            self.list_of_simulations = S1.list_of_simulations + S2.list_of_simulations
        else:
            simlog.error(
                "Simulation compositions can only be performed with other simulations",
            )
        self.results: Any = None
        self.fres: Any = None
        self.base_sim: Any = self.list_of_simulations[0]

    def __add__(self, other: Any) -> SimulationComposition:
        return SimulationComposition(self, other)

    def __setattr__(self, name: str, value: Any) -> None:
        white_list = ["list_of_simulations", "results", "base_sim", "fres"]
        multi_cast_parameters = ["duration"]
        broad_cast_parameters = ["level", "rate_type", "plot_type", "repetitions"]
        double_cast_parameters = ["simulation_method", "volume", "method"]

        if name in double_cast_parameters:
            # Broadcast if single value
            if isinstance(value, (str, int, float, Quantity)):
                for sim in self:
                    sim.__dict__["parameters"][name] = value
            else:
                # Multicast if list
                try:
                    if not len(self) == len(value):
                        raise SystemExit
                except Exception:
                    simlog.error(
                        f"The parameter {name} was assigned non-accepted type.",
                    )

                # Don't add directly to __dict__; volume/duration
                # changes are checked in individual sim setattr
                for par, sim in zip(value, self, strict=False):
                    if name == "volume":
                        sim.volume = par
                    elif name == "duration":
                        sim.duration = par
                    else:
                        sim.__dict__["parameters"][name] = par

        elif name in multi_cast_parameters:
            try:
                if not len(self) == len(value):
                    raise SystemExit
            except Exception:
                simlog.error(
                    "From 2.4.4 duration must be assigned to each "
                    "simulation individually or a list with all "
                    "durations must be assigned to the "
                    "concatenated simulation",
                )

            # Don't add directly to __dict__; volume/duration
            # changes are checked in individual sim setattr
            for par, sim in zip(value, self, strict=False):
                if name == "volume":
                    sim.volume = par
                elif name == "duration":
                    sim.duration = par
                else:
                    sim.__dict__["parameters"][name] = par
        elif name in broad_cast_parameters:
            for sim in self:
                sim.__dict__["parameters"][name] = value
        else:
            if name in white_list:
                self.__dict__[name] = value
            else:
                self.base_sim.__setattr__(name, value)

    def __getattr__(self, item: str) -> Any:
        if item == "plot_config":
            self.base_sim.__dict__["plot_flag"] = True
            return self.base_sim

    def compile(self, verbose: bool = True) -> str | None:
        str = ""
        for sim in self.list_of_simulations:
            str += sim.compile(verbose)

        self._compile_multi_simulation()

        for sim in self.list_of_simulations:
            if sim == self.base_sim:
                continue

            self.base_sim._list_of_models += sim._list_of_models
            self.base_sim._list_of_parameters += sim._list_of_parameters

        self.base_sim._assemble_multi_simulation_structure()

        if str != "":
            return str
        return None

    def _check_all_sims_compilation(self) -> None:
        for sim in self.list_of_simulations:
            if sim._species_for_sbml is None:
                sim.compile(verbose=False)

    # This run is for the multiple simulations
    def run(
        self,
        duration: Any = None,
        volume: Any = None,
        dimension: Any = None,
        repetitions: int | None = None,
        level: int | None = None,
        simulation_method: str | None = None,
        start_time: float | None = None,
        r_tol: float | None = None,
        a_tol: float | None = None,
        seeds: list[int] | None = None,
        step_size: float | None = None,
        jobs: int | None = None,
        unit_x: Any = None,
        unit_y: Any = None,
        output_concentration: bool | None = None,
        output_event: bool | None = None,
        output_file: str | None = None,
        save_data: bool | None = None,
        plot_data: bool | None = None,
        rate_type: str | None = None,
        plot_type: str | None = None,
    ) -> None:
        """Run a concatenated simulation.

        Duration and volume must be iterables with one
        value per sub-simulation.
        """
        if level is not None:
            self.level = level

        self._check_all_sims_compilation()
        self._compile_multi_simulation()

        pr_manually_process_each_parameter(
            self,
            duration,
            volume,
            dimension,
            repetitions,
            level,
            simulation_method,
            start_time,
            r_tol,
            a_tol,
            seeds,
            step_size,
            jobs,
            unit_x,
            unit_y,
            output_concentration,
            output_event,
            output_file,
            save_data,
            plot_data,
            rate_type,
            plot_type,
        )

        multi_parameter_dictionary: dict[str, Any] = {}

        for sim in self.list_of_simulations:
            multi_parameter_dictionary = ps_unite_parameter_dictionaries(
                multi_parameter_dictionary, sim.model_parameters
            )

        self.base_sim.model_parameters = multi_parameter_dictionary

        for sim in self.list_of_simulations:
            if sim == self.base_sim:
                continue

            self.base_sim._list_of_models += sim._list_of_models
            self.base_sim._list_of_parameters += sim._list_of_parameters

        self.base_sim.run()
        self.results = self.base_sim.results
        self.fres = self.base_sim.fres

    def plot_deterministic(self, *species: Any) -> None:
        self.base_sim.plot_deterministic(*species)

    def plot_stochastic(self, *species: Any) -> None:
        self.base_sim.plot_stochastic(*species)

    def plot(self, *species: Any) -> None:
        self.base_sim.plot(*species)

    def plot_raw(self, parameters_or_file: str | dict[str, Any]) -> None:
        self.base_sim.plot_raw(parameters_or_file)

    def add_plot_params(self, *args: Any, **kwargs: Any) -> None:
        for a in args:
            if isinstance(a, dict):
                for par in a:
                    self.base_sim.plot_parameters[par] = a[par]

        for key in kwargs:
            self.base_sim.plot_parameters[key] = deepcopy(kwargs[key])

    def generate_sbml(self, compose: bool = False) -> list[str]:
        """
        Generates a string with an SBML model from a respective MobsPy model
        :param compose: (bool) Join composite simulations into a single sbml
        """

        self._check_all_sims_compilation()
        self._compile_multi_simulation()

        for sim in self.list_of_simulations:
            if sim == self.base_sim:
                continue

            self.base_sim._list_of_models += sim._list_of_models
            self.base_sim._list_of_parameters += sim._list_of_parameters

        return self.base_sim.generate_sbml(compose=compose)

    def generate_antimony(
        self, compose: bool = False, model_name: str | None = None
    ) -> list[str]:
        """
        Generates a string with an Antimony model from a respective MobsPy model
        :param compose: (bool) Join composite simulations into a single sbml
        """
        self._check_all_sims_compilation()
        self._compile_multi_simulation()

        for sim in self.list_of_simulations:
            if sim == self.base_sim:
                continue

            self.base_sim._list_of_models += sim._list_of_models
            self.base_sim._list_of_parameters += sim._list_of_parameters

        return self.base_sim.generate_antimony(compose=compose, model_name=model_name)

    def to_dataframe(self) -> None:
        self.base_sim.to_dataframe()

    @classmethod
    def is_simulation(cls) -> bool:
        return True
