"""SimulationComposition: concatenation of multiple Simulation objects."""

from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING
from typing import Any as TypingAny

from pint import Quantity

from mobspy.data_handler.time_series_object import SimulationResults
from mobspy.exceptions import SimulationError
from mobspy.modules.reactions import Reacting_Species
from mobspy.modules.species import Species
from mobspy.parameter_scripts.parameter_reader import (
    manually_process_each_parameter as pr_manually_process_each_parameter,
)
from mobspy.parameter_scripts.parametric_sweeps import (
    unite_parameter_dictionaries as ps_unite_parameter_dictionaries,
)
from mobspy.simulation_config import PlotConfig

if TYPE_CHECKING:
    from collections.abc import Generator

    from mobspy.simulation import Simulation


class SimulationComposition:
    """Concatenation of multiple Simulation objects via the ``+`` operator.

    Represents a sequential pipeline where the end state of one simulation
    becomes the initial state of the next.
    """

    def update_model(self, *args: TypingAny) -> None:
        """Delegate model updates to the base simulation."""
        self.base_sim.update_model(*args)

    def _compile_multi_simulation(self) -> None:
        """Validate shared species across simulations.

        Ensures consistent characteristics.
        """
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
                            raise SimulationError(
                                f"Species {spe1.get_name()} "
                                "was modified through "
                                "simulations. \n" + "Although reactions can be "
                                "removed, the characteristics "
                                "inherited must remain the same"
                            )

    def __len__(self) -> int:
        return len(self.list_of_simulations)

    def __iter__(self) -> Generator[Simulation, None, None]:
        yield from self.list_of_simulations

    def __init__(
        self,
        S1: Simulation | SimulationComposition,  # noqa: N803
        S2: Simulation | SimulationComposition,  # noqa: N803
    ) -> None:
        from mobspy.simulation import Simulation  # noqa: PLC0415

        if isinstance(S1, Simulation) and isinstance(S2, Simulation):
            self.list_of_simulations = [S1, S2]
        elif isinstance(S1, SimulationComposition) and isinstance(S2, Simulation):
            self.list_of_simulations = [*S1.list_of_simulations, S2]
        elif isinstance(S1, Simulation) and isinstance(S2, SimulationComposition):
            self.list_of_simulations = [S1, *S2.list_of_simulations]
        elif isinstance(S1, SimulationComposition) and isinstance(
            S2, SimulationComposition
        ):
            self.list_of_simulations = S1.list_of_simulations + S2.list_of_simulations
        else:
            raise SimulationError(
                "Simulation compositions can only be performed with other simulations"
            )
        self.results: SimulationResults | dict[str, TypingAny] | None = None
        self.fres: SimulationResults | dict[str, TypingAny] | None = None
        self.base_sim = self.list_of_simulations[0]

    def __add__(
        self,
        other: Simulation | SimulationComposition,
    ) -> SimulationComposition:
        return SimulationComposition(self, other)

    _WHITE_LIST = frozenset(["list_of_simulations", "results", "base_sim", "fres"])
    _MULTI_CAST_PARAMETERS = frozenset(["duration"])
    _BROAD_CAST_PARAMETERS = frozenset(
        ["level", "rate_type", "plot_type", "repetitions"]
    )
    _DOUBLE_CAST_PARAMETERS = frozenset(["simulation_method", "volume", "method"])

    def __setattr__(self, name: str, value: TypingAny) -> None:
        if name in self._DOUBLE_CAST_PARAMETERS:
            if isinstance(value, (str, int, float, Quantity)):
                for sim in self:
                    sim._set_parameter(name, value)
            else:
                self._multicast_parameter(name, value)
        elif name in self._MULTI_CAST_PARAMETERS:
            self._multicast_parameter(name, value, require_list=True)
        elif name in self._BROAD_CAST_PARAMETERS:
            for sim in self:
                sim._set_parameter(name, value)
        elif name in self._WHITE_LIST:
            self.__dict__[name] = value
        else:
            self.base_sim.__setattr__(name, value)

    def _multicast_parameter(
        self,
        name: str,
        value: TypingAny,
        *,
        require_list: bool = False,
    ) -> None:
        """Distribute a list of values across child simulations."""
        try:
            value_len = len(value)
        except TypeError as e:
            msg = (
                f"From 2.4.4 {name} must be assigned to each simulation "
                "individually or as a list"
                if require_list
                else f"The parameter {name} was assigned non-accepted type."
            )
            raise SimulationError(msg) from e

        if value_len != len(self):
            raise SimulationError(
                f"The parameter {name} list length must match "
                f"the number of simulations ({len(self)})."
            )

        for par, sim in zip(value, self, strict=False):
            if name == "volume":
                sim.volume = par
            elif name == "duration":
                sim.duration = par
            else:
                sim._set_parameter(name, par)

    @property
    def plot_config(self) -> PlotConfig:
        """Access plot configuration via the base simulation."""
        return self.base_sim.plot_config

    def compile(self, verbose: bool = True) -> str | None:
        """Compile all child simulations and merge their models."""
        result_str = ""
        for sim in self.list_of_simulations:
            compiled = sim.compile(verbose)
            if compiled is not None:
                result_str += compiled

        self._compile_multi_simulation()

        for sim in self.list_of_simulations:
            if sim == self.base_sim:
                continue

            self.base_sim._list_of_models += sim._list_of_models
            self.base_sim._list_of_parameters += sim._list_of_parameters

        self.base_sim._assemble_multi_simulation_structure()

        if result_str != "":
            return result_str
        return None

    def _check_all_sims_compilation(self) -> None:
        for sim in self.list_of_simulations:
            if sim._species_for_sbml is None:
                sim.compile(verbose=False)

    def run(  # noqa: PLR0913
        self,
        duration: TypingAny = None,
        volume: TypingAny = None,
        dimension: int | None = None,
        repetitions: int | None = None,
        level: int | None = None,
        simulation_method: str | None = None,
        start_time: float | None = None,
        r_tol: float | None = None,
        a_tol: float | None = None,
        seeds: list[int] | None = None,
        step_size: float | None = None,
        jobs: int | None = None,
        unit_x: TypingAny = None,
        unit_y: TypingAny = None,
        output_concentration: bool | None = None,
        output_event: bool | None = None,
        output_file: str | None = None,
        save_data: bool | None = None,
        plot_data: bool | None = None,
        rate_type: str | None = None,
        plot_type: str | None = None,
    ) -> None:
        """Run a concatenated simulation with multiple simulation objects.

        Some inputs are different here: duration and volume must receive any
        iterable with the value for each simulation.

        Args:
            duration: Duration of a simulation.
            volume: Volume of the simulation, if none given 1 liter is used.
            dimension: Spatial dimension.
            repetitions: Number of times to repeat.
            level: 0 - only error messages, 3 - errors, warnings, compilation info.
            simulation_method: Stochastic, deterministic, direct_method.
            start_time: Simulation will only display data after the start time.
            r_tol: Relative tolerance.
            a_tol: Absolute tolerance.
            seeds: Seeds for stochastic simulation.
            step_size: Time step-size.
            jobs: Number of jobs.
            unit_x: Unit of the time x-axis.
            unit_y: Unit of the y-axis.
            output_concentration: Outputs concentration instead of counts.
            output_event: When an event happens, adds the data point to the results.
            output_file: Name of the file.
            save_data: Save data or not.
            plot_data: Plot data or not.
            rate_type: Stochastic or deterministic rate expression.
            plot_type: Stochastic or deterministic style of MobsPy plot.
        """
        if level is not None:
            self.level = level

        self._check_all_sims_compilation()
        self._compile_multi_simulation()

        pr_manually_process_each_parameter(
            self,
            duration=duration,
            volume=volume,
            dimension=dimension,
            repetitions=repetitions,
            level=level,
            simulation_method=simulation_method,
            start_time=start_time,
            r_tol=r_tol,
            a_tol=a_tol,
            seeds=seeds,
            step_size=step_size,
            jobs=jobs,
            unit_x=unit_x,
            unit_y=unit_y,
            output_concentration=output_concentration,
            output_event=output_event,
            output_file=output_file,
            save_data=save_data,
            plot_data=plot_data,
            rate_type=rate_type,
            plot_type=plot_type,
        )

        multi_parameter_dictionary: dict[str, TypingAny] = {}

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

    def plot_deterministic(self, *species: str | Species | Reacting_Species) -> None:
        """Plot deterministic results via the base simulation."""
        self.base_sim.plot_deterministic(*species)

    def plot_stochastic(self, *species: str | Species | Reacting_Species) -> None:
        """Plot stochastic results via the base simulation."""
        self.base_sim.plot_stochastic(*species)

    def plot(self, *species: str | Species | Reacting_Species) -> None:
        """Plot results using the default plot type."""
        self.base_sim.plot(*species)

    def plot_raw(self, parameters_or_file: str | dict[str, TypingAny]) -> None:
        """Plot results with raw user-supplied parameters."""
        self.base_sim.plot_raw(parameters_or_file)

    def add_plot_params(self, *args: TypingAny, **kwargs: TypingAny) -> None:
        """Merge additional plotting parameters into the composition."""
        for a in args:
            if isinstance(a, dict):
                for par in a:
                    self.base_sim.plot_parameters[par] = a[par]

        for key in kwargs:  # noqa: PLC0206
            self.base_sim.plot_parameters[key] = deepcopy(kwargs[key])

    def generate_sbml(self, compose: bool = False) -> list[str]:
        """Generate SBML model strings from a composed MobsPy model.

        Args:
            compose: Join composite simulations into a single sbml.
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
        self,
        compose: bool = False,
        model_name: str | None = None,
    ) -> list[str]:
        """Generate Antimony model strings from a composed MobsPy model.

        Args:
            compose: Join composite simulations into a single sbml.
            model_name: Optional name for the Antimony model.
        """
        self._check_all_sims_compilation()
        self._compile_multi_simulation()

        for sim in self.list_of_simulations:
            if sim == self.base_sim:
                continue

            self.base_sim._list_of_models += sim._list_of_models
            self.base_sim._list_of_parameters += sim._list_of_parameters

        return self.base_sim.generate_antimony(compose=compose, model_name=model_name)

    def to_dataframe(self) -> TypingAny:
        """Convert composition results to a pandas DataFrame."""
        self.base_sim.to_dataframe()

    @classmethod
    def is_simulation(cls) -> bool:
        """Return True to identify this class as a simulation."""
        return True
