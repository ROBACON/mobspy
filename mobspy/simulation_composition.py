"""SimulationComposition: concatenation of multiple Simulation objects."""

from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING
from typing import Any as TypingAny

from pint import Quantity

from mobspy.exceptions import SimulationError
from mobspy.model_generation import ModelGenerationMixin
from mobspy.params.parameter_reader import (
    manually_process_each_parameter as pr_manually_process_each_parameter,
)
from mobspy.plotting import PlottingMixin
from mobspy.results.io import save_results
from mobspy.results.processing import process_results
from mobspy.results.time_series import SimulationResults
from mobspy.runtime import build_execution_plan
from mobspy.simulation_config import PlotConfig
from mobspy.types import CompiledModel, ExecutionPlan, SimulationBackend

if TYPE_CHECKING:
    from collections.abc import Generator

    from mobspy.simulation import Simulation


class SimulationComposition(ModelGenerationMixin, PlottingMixin):
    """Concatenation of multiple Simulation objects via the ``+`` operator.

    Represents a sequential pipeline where the end state of one simulation
    becomes the initial state of the next.
    """

    fres: SimulationResults | dict[str, TypingAny]
    _plan: ExecutionPlan | None
    _parameter_list_of_dic: list[dict[str, TypingAny]]

    def update_model(self, *args: TypingAny) -> None:
        """Delegate model updates to the base simulation."""
        self.base_sim.update_model(*args)

    def delete(self) -> None:
        """Release this composition's results without changing its children."""
        self.results = {}
        self.fres = {}
        self._plan = None
        self.sbml_data_list = []

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
        S1: Simulation | SimulationComposition,  # noqa: N803  # legacy name
        S2: Simulation | SimulationComposition,  # noqa: N803  # legacy name
    ) -> None:
        from mobspy.simulation import Simulation  # noqa: PLC0415  # circular import

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
        self.results: SimulationResults | dict[str, TypingAny] = {}
        self.fres = {}
        self.base_sim = self.list_of_simulations[0]
        self._plan = None
        self.sbml_data_list = []
        self._parameter_list_of_dic = []
        self._plot_parameters = PlotConfig(deepcopy(self.base_sim.plot_parameters))

    def __add__(
        self,
        other: Simulation | SimulationComposition,
    ) -> SimulationComposition:
        return SimulationComposition(self, other)

    _WHITE_LIST = frozenset(
        [
            "list_of_simulations",
            "results",
            "base_sim",
            "fres",
            "_plan",
            "sbml_data_list",
            "_parameter_list_of_dic",
            "_plot_parameters",
        ]
    )
    _MULTI_CAST_PARAMETERS = frozenset(["duration"])
    _BROAD_CAST_PARAMETERS = frozenset(
        ["level", "rate_type", "plot_type", "repetitions"]
    )
    _DOUBLE_CAST_PARAMETERS = frozenset(["simulation_method", "volume", "method"])

    def __setattr__(self, name: str, value: TypingAny) -> None:
        if name == "plot_parameters":
            self._plot_parameters = PlotConfig(deepcopy(value))
        elif name in self._DOUBLE_CAST_PARAMETERS:
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

        for par, sim in zip(value, self, strict=True):
            if name == "volume":
                sim.volume = par
            elif name == "duration":
                sim.duration = par
            else:
                sim._set_parameter(name, par)

    @property
    def plot_config(self) -> PlotConfig:
        """Access this composition's plot configuration."""
        return self._plot_parameters

    @property
    def _backend(self) -> SimulationBackend:
        return self.base_sim._backend

    @property
    def _species_for_sbml(self) -> TypingAny:
        return self.base_sim._species_for_sbml

    @property
    def _list_of_models(self) -> list[CompiledModel]:
        return [model for sim in self for model in sim._list_of_models]

    @property
    def _list_of_parameters(self) -> list[TypingAny]:
        return [sim._runtime_parameters for sim in self]

    @property
    def plot_parameters(self) -> dict[str, TypingAny]:
        return self._plot_parameters

    def _assemble_multi_simulation_structure(self) -> None:
        self._check_all_sims_compilation()
        self._compile_multi_simulation()
        self._plan = build_execution_plan(
            [sim.compiled_model for sim in self], self._list_of_parameters
        )
        self.sbml_data_list = [
            [
                model.to_compiled_model(dict(model.species), dict(model.mappings))
                for model in chain
            ]
            for chain in self._plan.models
        ]
        self._parameter_list_of_dic = list(self._plan.parameter_values)

    def compile(self, verbose: bool = True) -> str | None:
        """Compile each child and build an owned execution plan."""
        summaries = [sim.compile(verbose) for sim in self]
        self._assemble_multi_simulation_structure()
        return "".join(summary for summary in summaries if summary) or None

    def _check_all_sims_compilation(self) -> None:
        for sim in self:
            sim._ensure_compiled()

    def run(  # noqa: PLR0913  # complex function signature
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
    ) -> SimulationResults:
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

        self._assemble_multi_simulation_structure()
        if self._plan is None:
            raise SimulationError("Composition has no execution plan")
        jobs = self.base_sim.set_job_number(self.base_sim._runtime_parameters)
        raw_results = self._backend.run(self._plan, jobs=jobs)
        self.results, self.fres = process_results(
            raw_results, self._plan, self._list_of_parameters
        )
        parameters = self.base_sim._runtime_parameters
        for name in ("unit_x", "unit_y", "output_concentration"):
            self.plot_parameters[name] = parameters[name]
        if parameters["save_data"]:
            self.save_data()
        if parameters["plot_data"]:
            self.plot()
        return self.results

    def save_data(self, file: str | None = None) -> None:
        """Save this composition's results without changing its children."""
        save_results(
            self.results,
            file,
            self.base_sim._runtime_parameters.get("absolute_output_file"),
        )

    def generate_sbml(self, compose: bool = False) -> list[str]:
        self._assemble_multi_simulation_structure()
        return super().generate_sbml(compose=compose)

    def generate_antimony(
        self, compose: bool = False, model_name: str | None = None
    ) -> list[str]:
        self._assemble_multi_simulation_structure()
        return super().generate_antimony(compose=compose, model_name=model_name)

    def to_dataframe(self) -> TypingAny:
        """Convert composition results to a pandas DataFrame."""
        if not isinstance(self.results, SimulationResults):
            raise SimulationError("Run the composition before accessing results")
        return self.results.return_pandas()

    @classmethod
    def is_simulation(cls) -> bool:
        """Return True to identify this class as a simulation."""
        return True
