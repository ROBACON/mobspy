"""Main MobsPy module.

Contains the Simulation class which is responsible for
simulating a Model.
"""

from __future__ import annotations

import logging
from copy import deepcopy
from dataclasses import fields
from json import dump as json_dump
from json import load as json_load
from os.path import splitext as os_path_splitext
from typing import TYPE_CHECKING
from typing import Any as TypingAny

import joblib
from pint import Quantity

from mobspy.constants import DOT_SEPARATOR
from mobspy.data_handler.process_result_data import (
    convert_data_to_desired_unit as dh_convert_data_to_desired_unit,
)
from mobspy.data_handler.process_result_data import (
    extract_time_and_volume_list as dh_extract_time_and_volume_list,
)
from mobspy.data_handler.time_series_object import (
    MobsPyList_of_TS,
    MobsPyTimeSeries,
)
from mobspy.event_handling import EventHandlingMixin
from mobspy.exceptions import (
    CompilationError,
    MobsPyError,
    ParameterError,
    ReactionError,  # noqa: F401
    SimulationError,
    ValidationError,
)
from mobspy.mobspy_logging import get_logger
from mobspy.model_generation import ModelGenerationMixin
from mobspy.modules.any_species import (
    Any as Any,
)
from mobspy.modules.assignments_implementation import (
    Assign,
)
from mobspy.modules.compiler import Compiler
from mobspy.modules.logic_operators import (
    MetaSpeciesLogicResolver as lop_MetaSpeciesLogicResolver,
)
from mobspy.modules.meta_class import (
    BaseSpecies,
    List_Species,
    ListSpecies,
    New,
    Reacting_Species,
    Species,
    Zero,
)
from mobspy.modules.mobspy_expressions import u
from mobspy.modules.mobspy_parameters import (
    Internal_Parameter_Constructor as _ParameterConstructor,
)
from mobspy.modules.mobspy_parameters import (
    ModelParameters,
)
from mobspy.modules.model_unit_context import ModelUnitContext
from mobspy.modules.order_operators import (
    All,
    Default,
    Rev,
    Set,
)
from mobspy.modules.set_counts_module import set_counts
from mobspy.modules.species_utils import (
    create_orthogonal_vector_structure as mcu_create_orthogonal_vector_structure,
)
from mobspy.modules.unit_handler import (
    extract_length_dimension as uh_extract_length_dimension,
)
from mobspy.parameter_estimation_data_loader.data_loader import (
    Experimental_Data_Holder as pdl_Experimental_Data_Holder,
)
from mobspy.parameter_estimation_data_loader.parameter_estimation_scripts import (
    basiCO_parameter_estimation,
)
from mobspy.parameter_scripts.parameter_reader import (
    convert_time_parameters_after_compilation as pr_convert_time_parameters_after_compilation,  # noqa: E501
)
from mobspy.parameter_scripts.parameter_reader import (
    convert_volume_after_compilation as pr_convert_volume_after_compilation,
)
from mobspy.parameter_scripts.parameter_reader import (
    manually_process_each_parameter as pr_manually_process_each_parameter,
)
from mobspy.parameter_scripts.parameter_reader import (
    parameter_process as pr_parameter_process,
)
from mobspy.parameter_scripts.parameter_reader import (
    read_json as pr_read_json,
)
from mobspy.parameter_scripts.parametric_sweeps import (
    generate_all_sbml_models as ps_generate_all_sbml_models,
)
from mobspy.parameter_scripts.parametric_sweeps import (
    unite_parameter_dictionaries as ps_unite_parameter_dictionaries,
)
from mobspy.plot_params.default_plot_reader import (
    get_default_plot_parameters,
)
from mobspy.plotting import PlottingMixin
from mobspy.sbml_simulator.run import simulate as sbml_simulate
from mobspy.simulation_config import SimulationConfig
from mobspy.simulator_object.utils import (
    Simulation_Utils,
)
from mobspy.simulator_object.utils import (
    sim_remove_reaction as sof_sim_remove_reaction,
)
from mobspy.types import SimulationEventData, TimeSeriesDataDict

if TYPE_CHECKING:
    from collections.abc import Generator

    from mobspy.types import (
        CompiledModelDict,
        EventsForSbml,
        MappingsForSbml,
        ParametersForSbml,
        ParameterSweepList,
        ReactionsForSbml,
        SpeciesForSbml,
    )

__all__ = [
    "All",
    "Any",
    "Assign",
    "BaseSpecies",
    "Default",
    "ListSpecies",
    "ModelParameters",
    "New",
    "Rev",
    "Set",
    "Simulation",
    "SimulationComposition",
    "Zero",
    "basiCO_parameter_estimation",
    "set_counts",
    "simlog",
    "u",
]

_logger = get_logger(__name__)
simlog = _logger


class PlotConfigProxy:
    """Proxy object for setting plot parameters via sim.plot_config.param = value."""

    def __init__(self, plot_params: dict[str, TypingAny]) -> None:
        object.__setattr__(self, "_plot_params", plot_params)

    def __setattr__(self, name: str, value: TypingAny) -> None:
        self._plot_params[name] = value

    def __getattr__(self, name: str) -> TypingAny:
        return self._plot_params.get(name)


class Simulation(
    pdl_Experimental_Data_Holder,
    Simulation_Utils,
    EventHandlingMixin,
    ModelGenerationMixin,
    PlottingMixin,
):
    """Orchestrates compilation, execution, and result handling for a MobsPy model.

    Collects meta-species and reactions, compiles them via the Compiler into
    SBML, runs the simulation through BasiCO/COPASI, and stores the results.

    Examples:
        >>> from mobspy import *
        >>> A, B = BaseSpecies(['A', 'B'])
        >>> _ = A >> B[1]
        >>> A(10)  # doctest: +ELLIPSIS
        <...>
        >>> S = Simulation(A | B)
        >>> S.duration = 5
        >>> _ = S.compile(verbose=False)
        >>> S._is_compiled
        True
    """

    def __init__(
        self,
        model: Species | List_Species,
        reactions: set | None = None,
        names: dict | None = None,
        parameters: dict | None = None,
        plot_parameters: dict | None = None,
    ) -> None:
        """
        Constructor of the simulation object.

        Initialize a new simulation with the specified model and parameters.

        Args:
            model: Meta-species object or list of meta-species for modeling
            reactions: Optional set of reactions to include (None for all reactions)
            names: Optional dictionary of meta-species names in globals() format
            parameters: Optional dictionary of simulation parameters
            plot_parameters: Optional dictionary of plotting parameters

        Raises:
            ValidationError: If model contains invalid species types
            ParameterError: If required parameters are missing
        """
        super().__init__()

        # Event Variable Definitions
        self._event_time = 0
        self.previous_trigger = None
        self.current_event_count_data = []
        self.total_packed_events: list[SimulationEventData] = []
        self.number_of_context_comparisons = 0
        self.pre_number_of_context_comparisons = 0
        self._list_of_models: list[CompiledModelDict] = []
        self._list_of_parameters: list[TypingAny] = []
        self._context_not_active = True
        self._assigned_species_list: list[str] = []
        self._conditional_event = False
        self._end_condition = None
        self.model_parameters = {}
        self.sbml_data_list: ParameterSweepList = []
        self._parameter_list_of_dic: list[dict[str, TypingAny]] = []
        self._is_compiled = False
        self.dimension: int | None = None
        self.model_parameter_objects_dict: dict[str, TypingAny] | None = None

        # Must copy to avoid reference assignment
        # Change model to include linked species
        model_pre_link = List_Species(model)  # type: ignore[arg-type]
        model_pos_link: set[Species] = set()
        for spe in model_pre_link:
            model_pos_link.add(spe)
            model_pos_link = model_pos_link.union(spe._linked_species)
        self.model = List_Species(model_pos_link)

        self.names = names

        if not isinstance(model, Species) and not isinstance(model, List_Species):
            raise ValidationError(
                "Model must be formed only by Species objects "
                "or List_Species objects. "
                f"Received type {type(model)} with value {model}"
            )

        self.orthogonal_vector_structure = mcu_create_orthogonal_vector_structure(model)  # type: ignore[arg-type]
        if reactions is not None:
            self._reactions_set = set(reactions)
        else:
            self._reactions_set = set()
            for spe_object in self.model:
                for reference in spe_object.get_references():
                    self._reactions_set = self._reactions_set.union(
                        reference.get_reactions()
                    )

        self._species_counts = []
        for spe_object in self.model:
            for count in spe_object.get_quantities():
                self._species_counts.append(
                    {
                        "object": spe_object,
                        "characteristics": count["characteristics"],
                        "quantity": count["quantity"],
                    }
                )

        if not parameters:
            self.parameters = SimulationConfig()  # type: ignore[assignment]
        else:
            config = SimulationConfig()
            config.update(parameters)
            self.parameters = config  # type: ignore[assignment]

        if not plot_parameters:
            self.plot_parameters = get_default_plot_parameters()

        # Other needed things for simulating
        self.results: MobsPyList_of_TS | dict[str, TypingAny] = {}
        self.fres: MobsPyList_of_TS | dict[str, TypingAny] = {}
        self.default_order = Default

        self._species_for_sbml: SpeciesForSbml | None = None
        self._reactions_for_sbml: ReactionsForSbml | None = None
        self._parameters_for_sbml: ParametersForSbml | None = None
        self._mappings_for_sbml: MappingsForSbml | None = None
        self._events_for_sbml: EventsForSbml | None = None
        self._model_context: ModelUnitContext | None = None
        self.model_string = ""

    def _set_parameter(self, name: str, value: TypingAny) -> None:
        """Set a simulation parameter directly, bypassing __setattr__."""
        self.__dict__["parameters"][name] = value

    def compile(self, verbose: bool = True) -> str | None:
        """
        Compile the chemical reaction network model into executable SBML format.

        This method processes the meta-species, reactions, and parameters to generate
        an SBML representation that can be executed by the simulation engine.

        Args:
            verbose (bool): If True, print detailed compilation
                information. Defaults to True.

        Returns:
            Optional[str]: The compiled model string in SBML format, or None if
                          compilation failed or no model was generated.

        Raises:
            CompilationError: If the model contains syntax errors or invalid constructs.
            ParameterError: If required parameters are missing or invalid.
            ValidationError: If the model structure is invalid.
        """
        try:
            if self.dimension is None:
                if isinstance(self.volume, Quantity):
                    self.dimension = uh_extract_length_dimension(
                        str(self.volume.dimensionality), self.dimension
                    )
                else:
                    self.dimension = 3

            # Set log level based on parameters
            log_level = self.parameters.get("level", logging.INFO)
            _logger.set_log_level(log_level)

            pr_parameter_process(self.parameters)  # type: ignore[arg-type]
            if self.parameters["method"] is not None:
                self.parameters["simulation_method"] = self.parameters["method"]

            if self.parameters["simulation_method"].lower() == "deterministic":
                self.plot_parameters["simulation_method"] = "deterministic"
            elif self.parameters["simulation_method"].lower() == "stochastic":
                self.plot_parameters["simulation_method"] = "stochastic"
            else:
                raise ParameterError(
                    "Invalid simulation method: "
                    f"{self.parameters['simulation_method']}. "
                    "Must be 'deterministic' or 'stochastic'"
                )

            # Pass end condition to dict parameters
            # It is stored outside of parameters to keep
            # the parameters serializable. However, it is
            # necessary for compilation so it is passed here.

            self.parameters["_end_condition"] = self._end_condition

            _model_context = ModelUnitContext.from_simulation(
                volume=self.parameters["volume"],
                duration=self.parameters["duration"],
                species_counts=self._species_counts,
                dimension=self.dimension,
            )

            _result = Compiler.compile(
                self.model,
                reactions_set=self._reactions_set,
                species_counts=self._species_counts,
                orthogonal_vector_structure=self.orthogonal_vector_structure,
                volume=self.parameters["volume"],
                dimension=self.dimension,
                type_of_model=self.parameters.get("rate_type") or "stochastic",
                verbose=verbose,
                event_dictionary=self.total_packed_events,
                continuous_sim=self.parameters["_continuous_simulation"],
                ending_condition=self.parameters["_end_condition"],
                skip_expression_check=self.parameters["skip_expression_check"],
                parameter_context=dict(_ParameterConstructor.parameter_stack),
                model_context=_model_context,
            )
        except MobsPyError:
            raise
        except (TypeError, ValueError, KeyError, AttributeError) as e:
            raise CompilationError(f"Model compilation failed: {e!s}") from e

        self._species_for_sbml = _result.species_for_sbml
        self._reactions_for_sbml = _result.reactions_for_sbml
        self._parameters_for_sbml = _result.parameters_for_sbml
        self._mappings_for_sbml = _result.mappings_for_sbml
        self.model_string = _result.model_str
        self._events_for_sbml = _result.events_for_sbml
        self._assigned_species_list = _result.assigned_species
        self.model_parameters = _result.parameters_used
        self.model_parameter_objects_dict = _result.parameter_object_dict
        self._assignments_for_sbml = _result.assignments_for_sbml
        self._has_mole = _result.has_mole
        self._model_context = _result.model_context

        # The volume is converted to the proper unit at the compiler level
        self.parameters["volume"] = self._parameters_for_sbml["volume"][0]
        self.mappings = deepcopy(self._mappings_for_sbml)

        self.all_species_not_mapped = {}
        for key in self._species_for_sbml:
            self.all_species_not_mapped[key.replace(DOT_SEPARATOR, ".")] = (
                self._species_for_sbml[key]
            )

        compiled_model = _result.to_compiled_model_dict(
            species_not_mapped=self.all_species_not_mapped,
            mappings=self.mappings,
        )
        self._list_of_models += [compiled_model]

        self._list_of_parameters = [self.parameters]

        self._is_compiled = True

        if self.model_string == "":
            return None

        return self.model_string

    def _assemble_multi_simulation_structure(self) -> None:
        data_for_sbml_construction: ParameterSweepList
        data_for_sbml_construction, parameter_list_of_dic = ps_generate_all_sbml_models(
            self.model_parameters, self._list_of_models
        )
        self.sbml_data_list = data_for_sbml_construction
        self._parameter_list_of_dic = parameter_list_of_dic

    def _process_run_parameters(
        self,
        *,
        duration: float | Quantity | None = None,
        volume: float | Quantity | None = None,
        dimension: int | None = None,
        repetitions: int | None = None,
        level: int | None = None,
        simulation_method: str | None = None,
        rate_type: str | None = None,
        plot_type: str | None = None,
        start_time: float | None = None,
        r_tol: float | None = None,
        a_tol: float | None = None,
        seeds: list[int] | None = None,
        step_size: float | None = None,
        jobs: int | None = None,
        unit_x: Quantity | None = None,
        unit_y: Quantity | None = None,
        output_concentration: bool | None = None,
        output_event: bool | None = None,
        output_file: str | None = None,
        save_data: bool | None = None,
        plot_data: bool | None = None,
    ) -> None:
        """Process and apply run-time parameter overrides, then ensure compilation."""
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

        if level is not None:
            self.level = level

        # Base case - If there are no events we compile the model here
        if self._species_for_sbml is None:
            self.compile(verbose=False)

        self._assemble_multi_simulation_structure()

    def _execute_simulations(self) -> tuple[list[TypingAny], int]:
        """Run SBML simulations via joblib and return raw results with job count."""
        jobs = self.set_job_number(self.parameters)  # type: ignore[arg-type]

        def simulation_function(x: TypingAny) -> TypingAny:
            """Run a single SBML simulation batch."""
            return sbml_simulate(jobs, self._list_of_parameters, x)

        results: list[TypingAny] = list(  # pyright: ignore[reportArgumentType]
            joblib.Parallel(n_jobs=jobs, prefer="threads")(
                joblib.delayed(simulation_function)(sbml)
                for sbml in self.sbml_data_list
            )
        )

        return results, jobs

    def _convert_and_store_results(
        self,
        raw_results: TypingAny,
        jobs: int,
    ) -> None:
        """Convert raw time-series data to desired units and store in self.results."""
        # Auto-set unit_y when user used molar units.
        # When model_context is active and substance is molar, data from COPASI
        # is already in moles - no auto-conversion needed.
        _substance_is_molar = (
            self._model_context is not None and self._model_context.substance_is_molar
        )
        if not _substance_is_molar:
            if (
                self.parameters["unit_y"] is None
                and self.parameters["output_concentration"]
                and self._has_mole
            ):
                self.parameters["unit_y"] = 1 * u.unit_registry_object.molar
            elif (
                self.parameters["unit_y"] is None
                and not self.parameters["output_concentration"]
                and self._has_mole
            ):
                self.parameters["unit_y"] = 1 * u.unit_registry_object.mol

        # Volume list and time list are to convert into concentrations
        # This section also checks if the output_concentration parameter is valid
        volume_list, time_list, flag_concentration = dh_extract_time_and_volume_list(
            self._list_of_parameters
        )
        _unit_y = self.parameters["unit_y"]
        tcb = (
            _unit_y is not None and "[length]" not in _unit_y.dimensionality  # type: ignore[union-attr,attr-defined]
        )
        if not flag_concentration or tcb:
            self.parameters["output_concentration"] = False

        def convert_one_ts_to_desired_unit(unconverted_data: TypingAny) -> TypingAny:
            """Convert a single time series to the requested units."""
            return dh_convert_data_to_desired_unit(
                unconverted_data,
                time_list,
                volume_list,
                self.parameters["unit_x"],
                self.parameters["unit_y"],
                self.parameters["output_concentration"],
                model_context=self._model_context,
            )

        def convert_all_ts_to_correct_format(
            single_ts: TypingAny,
            parameters: TypingAny,
            unit_convert: bool = False,
        ) -> TypingAny:
            """Wrap a time series into a MobsPyTimeSeries object."""
            data_dict = TimeSeriesDataDict(
                data=convert_one_ts_to_desired_unit(single_ts)
                if unit_convert
                else single_ts,
                params=self.parameters,
                models=self._list_of_models,
            )
            return MobsPyTimeSeries(data_dict, parameters)

        flatt_ts: list[tuple[TypingAny, TypingAny]] = []
        if self._parameter_list_of_dic:
            for r, params in zip(raw_results, self._parameter_list_of_dic):
                for ts in r:  # pyright: ignore[reportOptionalIterable]
                    flatt_ts.append((ts, params))
        else:
            for r in raw_results:
                for ts in r:  # pyright: ignore[reportOptionalIterable]
                    flatt_ts.append((ts, {}))

        ta = self.parameters["unit_x"] is not None
        tb = self.parameters["unit_y"] is not None
        tc = self.parameters["output_concentration"] if flag_concentration else False

        if ta or tb or tc:
            all_processed_data = joblib.Parallel(n_jobs=jobs, prefer="threads")(
                joblib.delayed(convert_all_ts_to_correct_format)(ts, params, True)
                for ts, params in flatt_ts
            )
        else:
            all_processed_data = joblib.Parallel(n_jobs=jobs, prefer="threads")(
                joblib.delayed(convert_all_ts_to_correct_format)(ts, params, False)
                for ts, params in flatt_ts
            )

        self.results = MobsPyList_of_TS(
            all_processed_data,  # pyright: ignore[reportArgumentType]
            self.model_parameter_objects_dict,  # pyright: ignore[reportArgumentType]
        )
        self.fres = MobsPyList_of_TS([all_processed_data[0]], None, True)  # pyright: ignore[reportArgumentType, reportIndexIssue]

    def run(
        self,
        duration: float | Quantity | None = None,
        volume: float | Quantity | None = None,
        dimension: int | None = None,
        repetitions: int | None = None,
        level: int | None = None,
        simulation_method: str | None = None,
        rate_type: str | None = None,
        plot_type: str | None = None,
        start_time: float | None = None,
        r_tol: float | None = None,
        a_tol: float | None = None,
        seeds: list[int] | None = None,
        step_size: float | None = None,
        jobs: int | None = None,
        unit_x: Quantity | None = None,
        unit_y: Quantity | None = None,
        output_concentration: bool | None = None,
        output_event: bool | None = None,
        output_file: str | None = None,
        save_data: bool | None = None,
        plot_data: bool | None = None,
    ) -> None:
        """
        Execute the simulation with specified parameters.

        This method runs the chemical reaction network simulation using the provided
        parameters. Parameters not explicitly set will use previously configured values.

        Args:
            duration: Duration of the simulation (time units or Quantity)
            volume: Volume of the simulation system (liter units or Quantity)
            dimension: Spatial dimension of the simulation (0, 1, 2, or 3)
            repetitions: Number of times to repeat the simulation
            level: Logging level (0=errors only, 3=verbose)
            simulation_method: Simulation method ('stochastic', 'deterministic', etc.)
            rate_type: Rate expression type ('stochastic' or 'deterministic')
            plot_type: Type of plot to generate ('stochastic' or 'deterministic')
            start_time: Time to start displaying results
            r_tol: Relative tolerance for ODE solver
            a_tol: Absolute tolerance for ODE solver
            seeds: List of random seeds for stochastic simulations
            step_size: Time step size for simulation
            jobs: Number of parallel jobs for simulation
            unit_x: Unit for time axis
            unit_y: Unit for concentration/amount axis
            output_concentration: Whether to output concentration instead of counts
            output_event: Whether to output data points at event times
            output_file: Name of file to save results
            save_data: Whether to save simulation data
            plot_data: Whether to generate plots

        Raises:
            SimulationError: If simulation execution fails
            ParameterError: If invalid parameters are provided
            CompilationError: If model needs recompilation but fails
        """
        self._process_run_parameters(
            duration=duration,
            volume=volume,
            dimension=dimension,
            repetitions=repetitions,
            level=level,
            simulation_method=simulation_method,
            rate_type=rate_type,
            plot_type=plot_type,
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
        )

        raw_results, num_jobs = self._execute_simulations()
        self._convert_and_store_results(raw_results, num_jobs)

        if self.parameters["save_data"]:
            self.save_data()

        # Set common simulation and plot parameters
        self.plot_parameters["unit_x"] = self.parameters["unit_x"]
        self.plot_parameters["unit_y"] = self.parameters["unit_y"]
        self.plot_parameters["output_concentration"] = self.parameters[
            "output_concentration"
        ]

        if self.parameters["plot_data"]:
            methods_list = [x["plot_type"] for x in self._list_of_parameters]

            if len(self._parameter_list_of_dic) > 1:
                self.plot_parametric()
                return

            if "stochastic" in methods_list:
                self.plot_stochastic()
            else:
                self.plot_deterministic()
        return

    def save_data(self, file: str | None = None) -> None:
        """
        Save simulation result data to a file in JSON format.

        Args:
            file: Optional name of the file to save data to. If None, uses default name.
                  If provided without .json extension, it will be added.

        Raises:
            SimulationError: If simulation results are not available
            IOError: If file writing fails
        """
        if not hasattr(self, "results") or not self.results:
            raise SimulationError("No simulation results available to save")

        self._save_data(file=file)

    def _save_data(self, file: str | None = None) -> None:
        """
        Save results manually into file. Useful for jupyter notebook users.

        Args:
            file: Optional name of the file to create and save JSON data

        Raises:
            IOError: If file writing fails
            SimulationError: If results are not available
        """
        if not hasattr(self, "results") or not self.results:
            raise SimulationError("No simulation results available to save")

        try:
            if file is None:
                if "absolute_output_file" not in self.parameters:
                    raise ParameterError(
                        "No default output file specified in parameters"
                    )
                out_path = self.parameters["absolute_output_file"]
                with open(out_path, "w", encoding="utf-8") as f:
                    json_dump(self.results.to_dict(), f, indent=4)  # type: ignore[union-attr]
            else:
                # Add .json extension if not present
                if not file.endswith(".json"):
                    file += ".json"
                with open(file, "w", encoding="utf-8") as jf:
                    json_dump(self.results.to_dict(), jf, indent=4)  # type: ignore[union-attr]
                    _logger.info(f"Successfully saved simulation results to {file}")
        except OSError as e:
            raise SimulationError(f"Error saving data to file: {e!s}") from e
        except (TypeError, ValueError) as e:
            raise SimulationError(f"Error serializing simulation data: {e!s}") from e

    def _pack_data(self, time_series_data: TypingAny) -> None:
        """
        Packs data from multiple simulations or external data into one simulation object

        Args:
            time_series_data: Data to be packed in the simulation object.
        """
        self.packed_data.append(time_series_data)

    # Dealing with parameters
    def set_from_json(self, file_name: str) -> None:
        """Set simulation parameters from a JSON file.

        Only keys matching public simulation parameters are accepted.
        Unknown or internal keys are rejected.

        Args:
            file_name: Path to the JSON file.

        Raises:
            ParameterError: If the file contains unknown parameter keys.
        """
        with open(file_name, encoding="utf-8") as json_file:
            data = json_load(json_file)
            for key in data:
                if key not in self._SIMULATION_PARAMS:
                    raise ParameterError(f"Unknown parameter in JSON config: {key!r}")
                self.__setattr__(key, data[key])

    _INTERNAL_ATTRS: frozenset[str] = frozenset(
        {
            "default_order",
            "volume",
            "model",
            "names",
            "parameters",
            "model_string",
            "plot_parameters",
            "results",
            "_species_for_sbml",
            "_reactions_for_sbml",
            "_parameters_for_sbml",
            "_mappings_for_sbml",
            "mappings",
            "all_species_not_mapped",
            "event_times",
            "event_models",
            "event_count_dics",
            "_events_for_sbml",
            "total_packed_events",
            "species_initial_counts",
            "_event_time",
            "previous_trigger",
            "current_event_count_data",
            "current_condition",
            "current_event_trigger_data",
            "number_of_context_comparisons",
            "pre_number_of_context_comparisons",
            "_continuous_simulation",
            "initial_duration",
            "_reactions_set",
            "_list_of_models",
            "_list_of_parameters",
            "_context_not_active",
            "_species_counts",
            "_assigned_species_list",
            "_conditional_event",
            "_end_condition",
            "orthogonal_vector_structure",
            "model_parameters",
            "fres",
            "sbml_data_list",
            "_parameter_list_of_dic",
            "_is_compiled",
            "dimension",
            "experimental_data",
            "model_parameter_objects_dict",
            "_assignments_for_sbml",
            "_has_mole",
            "_model_context",
        }
    )

    _SIMULATION_PARAMS: frozenset[str] = frozenset(
        f.name for f in fields(SimulationConfig) if not f.name.startswith("_")
    )

    def __setattr__(self, name: str, value: TypingAny) -> None:
        # Internal attributes: store directly on instance
        if name in self._INTERNAL_ATTRS:
            self.__dict__[name] = value
            # Some attrs (e.g. volume) are also simulation parameters
            if name not in self._SIMULATION_PARAMS:
                return

        # Simulation parameters
        if name in self._SIMULATION_PARAMS:
            if self._is_compiled and name not in {"unit_x", "unit_y"}:
                value = pr_convert_time_parameters_after_compilation(
                    value, model_context=self._model_context
                )
            if self._is_compiled and name == "volume":
                value = pr_convert_volume_after_compilation(
                    self.dimension,
                    self._parameters_for_sbml,  # type: ignore[arg-type]
                    value,
                    model_context=self._model_context,
                )

            if name == "duration":
                if isinstance(value, bool):
                    raise SimulationError(
                        "MobsPy has received an invalid "
                        f"trigger type: {type(value)} \n"
                        + "Please make sure you are not using the operator == for "
                        + "creating event conditions \n"
                    )
                if isinstance(value, lop_MetaSpeciesLogicResolver):
                    self._set_parameter("_continuous_simulation", True)
                    self.__dict__["_end_condition"] = value
                    if not self.__dict__["parameters"].get(
                        "initial_conditional_duration"
                    ):
                        self._set_parameter("initial_conditional_duration", 1)
                    return

            self._set_parameter(name, value)
            return

        if name in self._INTERNAL_ATTRS:
            return

        raise ParameterError(f"Parameter {name} is not supported")

    def __getattribute__(self, item: str) -> TypingAny:
        if item == "results" and self.__dict__["results"] == {}:
            raise SimulationError(
                "The results were accessed before the execution of the simulation"
            )
        if item == "fres" and self.__dict__["fres"] == {}:
            raise SimulationError(
                "The results were accessed before the execution of the simulation"
            )
        return super().__getattribute__(item)

    @property
    def plot_config(self) -> PlotConfigProxy:
        """Access plot configuration. Usage: sim.plot_config.param = value."""
        return PlotConfigProxy(self.__dict__["plot_parameters"])

    def __getattr__(self, item: str) -> TypingAny:
        """Fallback attribute access: look up simulation parameters."""
        params = self.__dict__.get("parameters")
        if params is not None and item in params:
            return params[item]
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{item}'"
        )

    def configure_parameters(self, config: str | dict[str, TypingAny]) -> None:
        """
        Configure simulation parameters from json file or dictionary

        Args:
            config: Path to a json file or a dictionary of parameters.
        """
        config_dict = self.__config_parameters(config)
        new_config = SimulationConfig()
        new_config.update(config_dict)
        self.parameters = new_config  # type: ignore[assignment]

    def configure_plot_parameters(self, config: str | dict[str, TypingAny]) -> None:
        """
        Configure plot parameters from json file or dictionary

        Args:
            config: Path to a json file or a dictionary of parameters.
        """
        self.plot_parameters = self.__config_parameters(config)

    @staticmethod
    def __config_parameters(config: str | dict[str, TypingAny]) -> dict[str, TypingAny]:
        """Shared helper for configure_parameters and configure_plot_parameters."""
        if isinstance(config, str):
            if os_path_splitext(config)[1] != ".json":
                raise ParameterError("Wrong file extension")
            parameters_to_config: dict[str, TypingAny] = pr_read_json(config)
        elif isinstance(config, dict):
            parameters_to_config = config
        else:
            raise ParameterError("Parameters must be python dictionary or json file")
        return parameters_to_config

    def add_plot_params(self, *args: TypingAny, **kwargs: TypingAny) -> None:
        """Merge additional plotting parameters into this simulation."""
        for a in args:
            if isinstance(a, dict):
                for par in a:
                    self.base_sim.plot_parameters[par] = a[par]

        for key in kwargs:
            self.plot_parameters[key] = deepcopy(kwargs[key])

    def __add__(self, other: Simulation) -> SimulationComposition:
        """
        The add operator is used to concatenate simulations.
        """
        return SimulationComposition(self, other)

    def to_dataframe(self) -> TypingAny:  # -> pandas.DataFrame
        """
        Convert simulation results to a pandas DataFrame.

        Returns:
            pandas.DataFrame: DataFrame containing the simulation results

        Raises:
            SimulationError: If no results are available
            ImportError: If pandas is not available
        """
        if not hasattr(self, "results") or self.results is None:
            raise SimulationError(
                "Simulation results were accessed before a simulation was executed"
            )

        try:
            return self.results.return_pandas()  # type: ignore[union-attr]
        except (AttributeError, TypeError, ValueError) as e:
            raise ImportError(f"Failed to convert results to DataFrame: {e!s}") from e

    @classmethod
    def is_simulation(cls) -> bool:
        """Return True to identify this class as a simulation."""
        return True

    @classmethod
    def set_job_number(cls, params: dict[str, TypingAny]) -> int:
        """Determine the joblib job count from simulation parameters."""
        try:
            if params["jobs"] == 1:  # noqa: SIM108
                jobs = params["jobs"]
            else:
                jobs = params["jobs"]
        except KeyError:
            jobs = -1
        return int(jobs)

    def __sub__(self, other: TypingAny) -> Simulation:
        return sof_sim_remove_reaction(self, other, Simulation)  # type: ignore[no-any-return,return-value]

    def __rsub__(self, other: TypingAny) -> Simulation:
        return sof_sim_remove_reaction(other, self, Simulation)  # type: ignore[no-any-return,return-value]


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
        S1: Simulation | SimulationComposition,
        S2: Simulation | SimulationComposition,
    ) -> None:
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
        self.results: MobsPyList_of_TS | dict[str, TypingAny] | None = None
        self.fres: MobsPyList_of_TS | dict[str, TypingAny] | None = None
        self.base_sim = self.list_of_simulations[0]

    def __add__(
        self,
        other: Simulation | SimulationComposition,
    ) -> SimulationComposition:
        return SimulationComposition(self, other)

    def __setattr__(self, name: str, value: TypingAny) -> None:
        white_list = ["list_of_simulations", "results", "base_sim", "fres"]
        multi_cast_parameters = ["duration"]
        broad_cast_parameters = ["level", "rate_type", "plot_type", "repetitions"]
        double_cast_parameters = ["simulation_method", "volume", "method"]

        if name in double_cast_parameters:
            # Broadcast if single value
            if isinstance(value, (str, int, float, Quantity)):
                for sim in self:
                    sim._set_parameter(name, value)
            else:
                # Multicast if list
                try:
                    value_len = len(value)
                except TypeError as e:
                    raise SimulationError(
                        f"The parameter {name} was assigned non-accepted type."
                    ) from e
                if value_len != len(self):
                    raise SimulationError(
                        f"The parameter {name} list length must match "
                        f"the number of simulations ({len(self)})."
                    )

                for par, sim in zip(value, self, strict=False):
                    # volume/duration changes after compilation
                    # are checked in setattr
                    if name == "volume":
                        sim.volume = par
                    elif name == "duration":
                        sim.duration = par
                    else:
                        sim._set_parameter(name, par)

        elif name in multi_cast_parameters:
            try:
                value_len = len(value)
            except TypeError as e:
                raise SimulationError(
                    "From 2.4.4 duration must be assigned to "
                    "each simulation individually or a list "
                    "with all durations must be assigned to "
                    "the concatenated simulation"
                ) from e
            if value_len != len(self):
                raise SimulationError(
                    f"The parameter {name} list length must match "
                    f"the number of simulations ({len(self)})."
                )

            for par, sim in zip(value, self, strict=False):
                # volume/duration changes after compilation
                # are checked in setattr
                if name == "volume":
                    sim.volume = par
                elif name == "duration":
                    sim.duration = par
                else:
                    sim._set_parameter(name, par)
        elif name in broad_cast_parameters:
            for sim in self:
                sim._set_parameter(name, value)
        elif name in white_list:
            self.__dict__[name] = value
        else:
            self.base_sim.__setattr__(name, value)

    @property
    def plot_config(self) -> PlotConfigProxy:
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

    # This run is for the multiple simulations
    def run(
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
        """
        Runs a concatenated simulation with multiple
        simulation objects summed. Some inputs are different
        here, duration and volume must receive any iterable
        with the volume for each simulation.

        Args:
            duration: Duration of a simulation.
            volume: Volume of the simulation, if none given 1 liter is used.
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

        for key in kwargs:
            self.base_sim.plot_parameters[key] = deepcopy(kwargs[key])

    def generate_sbml(self, compose: bool = False) -> list[str]:
        """
        Generates a string with an SBML model from a respective MobsPy model

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
        """
        Generates a string with an Antimony model from a respective MobsPy model

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

        return self.base_sim.generate_antimony(compose=compose, model_name=model_name)

    def to_dataframe(self) -> TypingAny:
        """Convert composition results to a pandas DataFrame."""
        self.base_sim.to_dataframe()

    @classmethod
    def is_simulation(cls) -> bool:
        """Return True to identify this class as a simulation."""
        return True
