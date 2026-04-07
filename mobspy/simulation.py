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
from pathlib import Path
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
    MobsPyTimeSeries,
    SimulationResults,
)
from mobspy.event_handling import EventHandlingMixin
from mobspy.exceptions import (
    CompilationError,
    MobsPyError,
    ParameterError,
    SimulationError,
    ValidationError,
)
from mobspy.execution import generate_sbml_from_compiled, run_sbml
from mobspy.mobspy_logging import get_logger
from mobspy.model_generation import ModelGenerationMixin
from mobspy.modules.any_species import (
    Any,
)
from mobspy.modules.assignments_implementation import (
    Assign,
)
from mobspy.modules.compiler import compile_model
from mobspy.modules.declarations import snapshot_registry
from mobspy.modules.list_species import List_Species
from mobspy.modules.logic_operators import (
    MetaSpeciesLogicResolver as lop_MetaSpeciesLogicResolver,
)
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
from mobspy.modules.rate_builder import (
    hill,
    where,
)
from mobspy.modules.set_counts_module import set_counts
from mobspy.modules.species import Species
from mobspy.modules.species_constructors import (
    BaseSpecies,
    ListSpecies,
    New,
    Zero,
)
from mobspy.modules.species_utils import (
    create_orthogonal_vector_structure as mcu_create_orthogonal_vector_structure,
)
from mobspy.modules.unit_handler import (
    extract_length_dimension as uh_extract_length_dimension,
)
from mobspy.modules.unit_registry import u
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
from mobspy.plotting import PlottingMixin, plot_results
from mobspy.simulation_composition import SimulationComposition
from mobspy.simulation_config import PlotConfig, SimulationConfig
from mobspy.simulator_object.utils import (
    sim_remove_reaction as sof_sim_remove_reaction,
)
from mobspy.types import (
    ConcreteModel,
    Delta,
    RateValue,
    SimulationEventData,
    SimulationMethod,
    SpeciesArg,
    TimeSeriesDataDict,
)

if TYPE_CHECKING:
    from mobspy.modules.reactions import Reactions
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
    "Delta",
    "ListSpecies",
    "ModelParameters",
    "New",
    "RateValue",
    "Rev",
    "Set",
    "Simulation",
    "SimulationComposition",
    "SimulationMethod",
    "SimulationResults",
    "SpeciesArg",
    "Zero",
    "basiCO_parameter_estimation",
    "compile_model",
    "generate_sbml_from_compiled",
    "hill",
    "logger",
    "plot_results",
    "run_sbml",
    "set_counts",
    "u",
    "where",
]

_logger = get_logger(__name__)
logger = _logger


class Simulation(
    EventHandlingMixin,
    ModelGenerationMixin,
    PlottingMixin,
):
    """Orchestrates compilation, execution, and result handling for a MobsPy model.

    Collects meta-species and reactions, compiles them via the Compiler into
    SBML, runs the simulation through BasiCO/COPASI, and stores the results.

    Uses mixins for event handling, model generation, and plotting.

    Examples:
        >>> from mobspy import *
        >>> A, B = BaseSpecies(['A', 'B'])
        >>> _ = A >> B @ 1
        >>> A(10)
        Species('A')
        >>> S = Simulation(A | B)
        >>> S.duration = 5
        >>> _ = S.compile(verbose=False)
        >>> S._is_compiled
        True
    """

    def __init__(  # noqa: PLR0913
        self,
        model: Species | List_Species,
        reactions: set[Reactions] | None = None,
        names: dict[str, TypingAny] | None = None,
        parameters: dict[str, TypingAny] | None = None,
        plot_parameters: dict[str, TypingAny] | None = None,
        backend: TypingAny = None,
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
            backend: Simulation backend (default: SBMLBackend).
                Must conform to the ``SimulationBackend`` protocol.

        Raises:
            ValidationError: If model contains invalid species types
            ParameterError: If required parameters are missing
        """
        if backend is None:
            from mobspy.backends import SBMLBackend  # noqa: PLC0415

            self._backend = SBMLBackend()
        else:
            self._backend = backend
        self.experimental_data: TypingAny = None
        self._declarations = snapshot_registry()
        self._init_event_state()
        self._init_model(model, names)
        self._init_reactions(reactions)
        self._init_counts()
        self._init_config(parameters, plot_parameters)
        self._init_sbml_state()

    def _init_event_state(self) -> None:
        """Initialize event tracking and compilation state."""
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
        self.model_parameters: dict[str, TypingAny] = {}
        self.sbml_data_list: ParameterSweepList = []
        self._parameter_list_of_dic: list[dict[str, TypingAny]] = []
        self._is_compiled = False
        self.dimension: int | None = None
        self.model_parameter_objects_dict: dict[str, TypingAny] | None = None

    def _init_model(
        self,
        model: Species | List_Species,
        names: dict[str, TypingAny] | None,
    ) -> None:
        """Validate and expand the model with linked species."""
        if not isinstance(model, (Species, List_Species)):
            raise ValidationError(
                "Model must be formed only by Species objects "
                "or List_Species objects. "
                f"Received type {type(model)} with value {model}"
            )

        model_pre_link = List_Species(model)  # type: ignore[arg-type]
        model_pos_link: set[Species] = set()
        for spe in model_pre_link:
            model_pos_link.add(spe)
            model_pos_link = model_pos_link.union(spe._linked_species)
        self.model = List_Species(model_pos_link)
        self.names = names
        self.orthogonal_vector_structure = mcu_create_orthogonal_vector_structure(model)  # type: ignore[arg-type]

    def _init_reactions(self, reactions: set[Reactions] | None) -> None:
        """Collect reactions from explicit set or from the model registry.

        Filters the registry snapshot by species membership: only reactions
        where at least one participant is in the model (including references).
        """
        if reactions is not None:
            self._reactions_set = set(reactions)
        else:
            model_species_ids = self._model_species_ids()
            self._reactions_set = self._declarations.reactions_for_species(
                model_species_ids
            )

    def _model_species_ids(self) -> frozenset[int]:
        """Collect ids of all species in the model, including references."""
        ids: set[int] = set()
        for spe_object in self.model:
            for reference in spe_object.get_references():
                ids.add(id(reference))
        return frozenset(ids)

    def _init_counts(self) -> None:
        """Gather species counts from the registry or Species objects.

        The registry is the primary source (supports overwrite and
        reset semantics). Falls back to Species traversal when
        the registry has no counts for model species.
        """
        model_species = set(self.model)
        registry_counts = [
            ca for ca in self._declarations.counts if ca.species in model_species
        ]
        if registry_counts:
            self._species_counts = [
                {
                    "object": ca.species,
                    "characteristics": ca.characteristics,
                    "quantity": ca.quantity,
                }
                for ca in registry_counts
            ]
        else:
            # Legacy fallback
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

    def _init_config(
        self,
        parameters: dict[str, TypingAny] | None,
        plot_parameters: dict[str, TypingAny] | None,
    ) -> None:
        """Set simulation and plot configuration."""
        if not parameters:
            self.parameters = SimulationConfig()  # type: ignore[assignment]
        else:
            config = SimulationConfig()
            config.update(parameters)
            self.parameters = config  # type: ignore[assignment]

        if not plot_parameters:
            self.plot_parameters: dict[str, TypingAny] = PlotConfig()
        else:
            self.plot_parameters = PlotConfig(plot_parameters)

        self.results: SimulationResults | dict[str, TypingAny] = {}
        self.fres: SimulationResults | dict[str, TypingAny] = {}
        self.default_order = Default

    def _init_sbml_state(self) -> None:
        """Initialize SBML compilation output slots."""
        self._concrete_model: ConcreteModel | None = None
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

    # ------------------------------------------------------------------
    # Inlined from Experimental_Data_Holder
    # ------------------------------------------------------------------

    def load_experiment_data(self, data: TypingAny) -> None:
        """Store experimental data for parameter estimation.

        Args:
            data: List of dicts or a SimulationResults result.

        Raises:
            ValidationError: If the data format is invalid.
        """
        flag_jump_checks = isinstance(data, SimulationResults)
        if not isinstance(data, list) and not flag_jump_checks:
            raise ValidationError(
                "Data added must be in the format of list with"
                " each element being a dictionary "
                "with species names and time as keys"
                " or a MobsPy results object"
            )
        for e in data:
            if not isinstance(e, dict) and not flag_jump_checks:
                raise ValidationError(
                    "Data added must be in the format of list"
                    " with each element being a dictionary "
                    "with species names and time as keys"
                    " or a MobsPy results object"
                )
        self.experimental_data = data

    # ------------------------------------------------------------------
    # Inlined from Simulation_Utils
    # ------------------------------------------------------------------

    def update_model(self, *args: TypingAny) -> None:
        """Update species counts or parameters on an already-compiled model.

        Args:
            *args: Pairs of ``(name, value)`` to update.

        Raises:
            SimulationError: If the model has not been compiled yet
                or arguments are malformed.
        """
        if not self._list_of_models:
            raise SimulationError(
                "In .update_model method - \n"
                "The model was not compiled yet. The update_model"
                " method is reserved for simulations that "
                "have already been compiled"
            )
        _NAME_VALUE_PAIR_LEN = 2  # noqa: N806
        for arg in args:
            if len(arg) != _NAME_VALUE_PAIR_LEN:
                raise SimulationError(
                    "In .update_model method - \n"
                    "Please all parameters and species changes"
                    " must be in the format: \n"
                    "(name, value)"
                )
            self._update_from_compiler(arg)

    def _update_from_compiler(self, arg: TypingAny) -> None:
        """Dispatch a (name, value) update to either parameters or species."""
        from mobspy.modules.mobspy_parameters import (  # noqa: PLC0415
            Internal_Parameter_Constructor,
        )
        from mobspy.types import ConcreteSpeciesId as _CID  # noqa: PLC0415, N814

        try:
            is_species = arg[0].is_spe_or_reac()
        except AttributeError:
            is_species = False

        if isinstance(arg[0], Internal_Parameter_Constructor):
            self._update_parameter(arg)
        elif isinstance(arg[0], str):
            test_model = self._list_of_models[0]
            not_parameter = arg[0] not in test_model.parameters_for_sbml
            display_key = _CID.from_sbml_id(arg[0]).to_display()
            not_species = display_key not in test_model.species_for_sbml

            if not not_parameter:
                self._update_parameter(arg)
            if not not_species:
                self._update_species(arg)
            if not_species and not_parameter:
                raise SimulationError(
                    f"The string {arg[0]} was not found either in parameters or species"
                )
        elif is_species:
            self._update_species(arg)
        else:
            raise SimulationError("Unsupported argument type for model update")

    def _update_parameter(self, arg: TypingAny) -> None:
        """Update a parameter value across all compiled models."""
        try:
            iterable = iter(arg[1])
        except TypeError:
            iterable = False  # type: ignore[assignment]

        value_to_update = arg[1][0] if iterable else arg[1]
        parameter_str = arg[0] if isinstance(arg[0], str) else arg[0].get_name()

        for model in self._list_of_models:
            try:
                model.parameters_for_sbml[parameter_str] = (
                    value_to_update,
                    "dimensionless",
                )
            except KeyError as e:
                raise SimulationError(
                    f"The parameter named {parameter_str} was not found in the model"
                ) from e

        parameter_object = self.model_parameters[parameter_str].object
        parameter_object.update_value(arg[1])

        try:
            if not iterable:
                self.model_parameters[parameter_str].values = [parameter_object.value]
            else:
                self.model_parameters[parameter_str].values = parameter_object.value
        except KeyError as e:
            raise SimulationError(
                f"The parameter named {parameter_str} was not found in the model"
            ) from e

    def _update_species(self, arg: TypingAny) -> None:
        """Update species counts in the compiled model."""
        from mobspy.constants import ALL_CHAR as _ALL  # noqa: PLC0415
        from mobspy.modules.species_string_generator import (  # noqa: PLC0415
            construct_all_species_ids as sp_construct_all_species_ids,
        )
        from mobspy.modules.species_string_generator import (  # noqa: PLC0415
            construct_species_id as sp_construct_species_id,
        )
        from mobspy.modules.unit_handler import (  # noqa: PLC0415
            convert_counts as uh_convert_counts_fn,
        )

        volume = self.__dict__.get("volume", 1)
        dimension = self.__dict__["dimension"]
        model_context = getattr(self, "_model_context", None)

        spe_count = uh_convert_counts_fn(
            arg[1], volume, dimension, model_context=model_context
        )

        query = arg[0].get_query_characteristics()
        if _ALL in query:
            spe_ids = sp_construct_all_species_ids(
                arg[0], query, self.orthogonal_vector_structure
            )
            for sid in spe_ids:
                self._list_of_models[0].species_for_sbml[sid.to_sbml_id()] = spe_count
        else:
            spe_string = sp_construct_species_id(
                arg[0],
                query,
                self.orthogonal_vector_structure,
            ).to_sbml_id()
            self._list_of_models[0].species_for_sbml[spe_string] = spe_count

    def compile(self, verbose: bool = True) -> str | None:
        """
        Compile the chemical reaction network model into executable SBML format.

        This method processes the meta-species, reactions, and parameters to generate
        an SBML representation that can be executed by the simulation engine.

        Args:
            verbose (bool): If True, print detailed compilation
                information. Defaults to True.

        Returns:
            The compiled model summary string, or None if no reactions
            or species were defined (empty model).

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

            # MobsPy level: 0=errors, 1=+warnings, 2=+info, 3=+debug
            _LEVEL_MAP = {  # noqa: N806
                0: logging.ERROR,
                1: logging.WARNING,
                2: logging.INFO,
                3: logging.DEBUG,
            }
            mobspy_level = self.parameters.get("level", 3)
            log_level = _LEVEL_MAP.get(mobspy_level, logging.INFO)
            _logger.set_log_level(log_level)

            # Resolve model unit context BEFORE parameter processing
            # (parameter_process strips units from duration/volume)
            _model_context = ModelUnitContext.from_simulation(
                volume=self.parameters["volume"],
                duration=self.parameters["duration"],
                species_counts=self._species_counts,
                dimension=self.dimension,
            )

            pr_parameter_process(self.parameters, model_context=_model_context)  # type: ignore[arg-type]
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

            self.parameters["_end_condition"] = self._end_condition

            _result = compile_model(
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

        # Store the backend-agnostic IR
        self._concrete_model = _result.to_concrete_model()

        # Extract fields for backward compatibility
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

        compiled_model = self._concrete_model.to_compiled_model(
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

    def _process_run_parameters(self, **kwargs: TypingAny) -> None:
        """Process and apply run-time parameter overrides, then ensure compilation."""
        pr_manually_process_each_parameter(self, **kwargs)

        level = kwargs.get("level")
        if level is not None:
            self.level = level

        # Base case - If there are no events we compile the model here
        if self._species_for_sbml is None:
            self.compile(verbose=False)

        self._assemble_multi_simulation_structure()

    def _execute_simulations(self) -> tuple[list[TypingAny], int]:
        """Run simulations via the pluggable backend."""
        jobs = self.set_job_number(self.parameters)  # type: ignore[arg-type]
        results = self._backend.run(
            self.sbml_data_list, self._list_of_parameters, jobs=jobs
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
            flatt_ts = [
                (ts, params)
                for r, params in zip(
                    raw_results, self._parameter_list_of_dic, strict=False
                )
                for ts in r  # pyright: ignore[reportOptionalIterable]
            ]
        else:
            flatt_ts = [
                (ts, {})
                for r in raw_results
                for ts in r  # pyright: ignore[reportOptionalIterable]
            ]

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

        self.results = SimulationResults(
            all_processed_data,  # pyright: ignore[reportArgumentType]
            self.model_parameter_objects_dict,  # pyright: ignore[reportArgumentType]
        )
        self.fres = SimulationResults([all_processed_data[0]], None, True)  # pyright: ignore[reportArgumentType, reportIndexIssue]

    def run(  # noqa: PLR0913
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
    ) -> SimulationResults:
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
            elif "stochastic" in methods_list:
                self.plot_stochastic()
            else:
                self.plot_deterministic()
        return self.results  # type: ignore[return-value]

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
                with Path(out_path).open("w", encoding="utf-8") as f:
                    json_dump(self.results.to_dict(), f, indent=4)  # type: ignore[union-attr]
            else:
                # Add .json extension if not present
                if not file.endswith(".json"):
                    file += ".json"
                with Path(file).open("w", encoding="utf-8") as jf:
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
        with Path(file_name).open(encoding="utf-8") as json_file:
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
            "_declarations",
            "_backend",
            "_concrete_model",
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
    def plot_config(self) -> PlotConfig:
        """Access plot configuration. Usage: sim.plot_config.param = value."""
        result: PlotConfig = self.__dict__["plot_parameters"]
        return result

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
            if os_path_splitext(config)[1] != ".json":  # noqa: PTH122
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

        for key in kwargs:  # noqa: PLC0206
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
            jobs = params["jobs"]
        except KeyError:
            jobs = -1
        return int(jobs)

    def __sub__(self, other: TypingAny) -> Simulation:
        return sof_sim_remove_reaction(self, other, Simulation)  # type: ignore[no-any-return,return-value]

    def __rsub__(self, other: TypingAny) -> Simulation:
        return sof_sim_remove_reaction(other, self, Simulation)  # type: ignore[no-any-return,return-value]
