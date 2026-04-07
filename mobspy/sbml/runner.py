"""Execute compiled SBML models via BasiCO/COPASI with parallel run support."""

from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING, Any

from joblib import Parallel, delayed

import mobspy.sbml.builder as sbml_builder
from mobspy.constants import END_FLAG_SPECIES_NAME
from mobspy.exceptions import SimulationError
from mobspy.lazy_import import LazyImporter as ipm_LazyImporter
from mobspy.mobspy_logging import get_logger
from mobspy.types import ConcreteSpeciesId, SBMLModelData

if TYPE_CHECKING:
    import pandas as pd

    from mobspy.types import CompiledModelDict, SimParams

_logger = get_logger(__name__)
basico = ipm_LazyImporter("basico")


def simulate(
    jobs: int,
    list_of_params: list[SimParams],
    models: list[CompiledModelDict],
) -> list[dict[str, list[float]]] | None:
    """Run SBML models via BasiCO and return time-series data."""
    return job_execution(list_of_params, models, jobs)


def job_execution(
    params: list[SimParams],
    models: list[CompiledModelDict],
    jobs: int,
) -> list[dict[str, list[float]]] | None:
    """Execute simulation repetitions in parallel via joblib."""

    def __single_run(packed: int) -> dict[str, list[float]]:
        i = packed

        added_data: dict[str, list[float]] = {}
        reformatted_data: dict[str, list[float]] = {}
        for j, (sim_par, model) in enumerate(zip(params, models, strict=True)):
            # Generate SBML here
            if j > 0:
                sbml_str = __sbml_new_initial_values(
                    reformatted_data, model, sim_par, new_model=True
                )
            else:
                sbml_str = __sbml_new_initial_values({}, model, sim_par)

            end_condition_not_satisfied = True
            if sim_par["_continuous_simulation"]:
                duration = float(sim_par["initial_conditional_duration"])
            else:
                duration = float(sim_par["duration"])

            while end_condition_not_satisfied:
                basico_model = basico.model_io.load_model_from_string(sbml_str)
                data = __run_time_course(basico_model, duration, sim_par, i)

                reformatted_data = reformat_time_series(data)

                if sim_par["_continuous_simulation"]:
                    if reformatted_data[END_FLAG_SPECIES_NAME][-1] > 0:
                        reformatted_data = __filter_condition_event_time_data(
                            reformatted_data
                        )
                        end_condition_not_satisfied = False
                    else:
                        sbml_str = __sbml_new_initial_values(
                            reformatted_data, model, sim_par
                        )
                        duration = 2 * duration
                else:
                    end_condition_not_satisfied = False

                reformatted_data = __remap_species(
                    reformatted_data, model.mappings, model.species_for_sbml
                )
                added_data = __add_simulations_data(added_data, reformatted_data)

        return added_data

    parallel_data: list[dict[str, list[float]]] | None = Parallel(n_jobs=jobs)(  # pyright: ignore[reportAssignmentType]  # joblib returns Any
        delayed(__single_run)(i) for i in range(params[0]["repetitions"])
    )

    if not parallel_data:
        raise SimulationError(
            "The parallel model has not produced an output. "
            "Try adding ('sequential': True) to parameters"
        )

    return parallel_data


def __run_time_course(
    basico_model: Any,
    duration: float,
    params: SimParams,
    index: int,
) -> pd.DataFrame:
    params["simulation_method"] = params["simulation_method"].lower()
    if (params["_with_event"] or params["_continuous_simulation"]) and params[
        "simulation_method"
    ] == "stochastic":
        params["simulation_method"] = "directmethod"

    kargs: dict[str, Any] = {
        "model": basico_model,
        "method": params["simulation_method"],
        "start_time": params["start_time"],
        "r_tol": params["r_tol"],
        "a_tol": params["a_tol"],
        "output_event": params["output_event"],
    }

    if "seeds" in params and params["seeds"] is not None:
        kargs["use_seed"] = True
        kargs["seed"] = params["seeds"][index]

    if "step_size" in params and params["step_size"] is not None:
        kargs["automatic"] = False
        kargs["step_number"] = int(params["duration"] / params["step_size"])

    return basico.run_time_course(duration, **kargs)


def reformat_time_series(
    data: pd.DataFrame,
) -> dict[str, list[float]]:
    """Convert a BasiCO DataFrame to a plain dict of lists."""
    data_dict: dict[str, list[float]] = {"Time": data.index.tolist()}

    for key in data:
        data_dict[ConcreteSpeciesId.from_sbml_id(key).to_display()] = list(data[key])

    return data_dict


def __filter_condition_event_time_data(
    data: dict[str, list[float]],
) -> dict[str, list[float]]:
    new_data: dict[str, list[float]] = {}

    stop_index = len(data[END_FLAG_SPECIES_NAME]) - 1
    for i, e in enumerate(data[END_FLAG_SPECIES_NAME]):
        if e == 1:
            stop_index = i
            break

    for key in data:
        new_data[key] = data[key][: stop_index + 1]

    return new_data


def __sbml_new_initial_values(
    data: dict[str, list[float]],
    model: CompiledModelDict,
    sim_para: SimParams,
    new_model: bool = False,
) -> str:
    species_for_sbml = model.species_for_sbml

    # BasiCO returns concentrations (amount/volume) when
    # hasOnlySubstanceUnits=False. Multiply by compartment volume
    # to recover amounts for species_for_sbml (used as initialAmount).
    _vol = 1.0
    _mc = getattr(model, "model_context", None)
    if _mc is not None:
        _vol = getattr(_mc, "resolved_volume_magnitude", 1.0)

    check_list = ["stochastic", "directmethod"]
    for key in data:
        parts = key.split(".")
        sbml_key = ConcreteSpeciesId(
            base=parts[0], characteristics=tuple(parts[1:])
        ).to_sbml_id()
        if sbml_key not in species_for_sbml:
            continue

        if key == "Time":
            continue
        try:
            final_val = list(data[key])[-1] * _vol
            if sim_para["simulation_method"].lower() in check_list:
                species_for_sbml[sbml_key] = int(final_val)
            else:
                species_for_sbml[sbml_key] = final_val
        except KeyError:
            _logger.debug(
                "Species '%s' not found in species_for_sbml during "
                "initial value update; skipping",
                sbml_key,
            )

    if new_model:
        with contextlib.suppress(KeyError):
            species_for_sbml[END_FLAG_SPECIES_NAME] = 0

    # Extract model_context if available (for proper SBML unit declarations)
    model_context = getattr(model, "model_context", None)

    return sbml_builder.build(
        SBMLModelData(
            species_for_sbml=species_for_sbml,
            parameters_for_sbml=model.parameters_for_sbml,
            reactions_for_sbml=model.reactions_for_sbml,
            events_for_sbml=model.events_for_sbml,
            assignments_for_sbml=model.assignments_for_sbml,
        ),
        model_context=model_context,
    )


def __add_simulations_data(
    added_data: dict[str, list[float]],
    reformatted_data: dict[str, list[float]],
) -> dict[str, list[float]]:
    time_to_add = added_data["Time"][-1] if added_data else 0
    new_data: dict[str, list[float]] = {}

    for key in added_data:  # noqa: PLC0206  # iterating dict keys directly
        new_data[key] = added_data[key]

    __adjust_time_offsets(added_data, reformatted_data, time_to_add)
    already_added_keys = __merge_existing_keys(new_data, added_data, reformatted_data)
    __merge_new_keys(
        new_data, added_data, reformatted_data, already_added_keys, time_to_add
    )

    if time_to_add != 0:
        new_data["Time"] = added_data["Time"] + reformatted_data["Time"]
    else:
        new_data["Time"] = reformatted_data["Time"]

    return new_data


def __adjust_time_offsets(
    added_data: dict[str, list[float]],
    reformatted_data: dict[str, list[float]],
    time_to_add: float,
) -> None:
    """Remove repeated initial values and offset time in reformatted data."""
    if added_data and reformatted_data.get("Time") and reformatted_data["Time"][0] == 0:
        for key in list(reformatted_data):
            reformatted_data[key].pop(0)
    for i in range(len(reformatted_data.get("Time", []))):
        reformatted_data["Time"][i] = reformatted_data["Time"][i] + time_to_add


def __merge_existing_keys(
    new_data: dict[str, list[float]],
    added_data: dict[str, list[float]],
    reformatted_data: dict[str, list[float]],
) -> set[str]:
    """Merge species keys that exist in both added_data and reformatted_data."""
    already_added_keys: set[str] = set()
    for key in added_data:  # noqa: PLC0206  # iterating dict keys directly
        if key == "Time":
            continue
        try:
            new_data[key] = added_data[key] + reformatted_data[key]
            already_added_keys.add(key)
        except KeyError:
            dummy = added_data[key][-1]
            new_data[key] += [dummy for _ in reformatted_data["Time"]]
    return already_added_keys


def __merge_new_keys(
    new_data: dict[str, list[float]],
    added_data: dict[str, list[float]],
    reformatted_data: dict[str, list[float]],
    already_added_keys: set[str],
    time_to_add: float,
) -> None:
    """Merge species keys that only exist in reformatted_data."""
    for key in reformatted_data:  # noqa: PLC0206  # iterating dict keys directly
        if key == "Time" or key in already_added_keys:
            continue
        if time_to_add != 0:
            new_data[key] = [0 for _ in added_data["Time"]]
            new_data[key] = new_data[key] + reformatted_data[key]
        else:
            new_data[key] = reformatted_data[key]


def __remap_species(
    data: dict[str, Any],
    mapping: dict[str, list[str]],
    species_not_mapped: dict[str, float],
) -> dict[str, Any]:
    mapped_data: dict[str, Any] = {"Time": data["Time"]}
    T = range(len(data["Time"]))  # noqa: N806  # legacy DSL variable convention

    # copy over all unmapped ones
    for k in data:  # noqa: PLC0206  # iterating dict keys directly
        mapped_data[k] = data[k]

    dot_species_not_mapped: dict[str, float] = {}
    for key in species_not_mapped:  # noqa: PLC0206  # iterating dict keys directly
        dot_key = ConcreteSpeciesId.from_sbml_id(key).to_display()
        dot_species_not_mapped[dot_key] = species_not_mapped[key]

    # 1st pass with sum mappings
    for group in mapping:  # noqa: PLC0206  # iterating dict keys directly
        the_mapping = mapping[group]
        mapped_data[group] = {"runs": []}

        try:
            # check if is a list -> sum
            if isinstance(the_mapping, list):
                this_run: list[float] = []
                runs_not_returned_by_basico: dict[str, list[float]] = {}
                for t in T:
                    mapping_sum: float = 0
                    for spe in the_mapping:
                        try:
                            mapping_sum = mapping_sum + data[spe][t]
                        except KeyError:
                            mapping_sum = mapping_sum + dot_species_not_mapped[spe]
                            try:
                                runs_not_returned_by_basico[spe] += [
                                    dot_species_not_mapped[spe]
                                ]
                            except KeyError:
                                runs_not_returned_by_basico[spe] = [
                                    dot_species_not_mapped[spe]
                                ]
                    this_run.append(mapping_sum)
                for spe in runs_not_returned_by_basico:  # noqa: PLC0206
                    try:
                        mapped_data[spe] = runs_not_returned_by_basico[spe]
                    except KeyError:
                        mapped_data[spe] = [runs_not_returned_by_basico[spe]]
                mapped_data[group] = this_run

        except IndexError as e:
            raise SimulationError(
                f'run: remap_species: error when remapping "{the_mapping}".'
                + "Possible fix: All runs must have the same time"
            ) from e

        except TypeError as e:
            _logger.warning(
                "Copasi removes A >> A species from"
                " reaction calculations and does"
                " not provide an output"
            )
            _logger.warning(
                "Please check the output data to see if this is the problem"
            )
            raise SimulationError(
                "TypeError while mapping simulation results. "
                "This may be caused by A >> A identity reactions."
            ) from e

    return mapped_data
