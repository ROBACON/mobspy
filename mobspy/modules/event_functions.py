"""Build and validate SBML events from user-defined triggers and assignments."""

from __future__ import annotations

from re import split as re_split
from typing import TYPE_CHECKING, Any

from mobspy.constants import ALL_CHAR
from mobspy.exceptions import EventError
from mobspy.modules.unit_handler import convert_counts as uh_convert_counts
from mobspy.modules.unit_handler import convert_time as uh_convert_time
from mobspy.types import CompilationContext, Delta, EventData, SimulationEventData

if TYPE_CHECKING:
    from mobspy.types import EventsForSbml

from pint import Quantity

from mobspy.modules.mobspy_parameters import (
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)
from mobspy.modules.species_string_generator import (
    construct_all_species_ids as ssg_construct_all_species_ids,
)
from mobspy.modules.species_string_generator import (
    construct_species_id as ssg_construct_species_id,
)


def format_event_dictionary_for_sbml(  # noqa: PLR0913
    species_for_sbml: dict[str, int | float],
    event_list: list[SimulationEventData],
    characteristics_to_object: dict[str, Any],
    meta_species_to_simulate: Any,
    parameters_in_events: set[mp_Mobspy_Parameter],
    ctx: CompilationContext,
) -> tuple[EventsForSbml, set[str]]:
    """
    Creates events_for_sbml dictionary for sbml file construction
    on the sbml_simulator/sbml_writer.py

    Args:
        species_for_sbml: {'species_string': count, ....}.
        event_list:  [{'species':'meta_species_object', 'characteristics':['list of
            characteristics'], 'quantity': number or pint object}.
        characteristics_to_object: {'characteristic':'meta_species_object', .....}.
        volume: Simulation volume for unit conversion.
        dimension: Dimension of the system 2D, 3D, 4D, e.


    Returns:
        Event dictionary for the sbml file construction.
    """
    species_in_events: set[str] = set()

    _convert_trigger_counts(event_list, ctx)

    reformed_event_list = _build_reformed_event_list(
        event_list,
        characteristics_to_object,
        meta_species_to_simulate,
        parameters_in_events,
        ctx,
    )

    events_for_sbml = _assemble_events_for_sbml(
        reformed_event_list,
        species_for_sbml,
        species_in_events,
        parameters_in_events,
    )

    return events_for_sbml, species_in_events


def _convert_trigger_counts(
    event_list: list[SimulationEventData],
    ctx: CompilationContext,
) -> None:
    """Convert Quantity values inside event triggers to model units."""
    for ev in event_list:
        if ev.trigger != "true":
            if isinstance(ev.trigger, str):
                raise EventError(
                    "Event trigger must be a condition object, "
                    f"got string: {ev.trigger!r}"
                )
            for i, e in enumerate(ev.trigger.operation):
                if isinstance(e, Quantity):
                    ev.trigger.operation[i] = uh_convert_counts(
                        e,
                        ctx.volume,
                        ctx.dimension,
                        model_context=ctx.model_context,
                    )


def _build_reformed_event_list(
    event_list: list[SimulationEventData],
    characteristics_to_object: dict[str, Any],
    meta_species_to_simulate: Any,
    parameters_in_events: set[mp_Mobspy_Parameter],
    ctx: CompilationContext,
) -> list[SimulationEventData]:
    """Process events into reformed event list with resolved species keys."""
    reformed_event_list: list[SimulationEventData] = []

    for ev in event_list:
        if not ev.event_counts:
            continue

        converted_time: Any = ev.event_time
        if not isinstance(ev.event_time, mp_Mobspy_Parameter):
            result = uh_convert_time(ev.event_time, model_context=ctx.model_context)
            if result is not None:
                converted_time = result

        event_dictionary: dict[str, int | float | str] = {}

        _process_all_char_assignments(
            ev,
            event_dictionary,
            characteristics_to_object,
            ctx,
        )
        _process_specific_assignments(
            ev,
            event_dictionary,
            characteristics_to_object,
            parameters_in_events,
            ctx,
        )

        if isinstance(ev.trigger, str):
            reformed_event_list.append(
                SimulationEventData(
                    event_time=converted_time,
                    event_counts=event_dictionary,  # type: ignore[arg-type]
                    trigger=ev.trigger,
                )
            )
        else:
            for e in ev.trigger.operation:  # pyright: ignore[reportGeneralTypeIssues]
                if isinstance(e, dict) and e["object"] not in meta_species_to_simulate:
                    raise EventError(
                        f"Meta species {e['object']} was used"
                        " in an event but is not in the model"
                    )
            reformed_event_list.append(
                SimulationEventData(
                    event_time=converted_time,
                    event_counts=event_dictionary,  # type: ignore[arg-type]
                    trigger=ev.trigger.generate_string(
                        characteristics_to_object, to_sort=True
                    ),
                )
            )

    return reformed_event_list


def _process_all_char_assignments(
    ev: SimulationEventData,
    event_dictionary: dict[str, int | float | str],
    characteristics_to_object: dict[str, Any],
    ctx: CompilationContext,
) -> None:
    """Process event counts that use the ALL_CHAR wildcard."""
    for ec in ev.event_counts:
        if ALL_CHAR not in ec["characteristics"]:
            continue
        temp_char = set(ec["characteristics"])
        temp_char.remove(ALL_CHAR)
        species_ids = ssg_construct_all_species_ids(
            ec["species"],
            temp_char,
            characteristics_to_object,
        )
        for sid in species_ids:
            d = sid.to_sbml_id()
            if isinstance(ec["quantity"], Delta):
                raw = ec["quantity"].value
                converted = uh_convert_counts(
                    raw,
                    ctx.volume,
                    ctx.dimension,
                    model_context=ctx.model_context,
                )
                event_dictionary[d] = f"{d} + {converted}"
            elif not isinstance(ec["quantity"], str):
                event_dictionary[d] = uh_convert_counts(
                    ec["quantity"],
                    ctx.volume,
                    ctx.dimension,
                    model_context=ctx.model_context,
                )
            else:
                event_dictionary[d] = ec["quantity"]


def _process_specific_assignments(
    ev: SimulationEventData,
    event_dictionary: dict[str, int | float | str],
    characteristics_to_object: dict[str, Any],
    parameters_in_events: set[mp_Mobspy_Parameter],
    ctx: CompilationContext,
) -> None:
    """Process event counts that target specific characteristics."""
    for ec in ev.event_counts:
        if ALL_CHAR in ec["characteristics"]:
            continue
        dummy_key = ssg_construct_species_id(
            ec["species"],
            ec["characteristics"],
            characteristics_to_object,
        ).to_sbml_id()
        if isinstance(ec["quantity"], Delta):
            raw = ec["quantity"].value
            converted = uh_convert_counts(
                raw,
                ctx.volume,
                ctx.dimension,
                model_context=ctx.model_context,
            )
            event_dictionary[dummy_key] = f"{dummy_key} + {converted}"
        elif not isinstance(ec["quantity"], str):
            if isinstance(ec["quantity"], mp_Mobspy_Parameter):
                parameters_in_events.add(ec["quantity"])
                event_dictionary[dummy_key] = ec["quantity"].name
            else:
                event_dictionary[dummy_key] = uh_convert_counts(
                    ec["quantity"],
                    ctx.volume,
                    ctx.dimension,
                    model_context=ctx.model_context,
                )
        else:
            if ctx.parameter_exist:
                for token in re_split(r", |-|!|\*|\+|/|\)|\(| ", ec["quantity"]):
                    name = token.strip()
                    if name and name in ctx.parameter_exist:
                        parameters_in_events.add(ctx.parameter_exist[name])
            event_dictionary[dummy_key] = ec["quantity"]


def _assemble_events_for_sbml(
    reformed_event_list: list[SimulationEventData],
    species_for_sbml: dict[str, int | float],
    species_in_events: set[str],
    parameters_in_events: set[mp_Mobspy_Parameter],
) -> dict[str, EventData]:
    """Convert reformed events into the final events_for_sbml dict."""
    events_for_sbml: dict[str, EventData] = {}
    for i, event in enumerate(reformed_event_list):
        assignments: list[tuple[str, str | int | float]] = []
        for key in event.event_counts:
            if key in species_for_sbml:
                assignments.append((key, str(event.event_counts[key])))
                species_in_events.add(key)
            else:
                raise EventError(
                    f"Species '{key}' used in an event assignment "
                    "but it is not in the compiled model. "
                    "If using inheritance, list child species "
                    "explicitly in Simulation(), e.g. "
                    "Simulation(ChildA | ChildB) instead of "
                    "Simulation(Parent)."
                )

        assignments.sort()

        if isinstance(event.event_time, mp_Mobspy_Parameter):
            for par in event.event_time._parameter_set:
                parameters_in_events.add(par)

        events_for_sbml["e" + str(i)] = EventData(
            trigger=event.trigger,
            delay=str(event.event_time),
            assignments=assignments,
        )

    return events_for_sbml
