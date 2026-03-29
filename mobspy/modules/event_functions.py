"""Build and validate SBML events from user-defined triggers and assignments."""

from __future__ import annotations

from re import split as re_split
from typing import TYPE_CHECKING, Any

from mobspy.constants import ALL_CHAR, DOT_SEPARATOR
from mobspy.exceptions import EventError
from mobspy.modules.unit_handler import convert_counts as uh_convert_counts
from mobspy.types import EventData, SimulationEventData

if TYPE_CHECKING:
    from mobspy.modules.model_unit_context import ModelUnitContext
    from mobspy.types import EventsForSbml

from pint import Quantity

from mobspy.modules.mobspy_parameters import (
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)
from mobspy.modules.species_string_generator import (
    construct_all_combinations as ssg_construct_all_combinations,
)
from mobspy.modules.species_string_generator import (
    construct_species_char_list as ssg_construct_species_char_list,
)


def format_event_dictionary_for_sbml(
    species_for_sbml: dict[str, int | float],
    event_list: list[SimulationEventData],
    characteristics_to_object: dict[str, Any],
    volume: int | float,
    dimension: int,
    meta_species_to_simulate: Any,
    parameter_exist: dict[str, Any],
    parameters_in_events: set[mp_Mobspy_Parameter],
    model_context: ModelUnitContext | None = None,
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
    reformed_event_list: list[SimulationEventData] = []
    species_in_events: set[str] = set()

    # Convert count from triggers
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
                        volume,
                        dimension,
                        model_context=model_context,
                    )

    for ev in event_list:
        if not ev.event_counts:
            continue
        event_dictionary: dict[str, int | float | str] = {}

        # All assignments never take priority over specific assignments
        for ec in ev.event_counts:
            if ALL_CHAR not in ec["characteristics"]:
                continue

            temp_char = set(ec["characteristics"])
            temp_char.remove(ALL_CHAR)
            dummy = ssg_construct_all_combinations(
                ec["species"],
                temp_char,
                characteristics_to_object,
                symbol=DOT_SEPARATOR,
            )
            for d in dummy:
                if not isinstance(ec["quantity"], str):
                    event_dictionary[d] = uh_convert_counts(
                        ec["quantity"],
                        volume,
                        dimension,
                        model_context=model_context,
                    )
                else:
                    event_dictionary[d] = ec["quantity"]

        for ec in ev.event_counts:
            if ALL_CHAR in ec["characteristics"]:
                continue

            dummy_result = ssg_construct_species_char_list(
                ec["species"],
                ec["characteristics"],
                characteristics_to_object,
                symbol=DOT_SEPARATOR,
            )
            dummy_key: str = (
                dummy_result if isinstance(dummy_result, str) else str(dummy_result)
            )

            if not isinstance(ec["quantity"], str):
                if isinstance(ec["quantity"], mp_Mobspy_Parameter):
                    parameters_in_events.add(ec["quantity"])
                    event_dictionary[dummy_key] = ec["quantity"].name
                else:
                    event_dictionary[dummy_key] = uh_convert_counts(
                        ec["quantity"],
                        volume,
                        dimension,
                        model_context=model_context,
                    )
            else:
                if parameter_exist:
                    for token in re_split(r", |-|!|\*|\+|/|\)|\(| ", ec["quantity"]):
                        name = token.strip()
                        if name and name in parameter_exist:
                            parameters_in_events.add(parameter_exist[name])
                event_dictionary[dummy_key] = ec["quantity"]

        if isinstance(ev.trigger, str):
            reformed_event_list.append(
                SimulationEventData(
                    event_time=ev.event_time,
                    event_counts=event_dictionary,  # type: ignore[arg-type]
                    trigger=ev.trigger,
                )
            )
        else:
            for e in ev.trigger.operation:  # pyright: ignore[reportGeneralTypeIssues]
                if isinstance(e, dict):  # noqa: SIM102
                    if e["object"] not in meta_species_to_simulate:
                        raise EventError(
                            f"Meta species {e['object']} was used"
                            " in an event but is not in the model"
                        )
            reformed_event_list.append(
                SimulationEventData(
                    event_time=ev.event_time,
                    event_counts=event_dictionary,  # type: ignore[arg-type]
                    trigger=ev.trigger.generate_string(
                        characteristics_to_object, to_sort=True
                    ),
                )
            )

    events_for_sbml: dict[str, EventData] = {}
    for i, event in enumerate(reformed_event_list):
        assignments: list[tuple[str, str | int | float]] = []
        for key in event.event_counts:
            if key in species_for_sbml:
                assignments.append((key, str(event.event_counts[key])))
                species_in_events.add(key)
            else:
                raise EventError(
                    f"Species {key} used in an event assignment"
                    " but it is not in the model"
                )

        assignments.sort()

        if event.event_time:
            pass

        if isinstance(event.event_time, mp_Mobspy_Parameter):
            for par in event.event_time._parameter_set:
                parameters_in_events.add(par)

        events_for_sbml["e" + str(i)] = EventData(
            trigger=event.trigger,
            delay=str(event.event_time),
            assignments=assignments,
        )

    return events_for_sbml, species_in_events
