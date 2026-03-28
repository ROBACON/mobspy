from __future__ import annotations

from typing import TYPE_CHECKING, Any

from mobspy.exceptions import EventError
from mobspy.modules.unit_handler import convert_counts as uh_convert_counts
from mobspy.types import EventData, SimulationEventData

if TYPE_CHECKING:
    from mobspy.types import EventsForSbml

from pint import Quantity  # noqa: E402

from mobspy.modules.function_rate_code import (  # noqa: E402
    search_for_parameters_in_str as frc_search_for_parameters_in_str,
)
from mobspy.modules.mobspy_parameters import (  # noqa: E402
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)
from mobspy.modules.species_string_generator import (  # noqa: E402
    construct_all_combinations as ssg_construct_all_combinations,
)
from mobspy.modules.species_string_generator import (  # noqa: E402
    construct_species_char_list as ssg_construct_species_char_list,
)

# @TODO remove search parameters in string
# don't use it anymore - slowly deprecate this function


def format_event_dictionary_for_sbml(
    species_for_sbml: dict[str, Any],
    event_list: list[SimulationEventData],
    characteristics_to_object: dict[str, Any],
    volume: int | float,
    dimension: int,
    meta_species_to_simulate: Any,
    parameter_exist: dict[str, Any],
    parameters_in_events: set[Any],
) -> tuple[EventsForSbml, set[str]]:
    """
    Creates events_for_sbml dictionary for sbml file construction
    on the sbml_simulator/SBMLWriter.py

    :param species_for_sbml: (dict)
        {'species_string': count, ....}
    :param event_list:
        [{'species':'meta_species_object',
        'characteristics':['list of characteristics'],
        'quantity': number or pint object}
    :param characteristics_to_object: {'characteristic':'meta_species_object', .....}
    :param volume: Simulation volume for unit conversion
    :param dimension: Dimension of the system 2D, 3D, 4D, e
    :return: event dictionary for the sbml file construction
    :rtype: events = {'e': { 'trigger': 'true',
        'delay': '10', 'assignments': [('M','1'),]}}
    """
    reformed_event_list: list[SimulationEventData] = []
    species_in_events: set[str] = set()

    # Convert count from triggers
    for ev in event_list:
        if ev.trigger != "true":
            for i, e in enumerate(ev.trigger.operation):
                if isinstance(e, Quantity):
                    ev.trigger.operation[i] = uh_convert_counts(e, volume, dimension)

    for ev in event_list:
        if not ev.event_counts:
            continue
        event_dictionary: dict[str, Any] = {}

        # All assignments never take priority over specific assignments
        for ec in ev.event_counts:
            if "all$" not in ec["characteristics"]:
                continue

            temp_char = set(ec["characteristics"])
            temp_char.remove("all$")
            dummy = ssg_construct_all_combinations(
                ec["species"], temp_char, characteristics_to_object, symbol="_dot_"
            )
            for d in dummy:
                if type(ec["quantity"]) != str:  # noqa: E721
                    event_dictionary[d] = uh_convert_counts(
                        ec["quantity"], volume, dimension
                    )
                else:
                    event_dictionary[d] = ec["quantity"]

        for ec in ev.event_counts:
            if "all$" in ec["characteristics"]:
                continue

            dummy = ssg_construct_species_char_list(
                ec["species"],
                ec["characteristics"],
                characteristics_to_object,
                symbol="_dot_",
            )

            if type(ec["quantity"]) != str:  # noqa: E721
                if isinstance(ec["quantity"], mp_Mobspy_Parameter):
                    parameters_in_events.add(ec["quantity"])
                    event_dictionary[dummy] = ec["quantity"].name

                event_dictionary[dummy] = uh_convert_counts(
                    ec["quantity"], volume, dimension
                )
            else:
                if parameter_exist != {}:
                    frc_search_for_parameters_in_str(
                        ec["quantity"], parameter_exist, parameters_in_events
                    )
                event_dictionary[dummy] = ec["quantity"]

        if type(ev.trigger) == str:  # noqa: E721
            reformed_event_list.append(
                SimulationEventData(
                    event_time=ev.event_time,
                    event_counts=event_dictionary,
                    trigger=ev.trigger,
                )
            )
        else:
            for e in ev.trigger.operation:
                if type(e) == dict:  # noqa: SIM102, E721
                    if e["object"] not in meta_species_to_simulate:
                        raise EventError(
                            f"Meta species {e['object']} was used"
                            " in an event but is not in the model"
                        )
            reformed_event_list.append(
                SimulationEventData(
                    event_time=ev.event_time,
                    event_counts=event_dictionary,
                    trigger=ev.trigger.generate_string(
                        characteristics_to_object, to_sort=True
                    ),
                )
            )

    events_for_sbml: dict[str, EventData] = {}
    for i, event in enumerate(reformed_event_list):
        assignments: list[tuple[str, str]] = []
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
