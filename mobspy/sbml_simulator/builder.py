#!/usr/bin/env python3
"""
Picks and translates a python model into an
SBML readable by Copasi. New models must
be references here as well.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from libsbml import writeSBMLToString as sbml_writeSBMLToString

if TYPE_CHECKING:
    from mobspy.modules.model_unit_context import ModelUnitContext
    from mobspy.types import (
        AssignmentsForSbml,
        EventsForSbml,
        ParametersForSbml,
        ReactionsForSbml,
        SpeciesForSbml,
    )

from mobspy.sbml_simulator.SBMLWriter import create_model


def build(
    species: SpeciesForSbml | None,
    parameters: ParametersForSbml | None,
    reactions: ReactionsForSbml | None,
    events: EventsForSbml | None,
    assignments: AssignmentsForSbml | None,
    model_context: ModelUnitContext | None = None,
) -> str:
    """
    Constructs the sbml file for a model from the dictionary syntax for python sbml lib

    :param species: (dict) species as keys and counts as values
    :param parameters: (dict) parameter name and value
    :param reactions: (dict) reaction name and reaction in python sbml writer format
    :param events: (dict) event name and event in python sbml writer format
    :param assignments: (dict) assignments numbers and expressions
    :param model_context: optional unit context for proper SBML unit declarations

    :return: sbml_str (str) = string describing the model in sbml format
    """
    doc = create_model(species, parameters, reactions, events, assignments, model_context)

    # Convert sbml document into a string for basico
    sbml_str = sbml_writeSBMLToString(doc)
    return sbml_str
