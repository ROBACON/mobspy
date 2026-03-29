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

from mobspy.sbml_simulator.sbml_writer import create_model


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

    Args:
        species: Species as keys and counts as values.
        parameters: Parameter name and value.
        reactions: Reaction name and reaction in python sbml writer format.
        events: Event name and event in python sbml writer format.
        assignments: Assignments numbers and expressions.
        model_context: Optional unit context for proper SBML unit declarations.


    Returns:
        String describing the model in sbml format.
    """
    doc = create_model(
        species, parameters, reactions, events, assignments, model_context
    )

    # Convert sbml document into a string for basico
    sbml_str = sbml_writeSBMLToString(doc)
    return sbml_str  # type: ignore[no-any-return]
