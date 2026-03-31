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

from mobspy.sbml_simulator.sbml_writer import create_model
from mobspy.types import SBMLModelData


def build(
    model_data: SBMLModelData,
    *,
    model_context: ModelUnitContext | None = None,
) -> str:
    """Construct an SBML string from model data.

    Args:
        model_data: Species, parameters, reactions, events, and assignments.
        model_context: Optional unit context for proper SBML unit declarations.

    Returns:
        String describing the model in sbml format.
    """
    doc = create_model(
        model_data.species_for_sbml or None,
        model_data.parameters_for_sbml or None,
        model_data.reactions_for_sbml or None,
        model_data.events_for_sbml or None,
        model_data.assignments_for_sbml or None,
        model_context=model_context,
    )

    # Convert sbml document into a string for basico
    result: str = sbml_writeSBMLToString(doc)
    return result
