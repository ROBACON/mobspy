#!/usr/bin/env python3

"""
    Picks and translates a python model into an
    SBML readable by Copasi. New models must
    be references here as well.
"""


from mobspy.sbml_simulator.SBMLWriter import create_model
from libsbml import writeSBMLToString as sbml_writeSBMLToString


def build(species, parameters, reactions, events, assignments):
    """
    Constructs the sbml file for a model from the dictionary syntax for python sbml lib

    :param species: (dict) species as keys and counts as values
    :param parameters: (dict) parameter name and value
    :param reactions: (dict) reaction name and reaction in python sbml writer format
    :param events: (dict) event name and event in python sbml writer format
    :param assignments: (dict) assignments numbers and expressions

    :return: sbml_str (str) = string describing the model in sbml format
    """
    doc = create_model(species, parameters, reactions, events, assignments)

    # Convert sbml document into a string for basico
    sbml_str = sbml_writeSBMLToString(doc)
    return sbml_str




