#!/usr/bin/env python3

#####################################################
#                                                   #
#   Picks and translates a python model into an     #
#   SBML readable by Copasi. New models must        #
#   be references here as well.                     #
#                                                   #
#####################################################


from mobspy.sbml_simulator.SBMLWriter import create_model
import libsbml as sbml


def build(species, parameters, reactions, events={}):

    doc = create_model(species, parameters, reactions, events)

    # Convert sbml document into a string for basico
    sbml_str = sbml.writeSBMLToString(doc)
    return sbml_str




