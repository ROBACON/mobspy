#!/usr/bin/env python3

#####################################################
#                                                   #
#   Picks and translates a python model into an     #
#   SBML readable by Copasi. New models must        #
#   be references here as well.                     #
#                                                   #
#####################################################


import sys, inspect, os
from pathlib import Path
import logging

module_logger = logging.getLogger('root')

this_file_path = Path(inspect.getfile(inspect.currentframe()))
this_dir = this_file_path.absolute().parents[0]
py_dir = os.path.join(this_file_path.absolute().parents[3], 'python')
sys.path.append(str(this_dir))
sys.path.append(str(py_dir))

from utilities.SBMLWriter import create_model
import libsbml as sbml, codecs

import importlib
import pkgutil


# automatic import of models
MODULES_DIR = str(this_dir)
discovered_plugins = { name: importlib.import_module(name)
    for finder, name, ispkg
    in pkgutil.iter_modules([MODULES_DIR]) }


def build(params):
    # the model must return a dictionary with entries for 'species', 'parameters', 'reactions', and 'events'
    # see AB.py for an example  

    try:
        model_name = params['model']
        model = discovered_plugins[model_name]

    except KeyError as e:
        module_logger.error(f"builder : the autamatic discovery of the model {params['model']} did not work : {e}")
        return

    # build it
    model = model.build(params)

    # finally print the model
    doc = create_model(species= model['species'], parameters= model['parameters'], reactions= model['reactions'], events= model['events'],debug=True) #params['debug'])

    # Convert sbml document into a string for basico
    sbml_str = sbml.writeSBMLToString(doc)

    # file name matters and must match run.py
    # with codecs.open(params['cps_dir']+'model.sbml', "w", "utf-8") as file:
        #file.write(sbml.writeSBMLToString(doc))

    if 'mappings' in model.keys():
        return sbml_str, model['species'], model['mappings']
    else:
        return sbml_str, model['species'], None




