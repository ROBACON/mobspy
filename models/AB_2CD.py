import sys, os, inspect
from pathlib import Path
import logging

local_dir = Path(inspect.getfile(inspect.currentframe()))
root_dir, sim_dir = local_dir.absolute().parents[0], local_dir.absolute().parents[1]
sys.path.append(str(root_dir))
FILENAME = os.path.basename(__file__)

from pprint import pprint
import run
import json, math
from copy import deepcopy
import itertools
import numpy as np

module_logger = logging.getLogger('root')

VOLUME_ML = 1


def build(params):

    # potentially set volume
    global VOLUME_ML
    if 'volume_ml' in params.keys():
        VOLUME_ML = params['volume_ml']
        # print(f'Volume has been set to {VOLUME_ML} ml.')

    events = {}

    # Here we define the species used in the reaction
    # And assign their initial quantity
    species = {
        'A': params['init_A'],
        'B': params['init_B'],
        'C': params['init_C'],
        'D': params['init_D']
    }

    # Mapping between names and species
    mappings = {
        'reactant A': '2*[A]',
        'reactant B': ['B'],
        'product C': ['C'],
        'product D': ['D']
    }

    # potentially dilute to volume less than 1ml
    for s in species.keys():
        species[s] = int(species[s] * VOLUME_ML)

    # setting the reaction and the parameters
    # --- reaction parameters --- #
    parameters = {
        'res_1_rate': (params['res_1_rate'], 'per_min'),
    }

    reactions = {
        'res_1_transformation': {
            're': [(1, 'A'), (1, 'B')],
            'pr': [(2, 'C'), (1, 'D')],
            'kin': f"A * B * res_1_rate"
        },
    }

    return {'species': species, 'parameters': parameters, 'reactions': reactions, 'events': events,
            'mappings': mappings}