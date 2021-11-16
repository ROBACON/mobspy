#!/usr/bin/env python3

#####################################################
#                                                   #
#   E.coli Donor cells (for Fig 7)                  #
#                                                   #
#                                                   #
#####################################################

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

module_logger = logging.getLogger('root')
VOLUME_ML = 1


# antibiotics
def __resistant(ecoli: str):
    """
    Returns: if resistant to ADonor
    """
    return False


def __can_uptake_antibiotic(ecoli: str, params):
    """
    Returns: if can uptake antibiotic
    """
    # return True
    return (__antibiotic_slots(ecoli) < 1)


def __antibiotic(ecoli: str):
    """
    Returns: if antibiotic inside
    """
    return (__antibiotic_slots(ecoli) > 0)


def __antibiotic_slots(ecoli: str):
    """
    Returns: number of antibiotc slots taken

    e.g. 'E_ADonor0' -> 0
         'E_ADonor1' -> 1
    """
    parts = ecoli.split(f'_ADonor')
    parts2 = parts[1].split('_')
    return int(parts2[0])


def __is_resistant_to_antibiotics_inside(ecoli: str):
    return (not __antibiotic(ecoli)) or __resistant(ecoli)


def __ecoli_birthdeathrates(ecoli: str, params: dict):
    """
    In: ecoli as str
    Returns: birth rates (per resource), and
             death rate
    """
    global VOLUME_ML
    rates = {
        'birth': {
            'G': None,
            'AA': None
        },
        'death': None
    }
    init_glucose = int( params['init_glucose'] * VOLUME_ML )
    init_amino_acid = int( params['init_amino_acid'] * VOLUME_ML )
    
    # birth
    if __antibiotic(ecoli):
        if __is_resistant_to_antibiotics_inside(ecoli):
            rates['birth']['G']  = f"dup_rate * antibiotic_resistant_birth_penalty_factor / {init_glucose}"
            rates['birth']['AA'] = f"dup_rate * amino_acid_penalty_factor * antibiotic_resistant_birth_penalty_factor / {init_amino_acid}"
        else:
            rates['birth']['G']  = f"dup_rate * antibiotic_birth_penalty_factor / {init_glucose}"
            rates['birth']['AA'] = f"dup_rate * amino_acid_penalty_factor * antibiotic_birth_penalty_factor / {init_amino_acid}"

    else:
        # no antibiotic inside
        rates['birth']['G']  = f"dup_rate / {init_glucose}"
        rates['birth']['AA'] = f"dup_rate * amino_acid_penalty_factor / {init_amino_acid}"

    # death
    if __antibiotic(ecoli):
        if __is_resistant_to_antibiotics_inside(ecoli):
            rates['death'] = f"antibiotic_resistant_death_rate"

        else:
            rates['death'] = f"antibiotic_death_rate"

    else:
        rates['death'] = 'death_rate'

    return rates


def build(params):

    # potentially set volume
    global VOLUME_ML
    if 'volume_ml' in params.keys():
        VOLUME_ML = params['volume_ml']
        print(f'Donor: volume has been set to {VOLUME_ML} ml.') 


    # --- events --- #
    events = {}

    # --- init_food --- #
    species = {
        'G':            params['init_glucose'],
        'AA':           params['init_amino_acid'],
        'P':            params['init_phages'],
        'dead':         0,
        'E_ADonor0':    params['init_ecoli'],
        'E_ADonor1':    0,
        'ADonor':       params['init_antibiotics_ADonor'],
    }


    # --- species & mappings --- #

    # Various E. coli: (description is possibly outdated)
    # E_{infection_state}_{antibiotic_resistant_state}
    # where infection_state = {not_infected, early_infection, late_infection}
    # and antibiotic_resistant_state = {not_resistant, early_resistant, late_resistant}

    # keep each __identifier__ unique, since the code uses string replacement for reactions
    mappings = {
        'Resources':['G','AA'],
        'Phages':['P'],
        'Ecoli':['E_ADonor0', 'E_ADonor1'],
        'OD': '0.01 * [dead] + [Ecoli]'
    }

    # potentially dilute to volume less than 1ml  
    for s in species.keys():
        species[s] = int( species[s] * VOLUME_ML )


    # --- reaction parameters --- #
    # note that some parameters are hardcoded into reaction rates or elsewhere
    # instead of being listed here
    parameters = {
        'dup_rate':                         (params['dup_rate'],                        'per_min'),
        'amino_acid_penalty_factor':        (params['amino_acid_penalty_factor'],       'dimensionless'),
        'death_rate':                       (params['death_rate'],                      'per_min'),
        'sec_rate_G':                         (params['sec_rate_G'],                        'per_min'),
        'sec_rate_AA':                         (params['sec_rate_AA'],                        'per_min'),
        'antibiotic_birth_penalty_factor':              (params['antibiotic_birth_penalty_factor'],             'dimensionless'),
        'antibiotic_resistant_birth_penalty_factor':    (params['antibiotic_resistant_birth_penalty_factor'],   'dimensionless'),
        'antibiotic_death_rate':                        (params['antibiotic_death_rate'],                       'per_min'),
    }

    reactions = {}

    # PHAGE DECAY
    reactions[f'decay_P'] = {
        're': [(1, 'P')],
        'pr': [ ],
        'kin': f"P * {params['P_decay_rate']}"
    }

    for ecoli in ['E_ADonor0', 'E_ADonor1']:
        birthdeathrates = __ecoli_birthdeathrates(ecoli, params)

        # BIRTH events
        if __is_resistant_to_antibiotics_inside(ecoli):
            # normal duplication
            l_child = r_child = ecoli
            for res in ['G','AA']:
                # birth per resource
                reactions[f'birth_{ecoli}_{res}'] = {
                    're': [(1, ecoli), (3, res)],
                    'pr': [(1, l_child), (1, r_child)],
                    'kin': f"{ecoli} * {res} * {birthdeathrates['birth'][res]}"
                }
        else:
            # death with certain probability
            for res in ['G','AA']:
                reactions[f'antibiotic_death_on_birth_ADonor_{ecoli}_{res}'] = {
                    're': [(1, ecoli), (1, res)],
                    'pr': [(1, 'dead'), (1, res)],
                    'kin': f"{ecoli} * {res} * {params['antibiotic_death_on_duplication_prob_ADonor']} * {birthdeathrates['birth'][res]}"
                }

        # DEATH
        reactions[f'death_{ecoli}'] = {
            're': [(1, ecoli)],
            'pr': [(1, 'dead')],
            'kin': f"{ecoli} * {birthdeathrates['death']}"
        }

        # SECRETION
        # reactions[f'secrete_{ecoli}'] = {
        #    're': [(1, ecoli), ],
        #    'pr': [(1, ecoli), (1,'P'), ],
        #    'kin': f"sec_rate * {ecoli}",
        # }
        for res in ['G','AA']:
                # birth per resource
                reactions[f'secrete_{ecoli}_{res}'] = {
                    're': [(1, ecoli), (1, res)],
                    'pr': [(1, ecoli), (1, res), (1,'P'), ],
                    'kin': f"sec_rate_{res} * {ecoli} * {res} * {birthdeathrates['birth'][res]} / dup_rate",
                }

        # ANTIBIOTICS
        if __can_uptake_antibiotic(ecoli=ecoli, params=params):
            reactions[f'antibiotic_exposure_ADonor_{ecoli}'] = {
                're': [(1, ecoli),(1, 'ADonor')],
                'pr': [(1, ecoli.replace(f'_ADonor0',f'_ADonor1'))],
                'kin': f"ADonor * {ecoli} * {params['antibiotic_uptake_rate_ADonor']} / {VOLUME_ML}"
            }

    return {'species':species, 'parameters':parameters, 'reactions':reactions, 'events':events, 'mappings':mappings}


if __name__ == '__main__':
    params = parse.params(sys.argv[1])
