#!/usr/bin/env python3


#################################################
#                                               #
#   Merge data and compute statistical tests    #
#   for batch stochastic runs                   #
#                                               #
#################################################

import sys, os, inspect
from pathlib import Path

local_dir = Path(inspect.getfile(inspect.currentframe()))
root_dir, sim_dir = local_dir.absolute().parents[0], local_dir.absolute().parents[2]
sys.path.append(str(root_dir))
sys.path.append(str(sim_dir))

import math, numpy as np
import confidence_intervals
import logging

from statistics import median

module_logger = logging.getLogger('root')



#############################################################

def get_statistics(data,
                   statistics='avg',
                   statistics_params={'hull_CI': True}):
    """
    In:
        data:
        a dict of the form
        data['Time'] = [t0, t1, ...]
        data[species]['runs'] = [run0, run1, ...]

        with an optional own timebase as
        data[species]['Time'] = [run0, run1, ...]

        (this is used for external data which often comes with its own timebase)

        at the moment it also can have entries like (which are ignored):
        data['xy: Time'] = [t0, t1, ...]
      
        statistics:
        'avg', 'med', 'runs', ...

        statistics_params:
        parameters for the statistics

    Out:
        lines & shades of the form 
        lines_and_shades['Time'] = [t0, t1, ...]
        lines_and_shades[species] = {'x':      [t0, t1, ...]
                                     'lines':  [line0, line1, ...]
                                     'shades': [[shade0_upper_line, shade0_lower_line], ...]}

    """
    # generate function
    if statistics == 'avg':
        generate = lambda time, runs : stat_avg(time, runs, statistics_params)
    elif statistics == 'med':
        generate = lambda time, runs : stat_median(time, runs, statistics_params)
    elif statistics == 'runs':
        generate = lambda time, runs : stat_runs(time, runs, statistics_params)
    else:
        module_logger.error(f'Unkown statistics: {statistics}.')
        module_logger.error('Currently known: avg, med, runs')
        exit(1)

    # generate for species
    lines_and_shades = {}
    for species in data.keys():
        if 'Time' not in species:
            # NOTE: 'not in' is important here, since sometimes 'xy: Time' may occur in merged data. We skip this at the moment.
            
            # check if species has its own timebase
            if 'Time' in data[species].keys():
                # yes -> use it
                time = data[species]['Time']
            else:
                # no -> use standard timebase
                time = data['Time']

            # generate lines & shades for this species
            lines_and_shades[species] = generate(time, data[species]['runs'])

    return lines_and_shades



############# some statistics ################################################

def stat_avg(time, runs,
    statistics_params={'hull_CI': True}):
    """
    runs -> lines & shades for this statistics
    """
    if not all([ len(run) == len(time) for run in runs] ):
        module_logger.error('avg: not all runs in <runs> have the same length as <time>')
        print(f'Time has length {len(time)}.')
        print(f'The runs have lengths: {[ len(run) for run in runs]}')
        exit(1)

    avg = []
    CI_min = []
    CI_max = []
    for t in range(len(time)):
        a = [ run[t] for run in runs ]
        N = len(a)

        # average
        avg += [ sum(a)/N ]
        CI_min_new, CI_max_new = calc_CI(a, statistics_params)
        CI_min += [ CI_min_new ]
        CI_max += [ CI_max_new ]

    # only return shades if there are multiple runs
    if len(runs) > 1:
        return {'x':      time,
                'lines':  [ avg ],
                'shades': [ [CI_min, CI_max] ] }
    else:
        return {'x':      time,
                'lines':  [ avg ],
                'shades': [] }


def stat_median(time, runs,
    statistics_params={'hull_CI': True}):
    """
    runs -> lines & shades for this statistics
    """
    if not all([ len(run) == len(time) for run in runs] ):
        module_logger.error('avg: not all runs in <runs> have the same length as <time>')
        print(f'Time has length {len(time)}.')
        print(f'The runs have lengths: {[ len(run) for run in runs]}')
        exit(1)

    med = []
    CI_min = []
    CI_max = []
    for t in range(len(time)):
        a = [ run[t] for run in runs ]
        N = len(a)

        # median
        med += [ median(a)/N ]
        CI_min_new, CI_max_new = calc_CI(a, statistics_params)
        CI_min += [ CI_min_new ]
        CI_max += [ CI_max_new ]

    # only return shades if there are multiple runs
    if len(runs) > 1:
        return {'x':      time,
                'lines':  [ avg ],
                'shades': [ [CI_min, CI_max] ] }
    else:
        return {'x':      time,
                'lines':  [ avg ],
                'shades': [] }



def stat_runs(time, runs,
    statistics_params={}):
    """
    runs -> lines & shades for this statistics
    """
    if not all([ len(run) == len(time) for run in runs] ):
        module_logger.error('avg: not all runs in <runs> have the same length as <time>')
        print(f'Time has length {len(time)}.')
        print(f'The runs have lengths: {[ len(run) for run in runs]}')
        exit(1)

    return {'x':      time,
            'lines':  runs,
            'shades': [] }



############# some helpers ################################################

def calc_CI(a, statistics_params):
    if statistics_params['hull_CI']:
        CI_min, CI_max = min(a), max(a)
    elif not np.any(np.array(a)): #if a is only zeroes
        CI_min, CI_max = 0,0
    elif statistics_params['skewed_CI']:
        CI_min, CI_max = confidence_intervals.t_CI(a, confidence=statistics_params['confidence'], skewed=True)
    elif statistics_params['repetitions'] < 30:
        CI_min, CI_max = confidence_intervals.t_CI(a, confidence=statistics_params['confidence'], skewed=False)
    else:
        CI_min, CI_max = confidence_intervals.normal_CI(a, confidence=statistics_params['confidence'])
    
    return CI_min, CI_max