#!/usr/bin/env python3

#################################################
#                                               #
#   Primary script for running a simulation     #
#                                               #
#################################################

import sys, os, inspect
from pathlib import Path
import argparse
import logging

logging.basicConfig(level=logging.WARNING)
module_logging = logging.getLogger(__name__)
logging_array = [module_logging]

import matplotlib.pyplot as plt


local_dir = Path(inspect.getfile(inspect.currentframe()))
root_dir, sim_dir = local_dir.absolute().parents[0], local_dir.absolute().parents[1]
sys.path.append(str(root_dir))
sys.path.append(str(sim_dir))


from models import builder
from utilities import parse, pybash, util, log
from tools import rewrite_cps, optimize, init_solver
import math
import copy
import pickle
from time import time
from pprint import pprint
import contextlib
import joblib
from tqdm import tqdm
from joblib import Parallel, delayed
import tempfile
import shutil
import basico

import numpy as np


# path progress bar
@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


############################ PRIMARY SIMULATION FUNCTIONS ####################################


def simulate(params, verbose=False, plot_run=True, rescale_to_1ml=True):

    """
        This function coordinates the simulation by calling the necessary jobs
        In the future we hope to implement parallel cluster computing compatibility
    :param params: simulation parameters from the text file
    :param verbose: #TODO re setup verbose
    :param plot_run: #TODO ask about this
    :param rescale_to_1ml: rescale volume
    :return: COPASI simulated data
    """


    # Run in parallel or sequentially
    # If nothing is specified just run it in parallel
    try:
        if params["jobs"] == 1:
            logging_array[0].info("Running simulation sequentially")
            jobs = params["jobs"]
        else:
            logging_array[0].info("Running simulation in parallel")
            jobs = params["jobs"]
    except KeyError:
        logging_array[0].info("Running simulation in parallel")
        jobs = -1

    # TODO: If cluster compatibility is added I suggest here
    sbml_str, species, mappings = builder.build(params)
    data = job_execution(sbml_str, params, jobs)
    data = remap_species(data, params, mappings)

    # potentially print state
    if 'print_final_state' not in params.keys():
        params['print_final_state'] = False
    if 'print_initial_state' not in params.keys():
        params['print_initial_state'] = False
    if 'print_header' not in params.keys():
        params['print_header'] = False

    # potentially rescale data to 1ml
    # TODO: Why is this here?
    if rescale_to_1ml:
        if 'volume_ml' in params.keys():
            VOLUME_ML = params['volume_ml']
            for s in data.keys():
                if s != 'Time':
                    # real species
                    # print(f'rescale from: {data[s]}')
                    runs = data[s]['runs']
                    new_runs = [[value / VOLUME_ML for value in run] for run in runs]
                    data[s]['runs'] = new_runs
                    # print(f'rescale to  : {data[s]}')

    # potentially save data
    if params['save_data']:
        logging_array[0].debug("Saving data (reason: parameter <save_data>)")

        # TODO: Why none assignment before this?
        # plot_params = None
        # if plot_run and params['plot']:
            # load and also include plot params
            # plot_params = plot.get_plot_params(plot_json_filename=params['plot_json'])
        pickled = {'data': data,
                   'mappings': mappings,
                   'params': params}

        if not os.path.isdir(params['output_dir']):
            logging_array[0].info("Creating desired output directory: %s..." % (params['output_dir']))
            os.makedirs(params['output_dir'], exist_ok=True)

        # Save pickle data
        util.pickle_it(pickled, params['output_file'])

    else:
        logging_array[0].warning("NOT saving data (reason: parameter <save_data>)")

    if plot_run and params['plot']:
        # why 2 toggles for plot? One is if calling many times, such as during optimize, the other is for user control
        plot.plot_data_fromjsonfile(plot_json_filename=params['plot_json'], data=data)

    logging_array[0].info("Simulation is Over")

    return {'data': data, 'params': params, 'mappings': mappings}


def job_execution(sbml_str, params, jobs):

    # This is defined for parallelism purposes
    # THERE MUST BE ONE TEMP DIRECTORY FOR EVERY COPASI FILE
    # OTHERWISE PARALLELISM DOES NOT WORK - COPASI OVERWRITES THE OUTPUT
    def __single_run(packed):
        sbml_str, i = packed
        basico.model_io.load_model_from_string(sbml_str)

        data = basico.run_time_course(params['duration'])
        reformated_data = reformat_time_series(data)

        return reformated_data

    parallel_data = Parallel(n_jobs=jobs)(delayed(__single_run)((sbml_str, i)) for i in range(params['repetitions']))

    if not parallel_data:
        logging_array[0].error("Error: The parallel model has not produced an output.")
        logging_array[0].error("Try addding (sequential True bool) to parameters")
        exit(1)


    # We always call merge to keep the data in the format we want
    merged_data = merge(params, parallel_data)

    return merged_data


################################## HELPER FUNCTIONS ####################################


def reformat_time_series(data):

    data_dict = {'Time': data.index.tolist()}

    for key in data:
        data_dict[key] = list(data[key])

    return data_dict


def merge(params, data):
    """
    a basic merge of the data with possible resampling in time to save less data

    params: dict of parameters
    data:   Array of single runs: data[0], data[1], ...

    returns: Data merged in a dictionary with all time series from all experiments
                separated by species and mapping
    """

    for i in range(1,len(data)):
        assert (data[i-1]['Time'] == data[i]['Time'])  # should be set at exact same times

    merged_data = {'Time': data[0]['Time']}

    # TODO: Ask them about this merger, why is it here?
    # Check before
    # generate 'Time'
    """
    if params['stats_stepsize'] is None or params['stats_stepsize'] == 0:
        # take time points from first run
        merged_data = {'Time': data[0]['Time']}
        time_points = data[0]['Time']
    else:
        # construct time points by resampling time
        merged_time, time_points = [], []
        for t in range(len(data[0]['Time'])):
            if t * params['sim_stepsize'] % params['stats_stepsize'] == 0:
                merged_time += [data[0]['Time'][t]]
                time_points += [t]
        merged_data = {'Time': merged_time}
    """

    for key in data[0].keys():
        if key not in ['Time']:
            merged_data[key] = {'runs': [data[i][key]
                                         for i in range(params["repetitions"])],
                                }

    return merged_data


def remap_species(data, params, mapping):
    """
    Takes the simulated species (data) and add ones defined by
    mapping (mapping).

    By default, a mapping is a list of species in which case the
    sum is taken: ['a', 'b', ...]

    A different behavior can be specified by a string, however:
    '[a] * 2 + 7 - [b]'.
    MIND: dont forget the '[]'.
    """

    mapped_data = {'Time': data['Time']}
    T = range(len(data['Time']))

    # copy over all unmapped ones
    for k in data.keys():
        mapped_data[k] = data[k]

    # 1st pass with sum mappings
    for group in mapping.keys():
        the_mapping = mapping[group]
        mapped_data[group] = {'runs': []}

        try:
            # check if is a list -> sum
            if type(the_mapping) is list:

                for run in range(len(data[k]['runs'])):
                    species = the_mapping
                    mapped_data[group]['runs'].append([sum([data[species[i]]['runs'][run][t]
                                                            for i in range(len(species))]) for t in T])

        except Exception as e:
            logging_array[0].error(f'run: remap_species: error when remapping "{the_mapping}".')
            logging_array[0].error('Possible fix: All runs must have the same time')
            logging_array[0].error(e)
            exit(1)
    # save as data
    data = mapped_data

    # 2nd pass with function mappings
    for group in mapping.keys():
        the_mapping = mapping[group]

        try:

            # check if is str -> treat as formula
            if type(the_mapping) is str:

                for run in range(len(data[k]['runs'])):

                    mapped_data[group]['runs'].append([])

                    # rewrite: [S] -> data['S']['runs'][run][t]
                    the_mapping_new = the_mapping.replace("[", "data['").replace("]", "']['runs'][run][t]")
                    # map with local scope
                    for t in T:
                        locs = locals()
                        mapped_data[group]['runs'][run].append(eval(the_mapping_new, locs))

        except Exception as e:
            logging_array[0].error(f'run: remap_species: error when remapping "{the_mapping}".')
            logging_array[0].error(e)
            exit(1)

    return mapped_data

########################################################################################
####################################### MAIN FUNCTION AND API CALL #####################


def run_api(text_file, log_file=None):

    """
        This runs the simulation through an api call
        :param text_file: text file with the simulations descriptions
        :return: data resulting from simulation, also saves the simulation on output file if specified
    """

    # Call our lovely parser function to read the txt files
    params = parse.params(text_file)

    # Return copasi simulated data
    return simulate(params)


if __name__ == "__main__":
    """
        main is used to run the program through the terminal
        just type python3 run.py PARAMETER_FILE_NAME
        other arguments are available check read me to understand them
    """

    parser = argparse.ArgumentParser(description='Simulation framework for CRNs.')
    parser.add_argument('PARAM_FILE',
                        nargs=1,
                        type=str,
                        help='the parameter file')
    parser.add_argument("-log", "--log",
                        default="info",
                        help=(
                            "Provide logging level. "
                            "Example --log debug', default=info")
                        )

    # Get arguments
    args = parser.parse_args()


    # Create logging file
    log_file = "logs/last_logs.txt"
    logging_array[0].info("Starting log process")

    # Parse de parameter file and add log_file name
    params = parse.params(args.PARAM_FILE[0])

    # Calling the simulation coordinator function
    simulate(params)
