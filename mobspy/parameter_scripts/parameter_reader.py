"""
    This module is responsible for processing the parameters given to MobsPy before starting the simulation
"""
import json
import sys
from datetime import datetime
import os
from pathlib import Path
import numpy as np
import mobspy.simulation_logging.log_scripts as simlog
from pint import Quantity


def read_json(json_file_name):
    """
        Reads json file

        :param plot_json_filename: json file name

        :raise simlog.error: If the file was not able to be read

        :return: plot parameter dictionary
    """
    with open(json_file_name, 'r') as file:
        try:
            json_data = json.load(file)
        except json.decoder.JSONDecodeError as e:
            simlog.error('Error reading file')
            exit(1)

    return json_data


def __name_output_file(params):
    """
        Gives a name to the output file - just date time in case the user has not specified one

        :param params: (dict) Dictionary with simulation parameters
    """

    file_name = "r_"
    file_name += str(datetime.now()) + '.json'
    params["absolute_output_file"] = params["output_dir"] + file_name


def __check_stochastic_repetitions_seeds(params):
    """
        The list of seeds must be equal to the number of repetitions specified

        :param params (dict) = Dictionary with simulation parameters
    """
    if 'seeds' in params:
        try:
            if params['repetitions'] != len(params['seeds']):
                simlog.error('Seeds must be equal to the number of repetitions')
        except Exception:
            simlog.error('Parameter seeds must be a list')


def __convert_parameters_for_COPASI(params):
    """
        Converts parameters units to MobsPy standard units (basiCO needs seconds for simulation duration)

        :param params: (dict) = Dictionary with simulation parameters
    """
    for key, p in params.items():
        if isinstance(p, Quantity):
            if str(p.dimensionality) == '[time]':
                params[key] = p.to('second').magnitude
                continue


def parameter_process(params):
    __name_output_file(params)
    __check_stochastic_repetitions_seeds(params)
    __convert_parameters_for_COPASI(params)


