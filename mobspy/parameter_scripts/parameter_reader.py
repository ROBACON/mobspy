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
    In:
      plot_json_filename: json file name

    Returns:
      plot parameter dictionary
    """
    with open(json_file_name, 'r') as file:
        try:
            json_data = json.load(file)
        except json.decoder.JSONDecodeError as e:
            simlog.error('Error reading file')
            exit(1)

    return json_data


def __name_output_file(params):

    if params['output_dir'][0] == '/':
        params['output_dir'] = params['output_dir'][1:]

    try:
        main_directory = os.path.abspath(sys.modules['__main__'].__file__)
        save_dir = os.path.join(Path(main_directory).parent.absolute(), params['output_dir'])

        if params['output_file'] is None:
            file_name = "r_"
            file_name += str(datetime.now()) + '.json'
        else:
            file_name = params['output_file']

        params['output_absolute_directory'] = save_dir
        params['output_absolute_file'] = os.path.join(params['output_absolute_directory'], file_name)

    except:
        simlog.warning('Automatic data-saving setup failed. Please save manually')
        params['save_data'] = False


def __check_stochastic_repetitions_seeds(params):

    if 'seeds' in params:
        try:
            if params['repetitions'] != len(params['seeds']):
                simlog.error('Seeds must be equal to the number of repetitions')
        except Exception:
            simlog.error('Parameter seeds must be a list')


def __check_ode_repetitions(params):

    if params["simulation_method"].lower() == 'deterministic':
        params["repetitions"] = 1


def __convert_parameters_for_COPASI(params):
    for key, p in params.items():
        if isinstance(p, Quantity):
            if str(p.dimensionality) == '[time]':
                params[key] = p.to('second').magnitude
                continue


def parameter_process(params, mappings):
    __name_output_file(params)
    __check_stochastic_repetitions_seeds(params)
    __check_ode_repetitions(params)
    __convert_parameters_for_COPASI(params)


