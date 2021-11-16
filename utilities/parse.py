#!/usr/bin/env python3

#####################################################
#                                                    #
#    Low-level parser for reading in a tab-delimited #
#    parameter file and for reading out a            #
#    tab-delimited data file                         #    
#                                                    #
#####################################################

import sys, os, inspect, re
from pathlib import Path
import logging

logging.basicConfig(level=logging.WARNING)
module_logging = logging.getLogger(__name__)
logging_array = [module_logging]

this_file_path = Path(inspect.getfile(inspect.currentframe()))
this_dir = this_file_path.absolute().parents[0]
LIGHT_DIR = this_file_path.absolute().parents[1]
root_dir = this_file_path.absolute().parents[3]
sys.path.append(str(this_dir))
sys.path.append(str(root_dir))

from util import bool, rng, timestamp
from utilities.TSDReader import readtsd
import math
import pandas as pd


def params(param_file):

    """
    :param param_file: name of the parameter file to parse
    :param module_logger: logger object to keep track of simulation
    :return: simulation parameters in a dictionary, see read me for more detail
    """

    # add parameter folder
    # TODO: Check if file exists locally first
    if "/" not in param_file:
        param_file = os.path.join("params/", param_file)

    param_file = os.path.join(LIGHT_DIR, param_file)

    if not os.path.isfile(param_file):
        logging_array[0].error(f"Cannot find param file {param_file}. Check its path.")
        exit(1)

    logging_array[0].info(f"Starting simulation with param file: {param_file}")

    # EXPECT A TAB DELITED FILE: name\value\dtype
    if not os.path.isfile(param_file):
        logging_array[0].error(" Can't find parameter file at: " + str(param_file))
        assert False  # unknown file

    parameter_lists = text_separator(param_file)
    pdict = {}

    for param in parameter_lists:

        if len(param) not in [3]:
            logging_array[0].error("Problem parsing text file in parameter line  " + str(param))
            exit(1)

        if param[1][0] == '[': #ie is list
            param[1] = param[1].replace('[','').replace(']','')
            param[1] = param[1].split(',')

            is_list = True
        else:
            is_list = False

        if param[2] in ['str','string']:
            if not is_list:    val = param[1]
            else: val = [p for p in param[1]]
        elif param[2] == 'int':
            if not is_list: val = int(float(param[1]))
            else: val = [int(p) for p in param[1]]
        elif param[2] == 'float':
            if not is_list: val = float(param[1])
            else: val = [float(p) for p in param[1]]
        elif param[2] == 'bool':
            if not is_list: val = bool(param[1])
            else: val = [bool(p) for p in param[1]]
        elif param[2] == 'eval':
            if not is_list: val = eval(param[1])
            else: val = [eval(p) for p in param[1]]
        elif param[2] == 'exp':
            pieces = param[1].split('e')
            val = float(pieces[0]) * math.pow(10,float(pieces[1]))
        elif param[2] == 'csvcolumn':
            pieces = param[1].split('[')
            assert ( len(pieces) == 2 )
            assert ( pieces[1][-1] == ']' )
            filename = pieces[1][:-1]
            colname = pieces[0]
            # reead csv
            df = pd.read_csv(filename,
                sep=',',
                comment='#')
            val = df[colname][0]
        else:
            logging_array[0].error(f'Unknown datatype for parameter "{param[0]}".')
            exit(1)

        pdict[param[0]] = val

    # Get file name and timestamp for output
    pdict["file_name"] = param_file.split("/")[-1][:-4]
    pdict['timestamp'] = timestamp()

    # what about cps dir?
    pdict['cps_dir'] = os.path.join(LIGHT_DIR, 'cps_dir/')

    if pdict['output_dir'][-1] != "/":
        pdict['output_dir'] += "/"

    # Set up output parameter
    # Careful: the ouputfile name does not have an extension
    # TODO: CHANGE THIS CODE TO OUPUT OUTSIDE PROGRAM FOLDER
    if pdict['output_dir'].lower() == "default":
        pdict['output_dir'] = os.path.join(LIGHT_DIR, "output/" + pdict['file_name'])
    else:
        pdict['output_dir'] = os.path.join(LIGHT_DIR, pdict['output_dir'])

    if pdict['output_file'].lower() == "default":
        pdict['output_file'] = os.path.join(pdict['output_dir'], pdict["file_name"] + "_" + pdict['timestamp'] + ".pkl")
    else:
        pdict['output_file'] = os.path.join(pdict['output_dir'], pdict["output_file"] + ".pkl")

    if 'solver_log' not in pdict.keys():
        pdict['solver_log'] = None #TODO: should clean this, defaults to std.out is 'solver_log' is None

    if pdict['simulation_method'] in ['deterministic', 'Deterministic', 'LSODA']:
        pdict['repetitions'] = 1
        logging_array[0].warning("Since the simulation is deterministic only one repetition will be executed")

    return pdict


def tsd(params, root_light_dir=True):
    """
    Time is assumed to be FLOAT, all other also assumed to be FLOAT
    
    Returns:
    ...

    """
    data = None

    tsd_file = params['output_file']
    if root_light_dir:
        tsd_file = os.path.join(LIGHT_DIR, tsd_file)

    if params['solver'].lower() in ['copasi','either']:
        with open(tsd_file, 'r') as file:
            lines = file.readlines()
            for i in range(len(lines)):
                line = lines[i]
                line = line.replace('\n','')
                parts = line.split('\t')
                if i==0:
                    data = {part.replace('[','').replace(']',''):[] for part in parts}
                    header = [part.replace('[','').replace(']','') for part in parts]
                else:
                    for j in range(len(parts)):
                        part = float(parts[j])
                        data[header[j]] += [part]

    elif params['solver'].lower() == 'ibiosim':
        with open(tsd_file, 'r') as file:
            header, data = readtsd(tsd_file, epsilon=0.0)
            data['Time'] = data['time']

    return data


def read_csv_target_file(target_file, transpose=False, reformat=True, time_units='minutes'):
    """
    If reformat returns data as:
    { 'Time': [t0, t1, ..]
      species0 : {'runs': [run0, run1, run2, ..]},
      species1 : {'runs': [run0, run1, run2, ..]}, ... }
    """

    if time_units in ['hour', 'hours']:
        time_multiplier = 60
    elif time_units in ['minutes','min']:
        time_multiplier = 1
    else:
        print("\nparse: read_csv_target_file: unrecognized time_units: ",time_units,'\n')
        assert(False)

    if not transpose:
        with open(target_file) as f:
            lines = f.readlines()
            lines = [lines[i].replace('\n','') for i in rng(lines) if not lines[i][0] == '#']
            header = lines[0].split(',')
            lines = lines[1:]
            data = {k: [] for k in header}
            for i in range(len(lines)):
                line = lines[i].split(',')
                for j in range(len(header)):
                    if header[j] in ['Time','time']:
                        data[ 'Time' ] += [ float(line[j]) * time_multiplier ]
                    else:
                        data[ header[j] ] += [ float(line[j]) ]
    else:
        assert False, "fix me"
        """
        with open(target_file) as f:
            lines = f.readlines()
            lines = [lines[i].replace('\n','') for i in rng(lines) if not lines[i][0] == '#']
            data = {}
            for i in range(1,len(lines)):
                line = lines[i].split(',')
                line = [line[0]] + [float(line[i]) for i in range(1,len(line))]
                if line[0] in ['Time','time']:
                    line = [line[0]] + [line[i]*time_multiplier for i in range(1,len(line))]
                data[line[0]] = line[1:]
        """

    if reformat:
        re_data={ 'Time': data['Time'] }
        for key in data.keys():
            if key not in ['Time','time']:
                re_data[key] = { 'runs': [data[key]] }
        data = re_data

    return data


# This is the file separator function
# It is responsible for the interpretation of the text file
# Parses the elements into 3 strings
def text_separator(param_file):

    """
        File reader
        Input: file name or directory
        Output: lists of lists of triplet parameters
    """

    with open(param_file, "r") as file:
        parsed_lines = []

        # Discart first line
        file.readline()

        for line in file.readlines():

            # Ignore comments and empty lines
            ci = line.find("#")
            if ci != -1:
                line = line[0: ci]

            if not re.search("[a-zA-Z0-9]+", line):
                continue

            # Pre treatment removal of whitespaces and stuff
            line = line.replace('\t', ' ')
            line = line.strip()
            line = line.strip("\n")

            # First and third parameters
            first_param = re.match("^\S+", line).group(0)
            third_param = re.match("^\S+", line[::-1]).group(0)[::-1]

            # Second parameter (middle one):
            line = line.replace(first_param, "")
            line = line.replace(third_param, "")

            isec = re.search("[a-zA-Z0-9()]", line).start()
            li = re.search("[a-zA-Z0-9()]", line[::-1]).start()
            fsec = len(line) - li
            sec_param = line[isec: fsec]

            parsed_lines.append([first_param, sec_param, third_param])

    return parsed_lines

######################################################################

