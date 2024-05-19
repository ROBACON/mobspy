"""
    mobspy.data_handler.process_result_data.py

    Handles converting the output data from a simulation into desired-units or concentration
"""
from mobspy.modules.mobspy_expressions import u
from scipy.constants import N_A
from copy import deepcopy
import mobspy.simulation_logging.log_scripts as simlog


def extract_time_and_volume_list(list_of_params):
    """
        This function extracts the list of durations from all concatenated simulations (or one for single simulation)
        It also extracts the respective volumes at each simulation. It does not work if their is a change in volume
        and a simulation without fixed duration

        Conditional simulations have a parameter called '_end_condition' in their dictionary

        :param list_of_params: list of parameters of all concatenated simulations
    """
    no_fixed_volume = False
    no_fixed_dur = False
    initial_volume = list_of_params[0]['volume']
    for par in list_of_params:
        if par['_end_condition'] is not None:
            no_fixed_dur = True

        if par['volume'] != initial_volume:
            no_fixed_volume = True

    flag_concentration = True
    if no_fixed_dur and len(list_of_params) > 1 and no_fixed_volume:
        flag_concentration = False
        simlog.warning('Could not resolve simulation volume due to multiple simulations with at least one with a'
                       ' conditional duration. The output will be printed in counts instead')

    volume_list = []
    previous_time = list_of_params[0]['duration']
    sim_time_list = [previous_time]
    if no_fixed_volume:
        for i, par in enumerate(list_of_params):
            volume_list.append(par['volume'])

            if i == 1:
                continue
            else:
                current_time = previous_time + list_of_params[i]['duration']
                sim_time_list.append(current_time)
                previous_time = current_time
    else:
        volume_list = [initial_volume]
        sim_time_list = [sum([par['duration'] for par in list_of_params])]

    return volume_list, sim_time_list, flag_concentration


def convert_data_to_desired_unit(data, time_list, volume_list,
                                 unit_x=None, unit_y=None, output_concentration=False):
    """Converts the simulation output data from the MobsPy standard units
    to the desired units specified by the user

    :param data: (dict) resulting data from a MobsPy simulation execution
    :param time_list: (list) list of times where the volume changes (in case of single simulation it's only one)
    :param volume_list: (list) list of volumes changes (in case of single simulation it's only one)
    :param unit_x: (str) unit that the user desires the time in, defaults to dimensionless (None)
    :param unit_y: (str) unit that the user desires the y axis to be in (Concentration or counts)
    :param output_concentration: (bool) decide if output should be a concentration or count
    :return: converted_data - input data converted to the desired units
    :rtype: (dict) Dictionary meta-species as key and run as value
    """
    ur = u.unit_registry_object
    converted_data = deepcopy(data)

    if unit_x is not None:
        new_time = []
        for time in data['Time']:
            quantity = time * ur.seconds
            new_time.append(quantity.to(unit_x).magnitude)
        converted_data['Time'] = new_time

    def multiply_data_by_factor(data, factor):
        for key in data:
            if key == 'Time':
                continue

            converted_data[key] = [count * factor for count in data[key]]

    if output_concentration:
        converted_data = convert_to_concentration(data, converted_data, volume_list, time_list)

    # Fix here - forgot concentration conversion
    if unit_y is not None:
        if 'mol' in str(unit_y):
            multiply_data_by_factor(converted_data, N_A ** -1)
            if output_concentration:
                factor = (1*ur.molar).to(unit_y).magnitude
                multiply_data_by_factor(converted_data, factor)
            else:
                factor = (1*ur.moles).to(unit_y).magnitude
                multiply_data_by_factor(converted_data, factor)
        else:
            if output_concentration:
                factor = (1/ur.l).to(unit_y).magnitude
                multiply_data_by_factor(converted_data, factor)

    return converted_data


def convert_to_concentration(data, converted_data, volume_list, time_list):
    """
        Converts output data from counts to concentration according to simulation volume

        :param data: simulation data
        :param converted_data: data converted to rquested units
        :param volume_list: list of volumes of all simulations (more than one if concatenated)
        :param time_list: list of durations of each simulation (to check for respective volume in results)
    """
    new_data = {}
    for key in data:
        new_data[key] = []

    current_volume = volume_list[0]
    k = 0
    for i in range(len(data['Time'])):

        if data["Time"][i] > time_list[k] and len(volume_list) > 1:
            k = k + 1
            current_volume = volume_list[k]

        for key in data:
            if key == 'Time':
                continue
            else:
                new_data[key].append(data[key][i]/current_volume)

    new_data['Time'] = converted_data['Time']

    return new_data


