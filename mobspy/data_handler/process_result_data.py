"""
    mobspy.data_handler.process_result_data.py

    Handles converting the output data from a simulation into desired-units or concentration
"""
from pint import UnitRegistry, Quantity
from scipy.constants import N_A
from copy import deepcopy

u = UnitRegistry()


def convert_data_to_desired_unit(data, unit_x=None, unit_y=None, output_concentration=False, volume=1):
    """Converts the simulation output data from the MobsPy standard units
    to the desired units specified by the user

    :param data: (dict) resulting data from a MobsPy simulation execution
    :param unit_x: (str) unit that the user desires the time in, defaults to dimensionless (None)
    :param unit_y: (str) unit that the user desires the y axis to be in (Concentration or counts)
    :param output_concentration: (bool) decide if output should be a concentration or count
    :param volume: (float) - volume set as parameter to the simulation
    :return: converted_data - input data converted to the desired units
    :rtype: (dict) Dictionary meta-species as key and run as value
    """
    u = UnitRegistry()
    converted_data = deepcopy(data)

    if unit_x is not None:
        new_time = []
        for time in data['Time']:
            quantity = time*u.seconds
            new_time.append(quantity.to(unit_x).magnitude)
        converted_data['Time'] = new_time

    def multiply_data_by_factor(data, factor):
        for key in data:
            if key == 'Time':
                continue

            converted_data[key] = [count*factor for count in data[key]]

    def convert_data(data, from_unit, to_unit):
        for key in data:
            if key == 'Time':
                continue

            converted_data[key] = [(count*from_unit).to(to_unit).magnitude for count in data[key]]

    if output_concentration:
        multiply_data_by_factor(converted_data, volume ** -1)

    if unit_y is not None:
        if 'mol' in unit_y:
            multiply_data_by_factor(converted_data, N_A ** -1)
            if output_concentration:
                convert_data(converted_data, u.molar, unit_y)
            else:
                convert_data(converted_data, u.moles, unit_y)

    return converted_data



