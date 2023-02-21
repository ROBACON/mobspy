"""
    mobspy.data_handler.process_result_data.py

    Handles converting the output data from a simulation into desired-units or concentration
"""
from pint import UnitRegistry, Quantity
from scipy.constants import N_A
from copy import deepcopy

u = UnitRegistry()


def convert_data_to_desired_unit(data, unit_x=None, unit_y=None, output_concentration=False, volume=1):
    """
        Converts the data from the MobsPy standard units to the desired units specified by the user

        Parameters:
            data (dict) - resulting data from a MobsPy simulation execution
                Accepts MobsPy results['data'] format

            unit_x (str) - (x axis unit - Time) unit that the user desires the time in
                Accepts strings with the unit name

            unit_y (str) - (y axis unit - Concentration or Counts) unit that the user desires the counts in
                Accepts strings with the unit name

            output_concentration (bool) - boolean - decide if output should be a concentration or count

            volume (float) - volume set as parameter to the simulation

        Returns:

            converted_data (dict) - input data converted to the desired units
                Output data in the MobsPy results['data'] format
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



