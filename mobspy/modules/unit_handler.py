"""
    This module is responsible for dealing with units in the modules directory level
    MobsPy standard units are Decimeter (Liter) - Second - Count
"""
from pint import UnitRegistry, Quantity
from scipy.constants import N_A
from copy import deepcopy
import os, sys
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.modules.unit_handler as uh


def convert_rate(quantity, reaction_order, dimension):
    """
        This function converts the rate from the users given unit to MobsPy standard units

        :param quantity: (int, float, Quantity) If it is a quantity object convert, otherwise it remains the same
        :param reaction_order: (int) number of reactants in the reaction, to check if the rate is in the correct unit
        :param dimension: (int) model's dimension (1D, 2D, 3D, ... )

        :param quantity: (int, float) converted unit into MobsPy standard units
    """
    volume_power = reaction_order - 1
    converted_quantity = deepcopy(quantity)

    # Dimension arrives as none but there is a rate
    if dimension is None and isinstance(quantity, Quantity) and reaction_order > 1:
        dimension = uh.extract_length_dimension(str(quantity.dimensionality), dimension,
                                                reaction_order)

    if isinstance(quantity, Quantity):
        try:
            if volume_power <= 0:
                if '[substance]' in str(quantity.dimensionality):
                    converted_quantity.ito(f'moles/seconds')
                    return converted_quantity.magnitude, dimension
                else:
                    converted_quantity.ito(f'1/seconds')
                    return converted_quantity.magnitude, dimension
            else:
                if '[substance]' in str(quantity.dimensionality):
                    converted_quantity.ito(f'decimeters ** {dimension*volume_power}/(moles ** {volume_power} * seconds)')
                    return converted_quantity.magnitude / (N_A ** volume_power), dimension
                else:
                    converted_quantity.ito(f'decimeters ** {dimension*volume_power}/seconds')
                    return converted_quantity.magnitude, dimension
        except Exception as e:
            simlog.error(str(e) + '\n' +
                         f'Problem converting rate {quantity} \n'
                         f'Is the rate in the form [volume]**{volume_power}/[time]?')

    else:
        return quantity, dimension


def convert_counts(quantity, volume, dimension):
    """
        This function converts the counts from the users given unit to MobsPy standard units. It also converts
        concentrations into counts

        :param quantity: (int, float, Quantity) If it is a quantity object convert, otherwise it remains the same
        :param volume: (int, float) volume in liters (converted beforehand)
        :param dimension: (int) model's dimension (1D, 2D, 3D, ... )

        :return: converted_quantity (int, float) = converted unit into MobsPy standard units
    """
    converted_quantity = deepcopy(quantity)

    if isinstance(quantity, Quantity):
        if '[length]' not in str(quantity.dimensionality) and '[substance]' not in quantity.dimensionality:
            simlog.error(f'The assigned quantity {quantity} is neither a count or concentration')

        try:
            if '[substance]' in quantity.dimensionality:
                if '[length]' in str(quantity.dimensionality):
                    converted_quantity.ito(f'moles/(decimeter ** {dimension})')
                    converted_quantity = converted_quantity*volume

                converted_quantity = converted_quantity.magnitude * N_A
            else:
                if '[length]' in str(quantity.dimensionality):
                    converted_quantity.ito(f'1/(decimeter ** {dimension})')
                    converted_quantity = converted_quantity*volume
                converted_quantity = converted_quantity.magnitude
        except Exception as e:
            simlog.error(str(e) + '\n' +
                         f'Problem converting rate {quantity} \n'
                         f'Is it really a count or concentration?')
    return converted_quantity


def check_dimension(dimension, value):
    """
        Checks for dimension consistency. It "stores" the first dimension it was given by returning it

        :param dimension: (int) model's dimension (1D, 2D, 3D ...)
        :param value: (int) dimension value being analysed

        :raise simlog.error: If dimensions are not consistent through the given units (units in 1D with 2D mixed)

        :return: dimension (int) = model's dimension (1D, 2D, 3D ...)
    """
    if dimension is None:
        dimension = int(value)
    else:
        if dimension != int(value):
            simlog.error('The dimensions are not consistent. There are at least two units given for different '
                         'dimension models')
    return dimension


def extract_length_dimension(unit_string, dimension, reaction_order=None):
    """
        Extracts the volume dimension from a Quantity object from Pint

        :param unit_string: (str) unit in str format
        :param dimension: (int) model's dimension (1D, 2D, 3D ...)
    """
    temp_list = unit_string.split()
    try:
        position = temp_list.index('[length]')
    except ValueError:
        position = -1

    if position == -1:
        return False
    if temp_list[position + 1] == '**':
        if reaction_order is None:
            dimension = check_dimension(dimension, temp_list[position + 2])
        else:
            temp_int = int(int(temp_list[position + 2])/(reaction_order - 1))
            dimension = check_dimension(dimension, temp_int)
    else:
        dimension = check_dimension(dimension, 1)
    return dimension


def convert_volume(volume, dimension):
    """
        Converts volume to decimetre**dimension

        :param volume: (int, float, Quantity) volume used in simulation

        :return: the converted volume in MobsPy units
    """
    if isinstance(volume, Quantity):
        dimension = extract_length_dimension(str(volume.dimensionality), dimension)
        volume = volume.to(f'decimeter ** {dimension}').magnitude
        return volume
    else:
        return volume


def convert_time(time):
    """
        Converts time into seconds

        :param time: (int, float, Quantity) any time used
    """
    if isinstance(time, Quantity):
        if str(time.dimensionality) == '[time]':
            return time.to('second').magnitude
    else:
        return time


if __name__ == '__main__':
    pass
