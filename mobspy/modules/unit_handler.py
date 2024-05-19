"""
    This module is responsible for dealing with units in the modules directory level
    MobsPy standard units are Decimeter (Liter) - Second - Count
"""
from pint import UnitRegistry, Quantity
from scipy.constants import N_A
from copy import deepcopy
import mobspy.simulation_logging.log_scripts as simlog


def convert_rate(quantity, reaction_order, dimension):
    """
        This function converts the rate from the users given unit to MobsPy standard units

        :param quantity: (int, float, Quantity) If it is a quantity object convert, otherwise it remains the same
        :param reaction_order: (int) number of reactants in the reaction, to check if the rate is in the correct unit
        :param dimension: (int) model's dimension (1D, 2D, 3D, ... )

        :param quantity: (int, float) converted unit into MobsPy standard units
    """
    volume_power = reaction_order - 1
    # For objects that cannot be deep copied
    try:
        converted_quantity = deepcopy(quantity)
    except:
        converted_quantity = quantity

    # Check to see if rate dimension is valid
    if dimension is None and isinstance(quantity, Quantity) and reaction_order > 1:
        dimension = extract_length_dimension(str(quantity.dimensionality), dimension,
                                             reaction_order)

    if isinstance(quantity, Quantity):
        try:
            if str(quantity.dimensionality) == '1 / [time]':
                converted_quantity = converted_quantity.convert(f'1/seconds')
                return converted_quantity.magnitude, dimension, True
            elif str(quantity.dimensionality) == '[substance] / [time]':
                converted_quantity = converted_quantity.convert(f'moles/seconds')
                return converted_quantity.magnitude * N_A, dimension, True
            elif '[substance]' in str(quantity.dimensionality):
                converted_quantity = converted_quantity.convert(
                    f'decimeters ** {dimension * volume_power}/(moles ** {volume_power} * seconds)')
                return converted_quantity.magnitude / (N_A ** volume_power), dimension, False
            else:
                converted_quantity = converted_quantity.convert(f'decimeters ** {dimension * volume_power}/seconds')
                return converted_quantity.magnitude, dimension, False
        except Exception as e:
            simlog.error(str(e) + '\n' +
                         f'Problem converting rate {quantity} \n'
                         f'Is the rate in the form [volume]**{volume_power}/[time]?')
    else:
        return quantity, dimension, False


def convert_counts(quantity, volume, dimension):
    """
        This function converts the counts from the users given unit to MobsPy standard units. It also converts
        concentrations into counts

        :param quantity: (int, float, Quantity) If it is a quantity object convert, otherwise it remains the same
        :param volume: (int, float) volume in liters (converted beforehand)
        :param dimension: (int) model's dimension (1D, 2D, 3D, ... )

        :return: converted_quantity (int, float) = converted unit into MobsPy standard units
    """
    try:
        converted_quantity = deepcopy(quantity)
    except:
        converted_quantity = quantity

    if isinstance(quantity, Quantity):

        if '[length]' not in str(quantity.dimensionality) and '[substance]' not in quantity.dimensionality\
                and str(quantity.dimensionality) != "dimensionless":
            simlog.error(f'The assigned quantity {quantity} is neither a count or concentration')
        elif str(quantity) == "dimensionless":
            return quantity.magnitude

        try:
            if '[substance]' in quantity.dimensionality:
                if '[length]' in str(quantity.dimensionality):
                    converted_quantity = converted_quantity.convert(f'moles/(decimeter ** {dimension})')
                    converted_quantity = converted_quantity * volume

                converted_quantity = converted_quantity.magnitude * N_A
            else:
                if '[length]' in str(quantity.dimensionality):
                    converted_quantity = converted_quantity.convert(f'1/(decimeter ** {dimension})')
                    converted_quantity = converted_quantity * volume
                converted_quantity = converted_quantity.magnitude
        except Exception as e:
            simlog.error(str(e) + '\n' +
                         f'Problem converting rate {quantity} \n'
                         f'Is it really a count or concentration?')
    return converted_quantity


def check_dimension(dimension, value, error_context=False):
    """
        Checks for dimension consistency. It "stores" the first dimension it was given by returning it

        :param dimension: (int) model's dimension (1D, 2D, 3D ...)
        :param value: (int) dimension value being analysed
        :param error_context: (bool or str) context of the error if dimensions are not consistent

        :raise simlog.error: If dimensions are not consistent through the given units (units in 1D with 2D mixed)

        :return: dimension (int) = model's dimension (1D, 2D, 3D ...)
    """
    if dimension is None:
        dimension = int(value)
    else:
        if dimension != int(value):
            message = 'The dimensions are not consistent. There are at least two units given for different ' \
                      'dimension models.'
            if error_context:
                message = message + '\n ' + error_context
            simlog.error(message)
    return dimension


def extract_length_dimension(unit_string, dimension, reaction_order=None, context=False):
    """
        Extracts the volume dimension from a Quantity object from Pint

        :param unit_string: (str) unit in str format
        :param dimension: (int) model's dimension (1D, 2D, 3D ...)
        :param reaction_order: (int) number of reactants in a reaction (for dimensional consistency in rates)
        :param context: (bool or str) context of the error if dimensions are not consistent
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
            dimension = check_dimension(dimension, temp_list[position + 2], context)
        else:
            temp_int = int(int(temp_list[position + 2]) / (reaction_order - 1))
            dimension = check_dimension(dimension, temp_int, context)
    else:
        dimension = check_dimension(dimension, 1, context)

    return dimension


def convert_volume(volume, dimension):
    """
        Converts volume to decimetre**dimension

        :param volume: (int, float, Quantity) volume used in simulation

        :return: the converted volume in MobsPy units
    """
    if isinstance(volume, Quantity):
        dimension = extract_length_dimension(str(volume.dimensionality), dimension)
        volume = volume.convert(f'decimeter ** {dimension}').magnitude
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
            return time.convert('second').magnitude
    else:
        return time


def time_convert_to_other_unit(time, other_unit):
    """
        Converts time into seconds

        :param time: (int, float, Quantity) any time used
    """
    if isinstance(time, Quantity):
        if str(time.dimensionality) == '[time]':
            return time.convert(other_unit).magnitude
    else:
        return time


if __name__ == '__main__':
    pass
