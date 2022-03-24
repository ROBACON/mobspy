from pint import UnitRegistry, Quantity
from scipy.constants import N_A
from copy import deepcopy
import os, sys
import mobspy.simulation_logging.log_scripts as simlog


def convert_rate(quantity, reaction_order):

    volume_power = reaction_order - 1
    converted_quantity = deepcopy(quantity)
    if isinstance(quantity, Quantity):
        try:
            if volume_power <= 0:
                converted_quantity.ito_base_units()
                converted_quantity.ito(f'1/seconds')
                return converted_quantity.magnitude

            if '[substance]' in quantity.dimensionality:
                converted_quantity.ito_base_units()
                converted_quantity.ito(f'liters ** {volume_power}/(moles ** {volume_power} * seconds)')
                return converted_quantity.magnitude / (N_A ** volume_power)
            else:
                converted_quantity.ito_base_units()
                converted_quantity.ito(f'liters ** {volume_power}/seconds')
                return converted_quantity.magnitude
        except Exception as e:
            print(e)
            simlog.error(f'Problem converting rate {quantity} \n'
                         f'Is the rate in the form [volume] ** (order - 1)/[time]?')

    else:
        return quantity


def convert_counts(quantity, volume):
    converted_quantity = deepcopy(quantity)
    if isinstance(quantity, Quantity):
        try:
            if '[substance]' in quantity.dimensionality:
                converted_quantity.ito_base_units()
                if '[length] ** 3' in str(quantity.dimensionality):
                    converted_quantity.ito('moles/liter')
                    converted_quantity = converted_quantity*volume

                converted_quantity = converted_quantity.magnitude * N_A
            else:
                converted_quantity.ito_base_units()
                if '[length] ** 3' in str(quantity.dimensionality):
                    converted_quantity.ito('moles/liter')
                    converted_quantity = converted_quantity*volume
                converted_quantity = converted_quantity.magnitude

        except Exception as e:
            print(e)
            simlog.error(f'Problem converting rate {quantity} \n'
                         f'Is it really a count? Concentrations must be converted')
    return converted_quantity


if __name__ == '__main__':
    pass
