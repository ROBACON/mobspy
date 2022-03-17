from pint import UnitRegistry, Quantity
from scipy.constants import N_A
from copy import deepcopy
import os, sys
import simulation_logging.log_scripts as simlog


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
                converted_quantity.ito(f'liters ** {volume_power}/(moles ** {reaction_order} * seconds)')
                return converted_quantity.magnitude * (N_A ** reaction_order)
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


def convert_counts(quantity):
    if isinstance(quantity, Quantity):
        if str(quantity.units) == 'mole':
            return quantity.magnitude * N_A
        else:
            return quantity.magnitude
    return quantity


if __name__ == '__main__':
    pass
