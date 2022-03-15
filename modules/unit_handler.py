from pint import UnitRegistry, Quantity
from scipy.constants import N_A
from copy import deepcopy
import os, sys
import simulation_logging.log_scripts as simlog


def convert_rate(quantity, reaction_order):
    converted_quantity = deepcopy(quantity)
    if isinstance(quantity, Quantity):
        try:
            if '[substance]' in quantity.dimensionality:
                converted_quantity.ito_base_units()
                converted_quantity.ito(f'liters ** {reaction_order}/(moles ** {reaction_order} * seconds)')
                converted_quantity = converted_quantity.magnitude * (N_A ** reaction_order)
            else:
                converted_quantity.ito_base_units()
                converted_quantity.ito(f'liters ** {reaction_order}/seconds')
                converted_quantity = converted_quantity.magnitude
        except Exception as e:
            print(e)
            simlog.error(f'Problem converting rate {quantity}')
        return converted_quantity
    else:
        return quantity


def convert_counts(quantity):
    if isinstance(quantity, Quantity):
        if str(quantity.units) == 'mole':
            return quantity.magnitude * N_A
    return quantity.magnitude


if __name__ == '__main__':
    pass
