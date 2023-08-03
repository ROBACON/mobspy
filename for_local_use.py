from mobspy import *
import numpy as np
import basico
import matplotlib.pyplot as plt
import os
import inspect
from mobspy.sbml_simulator.builder import build
import timeit


if __name__ == '__main__':

    # Test script
    x = MobsPyExpression('A', None, dimension=3,
                         count_in_model=True,
                         concentration_in_model=False,
                         count_in_expression=False,
                         concentration_in_expression=False)

    r = x ** 2.8

























