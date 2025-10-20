"""Simple tests for multiline reactions."""

import pytest
import numpy as np
from mobspy import BaseSpecies, Simulation

def test_numpy():
    expected = """
Species
A,10.0
B,0

Mappings
A :
A
B :
B

Parameters
volume,1

Reactions
reaction_0,{'re': [(1, 'A')], 'pr': [(1, 'B')], 'kin': 'A * 1.0'}
"""

    A, B = BaseSpecies()
    params = np.array([10.0, 1.0])

    # use numpy to init params
    A(params[0])
    A >> B [params[1]]
    
    MySim = Simulation(A | B)
    MySim.level = -1
    result = MySim.compile()    
    assert result == expected

if __name__ == "__main__":
    test_numpy()