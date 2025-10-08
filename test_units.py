"""Simple tests for units in MobsPy."""

import pytest
import yaml
from mobspy import BaseSpecies, Simulation, u


def test_def_by_str():
    with open("testconfig_units.yaml", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    A0 = u(config["A0"])
    rate_constant = u(config["rate_constant"])

    assert str(A0) == '1.0 / milliliter', str(A0)
    assert str(rate_constant) == '1.0 / minute', str(rate_constant)

    print('test_def_by_str: outside of function:', rate_constant)
    print('test_def_by_str: outside of function:', type(rate_constant))

    duration = 5 * u.min

    def rate_fn(a):
        print('test_def_by_str: in function:', rate_constant)
        print('test_def_by_str: in function:', type(rate_constant))

        assert str(A0) == '1.0 / milliliter', str(A0)
        assert str(rate_constant) == '1.0 / minute', str(rate_constant)

        return rate_constant

    A = BaseSpecies()
    A(A0)
    
    A >> 2 * A [rate_fn]

    MySim = Simulation(A)
    MySim.volume = 1 * u.mL
    MySim.compile()
    # MySim.run(duration=duration, unit_x=u.min, unit_y=1 / u.mL)    
        

def test_def_by_quantity():
    A0 = 1 / u.mL
    rate_constant = 1 / u.min
    
    assert str(A0) == '1.0 / milliliter', str(A0)
    assert str(rate_constant) == '1.0 / minute', str(rate_constant)

    print('test_def_by_quantity: outside of function:', rate_constant)
    print('test_def_by_quantity: outside of function:', type(rate_constant))

    duration = 5 * u.min

    def rate_fn(a):
        print('test_def_by_quantity: in function:', rate_constant)
        print('test_def_by_quantity: in function:', type(rate_constant))

        assert str(A0) == '1.0 / milliliter', str(A0)
        assert str(rate_constant) == '1.0 / minute', str(rate_constant)

        return rate_constant

    A = BaseSpecies()
    A(A0)
    
    A >> 2 * A [rate_fn]

    MySim = Simulation(A)
    MySim.volume = 1 * u.mL
    MySim.compile()
    # MySim.run(duration=duration, unit_x=u.min, unit_y=1 / u.mL)    


if __name__ == "__main__":
    test_def_by_str()
    test_def_by_quantity()