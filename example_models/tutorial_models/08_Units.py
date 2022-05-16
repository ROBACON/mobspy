from mobspy import *
"""
    This model was designed to explain how units can be added to a MobsPy project.
    The units are based on the pint module.
"""

# Defining base-species and give it some characteristics:
A = BaseSpecies(1)
A.a1, A.a2, A.a3, A.a4

# MobsPy adopts a Liter(decimeter)-Second-Count convention.
# We assume those units when a number is given without units.
# To add units use the u letter and the dot operator with the unit:
A.a1(1*u.mol)

# One can use concentrations to specify initial counts as well.
# They will be converted using the model volume:
A.a2(2*u.nanomolar)

# Rate constants can also use units.
# In this case the correctness of units for rate constants are checked:
# 1/time for one reactant, volume/time for 2 reactants, volume**2/time for 3 reactants ...
A.a1 + A.a2 >> A.a3 + A.a4 [1*u.meter**3/u.second]

# Functions that return rate constants can return units as well:
A >> A.a4 [lambda r1: 5/u.second if r1.a1 else 1/u.second]

# But functions that return strings have no support for units yet.
# They are implicitely supposed to return rates in terms of 1/second.

MySim = Simulation(A)

# To set the duration and the volume of the simulation with units we use:
MySim.volume = 1*u.nanolitre
MySim.duration = 3*u.hours

# The model compilation shows the values in MobsPy with standard units:
# Liter(decimeter)-Second-Count
print(MySim.compile())

"""
    A more advanced model for other unit types
"""

# In MobsPy one can also use 2-D concentrations and volumes.
# As an example one can refer to the For_The_Tress model in the application models.
# There, concentrations are expressed as Count/Area and rates are expressed in Area/Time.
# The volume would be an area in this case (Although the parameter is still named volume; it is a volume in a 2D-space).
# If an inconsistence in dimensions (e.g., mixture between 2-D and 3-D) is detected, the compiler will throw an error.

print()
print('Other Model')
print()

try:
    B = BaseSpecies(1)
    B(2*u.molar)
    MySim = Simulation(B)
    MySim.volume = 2*u.meter**2
    print(MySim.compile())
except SystemExit:
    pass



