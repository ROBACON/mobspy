import sys, os
from mobspy import *

"""
    Here we have a model with a Phage, Donor and Receptor with resources being considered
    Both the Donor and Receptor can reproduce with resources
    The Donor can produce Phages with resources
    Only the Receptor can be infected by the phages
"""
# Variable declaration
Resource, Mortal, Infectible = BaseSpecies()

# Resource Definition
AA, Glu = New(Resource)
AA(100), Glu(100)

# Donor and Phage Creation
Mortal >> Zero [0.1]
Donor, Phage = New(Mortal)
Donor(100)

# Here both Resources can react with the Donor due to inheritance
dup_rate = lambda _,resource: 0.2 if resource.is_a(AA) else 0.1
Donor + Resource >> 2*Donor [dup_rate]
Donor + Resource >> Donor + Resource + Phage [0.1]
Infectible.low_inf >> Infectible.high_inf [0.1]
Receptor = Mortal*Infectible
inf_rate = lambda receptor: 0.2 if receptor.high_inf else 0.1
Receptor.not_infected + Phage >> Receptor.early_infection [inf_rate]
Receptor.early_infection >> Receptor.late_infection [0.1]
# We separate the reproduction reactions, because Receptor will give birth to a low_inf one
Receptor + Resource >> Receptor.low_inf + Receptor [dup_rate]

MySim = Simulation(Donor | Receptor | Phage | AA | Glu)
print(MySim.compile())
