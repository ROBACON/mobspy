import sys, os
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *


# Variable declaration
Resource, Mortal, Infectible = BaseSpecies(3)

# Resource Definition
AA, Glu = New(Resource, 2)
AA(100), Glu(100)

# Donor and Phage Creation
Mortal >> Zero [0.1]
Donor, Phage = New(Mortal, 2)
Donor(100)

dup_rate = lambda _,resource: 0.2 if resource.is_a(AA) else 0.1
Donor + Resource >> 2*Donor [dup_rate]
Donor + Resource >> Donor + Resource + Phage [0.1]
Infectible.low_inf >> Infectible.high_inf [0.1]
Receptor = Mortal*Infectible
inf_rate = lambda receptor: 0.2 if receptor.high_inf else 0.1
Receptor.not_infected + Phage >> Receptor.early_infection [inf_rate]
Receptor.early_infection >> Receptor.late_infection [0.1]
Receptor + Resource >> Receptor.young + Receptor [dup_rate]

Simulation(Donor | Receptor | Phage | AA | Glu).compile()
