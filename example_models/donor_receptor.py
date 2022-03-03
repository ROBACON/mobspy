import sys, os
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *


# Variable declaration
Resource, Mortal, Infectible = BaseSpecies(3)

# Resource Definition
AA, Glu = New(Resource, 2)
AA(100)
Glu(100)

# Donor and Phage Creation
Mortal >> Zero [0.1]
Donor, Phage = New(Mortal, 2)
Donor(100)

Donor + Resource >> 2*Donor [lambda _,resource: 0.2 if IsReference(resource,AA) else 0.1]
Donor + Resource >> Donor + Resource + Phage [0.1]

# Receptor and infection
Infectible.low_inf >> Infectible.high_inf [0.1]
Receptor = Mortal*Infectible
Receptor.not_infected + Phage >> Receptor.early_infection [lambda receptor: 0.2 if receptor.high_inf else 0.1]
Receptor.early_infection >> Receptor.late_infection [0.1]
Receptor + Resource >> Receptor.young + Receptor [lambda _,resource: 0.2 if IsReference(resource,AA) else 0.1]

Simulation(Donor | Receptor | Phage | AA | Glu).compile()
