import sys, os
remove = os.getcwd().split('/')[:-1]
abs_path = '/'.join(remove)
sys.path.append(os.path.abspath(abs_path))

from mobspy import *

# TODO fix _dot_ for errors

# Variable declaration
Resource, Death, Age = BaseSpecies(3)

# Resource Definition
AA = Resource*New
Glu = Resource*New
AA(100)
Glu(100)

# Donor and Phage Creation
Death >> Zero [0.1]
Donor = Death*New
Phage = Death*New
Donor(100)

Donor + Resource >> 2*Donor [0.1]
Donor + Resource >> Donor + Resource + Phage [0.1]

# Receptor and infection
Age.young >> Age.old [0.1]
Receptor = Death*Age
Receptor.not_infected + Phage >> Receptor.early_infection [lambda rec: 0.2 if rec.old else 0.1]
Receptor.early_infection >> Receptor.late_infection [0.1]
Receptor + Resource >> Receptor.young + Receptor [0.1]

Simulation(Donor | Receptor | Phage | AA | Glu).compile()
Simulation.simulation_method = 'stochastic'
