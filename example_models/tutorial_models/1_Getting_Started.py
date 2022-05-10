"""
    To install MobsPy just use pip:
    pip install mobspy

    Another option is to download the github and use
    pip install -e 'downloaded file location'
"""

# Import MobsPy
from mobspy import *
"""
    To start of let us simulate a simple A + B >> 2*C + D reaction
"""

# To start of we define the meta-species that will be used in the reaction
# We use the BaseSpecies constructor with the number of desired meta-species
A, B, C, D = BaseSpecies(4)

# Now we write the reaction. It follows the CRN syntax
# With the >> operator being the distinction between reactants and products
# The counts are assigned using the call operator
# The rates are assigned using the [] (getitem) operator
# If a number is given as a rate, we assume mass action kinetics
A(200) + B(100) >> 2*C + D [0.1]

# Finally we construct the Simulation object by giving the meta-species we want to simulate
# Species are named according to their variable
My_Sim = Simulation(A | B | C | D)

# To configure the simulation just use the dot operator
# To check the names there is a read_me in the parameter folder
My_Sim.save_data = False
My_Sim.duration = 5
My_Sim.volume = 10
My_Sim.simulation_method = 'stochastic'
My_Sim.run()