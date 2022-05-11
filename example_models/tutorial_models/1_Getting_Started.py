"""
    To install MobsPy just use pip:
    pip install mobspy

    Another option is to download the github and use
    pip install -e 'downloaded file location'
"""

# Import MobsPy
from mobspy import *
"""
    To start off let us simulate a simple A + B >> 2*C + D reaction
"""

# First, we define the meta-species that will be used in the reaction.
# We use the BaseSpecies constructor with the number of desired meta-species.
A, B, C, D = BaseSpecies(4)

# Now we write the reaction. It follows the CRN syntax
# >> operator:
#   With the >> operator separating reactants on the left and products
#   on the right.
# Initial counts:
#   Initial counts are assigned using the call operator.
#   For example A(200) initializes species A with count 200.
# Reaction rate:
#   Each reaction has a rate that is assigned using the [] operator.
#   If a float is given as a rate, we assume mass action kinetics with the
#   float being the reaction rate constant.
A(200) + B(100) >> 2*C + D [0.1]

# This will create a reaction A + B -> 2*C + D
# with a rate equal to [A]*[B]*0.1.

# We next construct the Simulation object by giving it the meta-species we want to simulate.
# Species are named according to their variables.
My_Sim = Simulation(A | B | C | D)

# To configure the simulation just use the dot operator with the
# respective attribute.
# The available attributes are documented in a read_me in the parameter folder.
My_Sim.save_data = False
My_Sim.duration = 5
My_Sim.volume = 10
My_Sim.simulation_method = 'stochastic'

# Finally the simulation is run and results
# printed via
My_Sim.run()
