from mobspy import *
'''
    MobsPy allows Meta-Species to inherit reactions from other Meta-Species
    thus reducing the number of reactions needed to be programmed. 
    In this example we will go through the basics of inheritance.
'''

# We start by defining the BaseSpecies we will use.
Duplicator, Mortal, Eater, Food = BaseSpecies(4)

# Next, we define the reactions for each of the BaseSpecies.
# Zero represents no species.
Duplicator >> 2*Duplicator [1]
Eater + Food >> Eater [1]
Mortal >> Zero [1]

# Multiplication:
#   To define a bacteria species, we will use multiplication.
#   In MobsPy multiplication allows for the inheritance of reactions
#   from the factor-species.
# For example, we would like to define a Bacteria species that is also a Duplicator,
# an Eater, and a Mortal.
# This is done via
Bacteria = Duplicator*Eater*Mortal

# Another way to inherit is using the New call.
# Now New_Bacteria inherits from Bacteria, which means it is also a Duplicator, an Eater, and a Mortal.
New_Bacteria = New(Bacteria)

# The New operator may be used to generate several species at once as
# in this example.
Glucose, Amino = New(Food, 2)

# MobsPy will create all possible combinations for inheritance when defining the model,
# which means that both Bacteria and New_Bacteria will eat Glucose and Amino,
# although this was defined after the Bacteria species.

# Finally, we define the species we want to simulate.
# Since Duplicator, Eater, Mortal, and Food were just defined to construct the other meta-species,
# we don't pass them to Simulation.
MySim = Simulation(Bacteria | New_Bacteria | Glucose | Amino)

# Compile returns a model string that defines the resulting model.
# This allows one to see all the reactions, species, mappings, and parameters.
# The result is available in 2_Reaction_Inheritance_Model.txt
print(MySim.compile())

