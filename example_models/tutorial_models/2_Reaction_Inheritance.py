from mobspy import *
'''
    MobsPy allows Meta-Species to inherit reactions from other Meta-Species
    Thus reducing the number of reactions needed to be programmed 
    In this example we will go through the basics of inheritance
'''

# We start by defining the BaseSpecies we will use
Duplicator, Mortal, Eater, Food = BaseSpecies(4)

# Now we define the reactions for each of the BaseSpecies
# Zero represents no species
Duplicator >> 2*Duplicator [1]
Eater + Food >> Eater [1]
Mortal >> Zero [1]

# Now we define a bacteria
# In MobsPy multiplication allows for the inheritance of reactions
# Now a Bacteria is also a Duplicator, an Eater and a Mortal
# This means that the Bacteria will replace the species it inherits from to define a reaction
Bacteria = Duplicator*Eater*Mortal

# Another way to inherit is using the New call
# Now New_Bacteria inherits from Bacteria, which means it is also a Duplicator, an Eater and a Mortal
# And it will also replace them to define reactions
New_Bacteria = New(Bacteria)
Glucose, Amino = New(Food, 2)

# MobsPy will create all possible combinations for inheritance when defining the model
# Which means that both Bacteria and New_Bacteria will eat Glucose and Amino

# Finally we only add the species we want to simulate
# Since Duplicator, Eater, Mortal and Food were just defined to construct the other meta-species
# we don't pass them to the simulation
MySim = Simulation(Bacteria | New_Bacteria | Glucose | Amino)

# Compile returns a model string, it allows one to see all the reactions in the model,
# all species, all mappings and all parameters
# This resulting string is written on the 2_Reaction_Inheritance_Model.txt
print(MySim.compile())

