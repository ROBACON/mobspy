from mobspy import *

'''
    Here is a small introductory model to MobsPy 
    We will go through some tutorial examples to get started
    Let us start with the basic model
'''

# We start by declaring a set of base-species to be used for modeling
# Both reactants and products needs to be defined
A, B, C, D = BaseSpecies(4)

# Reactions are defined using the >> operator and rates are assigned using the brackets
# The call operator in a species is used to assign counts to them
A(200) + B(100) >> 2*C + D [0.001]

# Chose which species to simulate with the Simulation Object
# And set parameters using the . operator in the resulting object
# Standard simulation duration is 60 seconds
# Species are named according to their variable name. You can also name with the .name method
MySim = Simulation(A | B | C | D)
MySim.save_data = False
MySim.run()

'''
    Now a model with reaction inheritance
    In MobsPy meta-species inherit from the species used to create them
    Either by multiplying them or by using the New call
'''

# Let's define the base species and the deriving species
Base, S = BaseSpecies(2)
A, B, C = New(Base, 3)

# Now we define the reaction
Base + S >> 2*Base [1]
# Here all meta-species constructed from Base will replace base on the meta-reaction
# All possibilities will be constructed for all combinations of inheritance
# This means this meta reaction defines the following reactions
# A + S >> 2*A
# B + S >> 2*B
# C + S >> 2*C

# Since the Base was only defined for the reaction inheritance we do not pass it to the simulation
MySim = Simulation(A | B | C | S)
# Finally we compile the model to check it
print(MySim.compile())


'''
    Now a simple model to explain the query system
'''

# First we define 2 BaseSpecies
Age, Color = BaseSpecies(2)
# We add characteristics to them using the dot operator
Age.young, Age.old
Color.blue, Color.red, Color.yellow
# Obs: Characteristics can also be added dynamically when used for the first time in reactions

# We multiply Age and Color
Dummy = Age*Color
# Now Dummy contains the species Dummy.blue.young, Dummy.blue.old, ... (all possible combinations)

# So we say only old Dummy can die
Dummy.old >> Zero [1]
# Zero is a BaseSpecies that means to nothing
# With this syntax only the Dummies species with old in their characteristics will be a part of this reaction
# The query also works for inherited species

MySim = Simulation(Dummy)
print(MySim.compile())

'''
    Restraints
    MobsPy has some restraints for modeling
    Species can only share characteristics thorough inheritance (multiplication or new)
    If one adds the same characteristic to two different species the compiler will show an error
    In this case inheritance can easily correct this model
'''
try:
    A, B = BaseSpecies(2)
    A.a, A.b
    C = New(B)
    C.a, C.b
    MySim = Simulation(A | C)
    MySim.level = 0
    MySim.compile()
except SystemExit:
    pass













