from mobspy import *
"""
    Here we discuss another complex MobsPy feature the .c query
    This method queries over the value of a stored variable instead of it's name
"""

# We start by defining the base species that we will use
Color = BaseSpecies(1)

# Now we define a two lists of colors
light = ['blue', 'red', 'yellow', 'pink']
dark = ['dark_blue', 'dark_red', 'dark_yellow', 'dark_pink']

# Now we loop through the colors to add the reactions
# In this model all colors change to their dark version with the reaction
for lc, dc in zip(light, dark):
    Color.c(lc) >> Color.c(dc) [1]

# With the .c method we query over the characteristic value inside lc and dc
# Otherwise MobsPy would query over lc and dc 4 times

MySim = Simulation(Color)
print(MySim.compile())

