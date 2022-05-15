from mobspy import *
"""
    Here we discuss another complex MobsPy feature the .c query.
    This method queries over the value of a stored variable instead of it's name.
"""

# We start by defining the base species that we will use.
Color = BaseSpecies(1)

# Now we define a two lists of colors.
light = ['blue', 'red', 'yellow', 'pink']
dark = ['dark_blue', 'dark_red', 'dark_yellow', 'dark_pink']

# In the loop bellow we can see how the .c works:
#   We want to design transition reactions from the light colors to their dark versions.
#   Here, the variables lc and dc loop through the light colors and dark colors list.
#   At each loop iteration, lc is assigned the value of a color and
#   dc is set to the dark version of that color.
#   The .c method is then used to query for the string stored inside both variables.
#   The .c query can be used together with a standard query. The order does not mater.

for lc, dc in zip(light, dark):
    # With the .c method we now query over the characteristic value inside lc and dc.
    Color.c(lc) >> Color.c(dc) [1]

MySim = Simulation(Color)
print(MySim.compile())

