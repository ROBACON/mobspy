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

# In the loop bellow we can see how the .c works:
#   We want to design transition reactions from the light colors to their dark version
#   Here the variable lc and dc loop through the light colors and dark colors list
#   So at each iteration of the loop lc is granted the value of a color and
#   dc is given the dark version of that color
#   As lc is looping through the light list and dc is looping through the dark list (in their respective order)
#   Then the .c method is responsible for performing a query for the string stored inside the variable
#   Instead of the name of the variable like the previous queries
#
for lc, dc in zip(light, dark):
    Color.c(lc) >> Color.c(dc) [1]

# The .c query can be used with the standard query and the order does not mater

# With the .c method we query over the characteristic value inside lc and dc
# Otherwise MobsPy would query over lc and dc 4 times

MySim = Simulation(Color)
print(MySim.compile())

