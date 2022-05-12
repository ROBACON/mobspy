from mobspy import *

"""
    This example model was designed to talk about looping through species
    Since there are two ways, we leave it for the user to pick the best suited one
"""

"""
    In this model we loop through species in two different ways
"""


def first_loop():
    # Define a base species
    # And an inheritor through the New operator
    DNAStrains = BaseSpecies()
    DNACounter = New(DNAStrains, 1)

    # Here we loop through the species by assigning different characteristics
    for i in range(10):
        DNACounter.c(f"dna_{i}")
    print(DNACounter.get_characteristics())


"""
    The New Function can also receive a name instead of a number
    The name allows one to name the species with the string instead of it's variable
"""


def second_loop():
    DNAStrains = BaseSpecies()

    # In this model we use the New operator to create different meta-species
    # And store them in the DNAList variable using the | operator
    # Here since there are no variables to store each meta-species
    # the new function receives a string which will be the name of the meta-species
    # When it receives a string it only returns one species with the string as it's name
    DNAList = None
    for i in range(10):
        DummyDNA = New(DNAStrains,f"DNA_{i}")
        if DNAList is None:
            DNAList = DummyDNA
        else:
            DNAList = DummyDNA | DNAList

    # One can then loop through the DNAList using an iterator like a for loop
    for DNA in DNAList:
        print(DNA)


if __name__ == '__main__':
    first_loop()
    second_loop()

