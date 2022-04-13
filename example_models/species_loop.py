from mobspy import *

"""
    This example model was designed to talk about looping through species
    Since there are two ways, we leave it for the user to pick the best suited one
"""

"""
    In this model we loop through species in two different ways
    One by assigning a set of characteristics to a new dummy species
    And the other by creating a list of species
"""


def first_loop():
    DNAStrains = BaseSpecies()
    DNACounter = New(DNAStrains, 1)

    for i in range(10):
        DNACounter.c(f"dna_{i}")

    # Set of characteristics
    print(DNACounter.get_characteristics())


"""
    The New Function can also receive a name instead of a number
    The name allows one to name the species with the string instead of it's variable
"""


def second_loop():
    DNAStrains = BaseSpecies()

    DNAList = None
    for i in range(10):
        DummyDNA = New(DNAStrains,f"DNA_{i}")
        if DNAList is None:
            DNAList = DummyDNA
        else:
            DNAList = DummyDNA | DNAList

    for DNA in DNAList:
        print(DNA)


if __name__ == '__main__':
    first_loop()
    second_loop()

