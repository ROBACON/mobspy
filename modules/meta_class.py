import copy

import meta_class_utils


class Model:

    # Basic storage of defined variables for model
    Characteristics = set()
    Reactions_object_set = set()
    Species_set = set()
    Species_from_characteristic = {}
    Species_from_object = {}
    Last_rate = None

    # Model essence for Copasi sbml file
    Species = {}
    Mappings = {}
    Events = {}
    # Reactions = {}

    @classmethod
    def override_get_item(cls, object_to_return, item):
        cls.Last_rate = item
        return object_to_return

    @classmethod
    def compile(cls, *args):
        # TODO no events for now

        if type(args[0]) == list:
            args = args[0]

        '''
            Here we construct the species for Copasi using their objects
            
            First we extract all the species from all the given objects and add them to a set
            and for the species dictionary for the model
            
            And construct the Species_from_characteristic dictionary 
            This dictionary returns the species with an specific characteristic using it's key
            
            Also we construct the Species_from_object dictionary
            This dictionary returns the species from an specific object using the object as key
        '''
        for spe in args:
            lists_of_characteristics = meta_class_utils.extract_characteristics(spe)

            # s1 adds to species, s2
            s1, s2, s3 = meta_class_utils.combine_characteristics(spe, lists_of_characteristics)
            
            cls.Species_set = cls.Species_set.union(s1)
            meta_class_utils.unite_dictionaries(cls.Species_from_characteristic, s2)
            meta_class_utils.unite_dictionaries(cls.Species_from_object, s3)

        cls.Species = meta_class_utils.turn_set_into_0_value_dict(cls.Species_set)


class Reactions:

    @staticmethod
    def __create_reactants_string(list_of_reactants):
        reaction_string = ''
        for i, r in enumerate(list_of_reactants):
            if r[2] > 1:
                reaction_string += str(r[2]) + '*' + str(r[0])
            else:
                reaction_string += str(r[0])

            for j, characteristic in enumerate(r[1]):
                if j == 0:
                    reaction_string += '.'
                else:
                    reaction_string += '_'
                reaction_string += characteristic

            if i != len(list_of_reactants) - 1:
                reaction_string += ' '
                reaction_string += '+'
                reaction_string += ' '

        return reaction_string

    def __str__(self):
        reaction_string = ''
        return self.__create_reactants_string(self.reactants) + \
               ' -> ' + self.__create_reactants_string(self.products)

    def __getitem__(self, item):
        return Model.override_get_item(self, item)

    def __init__(self, reactants, products):
        self.reactants = reactants
        self.products = products

        self.rate = Model.Last_rate
        if Model.Last_rate is not None:
            Model.Last_rate = None

        Model.Reactions_object_set.add(self)


class Reactants:

    def __getitem__(self, item):
        return Model.override_get_item(self, item)

    # TODO change to dictionary
    def __init__(self, object_reference, characteristics, stoichiometry=1):
        self.list_of_reactants = [{'object': object_reference, 'characteristics': characteristics, 'stoichiometry': stoichiometry}]

    def __add__(self, other):
        self.list_of_reactants += other.list_of_reactants
        return self

    def __rshift__(self, other):

        if isinstance(other, Species):
            p = Reactants(other, '')
        else:
            p = other

        reaction = Reactions(self.list_of_reactants, p.list_of_reactants)
        return reaction


# TODO not being used yet
class Characteristics:

    def __getattr__(self, item):
        pass

    def __init__(self):
        self.characteristics = set()


# TODO Define rshift for species
class Species:

    entity_counter = 0

    def __str__(self):
        return self._name

    # Creation of reactions using entities ########################
    def __getitem__(self, item):
        return Model.override_get_item(self, item)

    def __rmul__(self, stoichiometry):
        if type(stoichiometry) == int:
            r = Reactants(self, '', stoichiometry)
        else:
            raise ValueError('Stoichiometry can only be an int')

    def __add__(self, other):
        r1 = Reactants(self, '')
        if isinstance(other, Reactants):
            r2 = other
        else:
            r2 = Reactants(other, '')
        return r1 + r2

    def __radd__(self, other):
        Species.__add__(self, other)

    # Adding characteristics to an Entity or Property ####################
    def __getattr__(self, item):

        characteristics_from_references = meta_class_utils.unite_characteristics(self.species_references)

        # TODO Redundant

    # Creation of other objects through multiplication ##################
    def __mul__(self, other):

        self.update_entity_counter()
        name = 'E' + str(self.entity_counter)
        new_entity = Species(name)
        new_entity.species_references = meta_class_utils.combine_references(self, other)
        new_entity.species_references.add(new_entity)

        return new_entity

    def __init__(self, name):
        self._name = name
        self.characteristics = set()
        self.species_references = {self}

    def name(self, name):
        self._name = name

    def get_name(self):
        return self._name


# Property Call to return several properties as called
def create_properties(number_of_properties=1):

    to_return = []
    for i in range(number_of_properties):
        Species.entity_counter += 1
        name = 'P' + str(Species.entity_counter)
        to_return.append(Species(name))

    if number_of_properties == 1:
        return to_return[0]
    else:
        return tuple(to_return)


if __name__ == '__main__':

    A, B, C, D = create_properties(4)
    A.young >> A.old [2.1]
    B.live >> B.dead  [3]
    C.sad >> C.happy [4]
    Human = A*B*C
    # Model = Human | A
    # Model.simulate()


