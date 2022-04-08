import mobspy.modules.meta_class_utils as mcu
from mobspy.modules.order_operators import *
from mobspy.modules.function_rate_code import *
import mobspy.modules.reaction_construction_nb as rc
import mobspy.modules.function_rate_code as frc
from pint import Quantity


# Easter Egg: I finished the first version on a sunday at the BnF in Paris
# If anyone is reading this, I highly recommend you study there, it is quite a nice place
class Compiler:
    """
        This is the simulator call
        It receives the constructed objects and keeps track of the number of number of base species created
        It does so for automatic naming before compilation
        After compilation the species receive the variable name
    """

    # Basic storage of defined variables for Compiler
    entity_counter = 0
    reactions_set = set()
    species_string_dict = {}
    ref_characteristics_to_object = {}
    sbml_species_for_query = {}
    last_rate = None

    extra_species_list = []

    @classmethod
    def get_extra_species_list(cls):
        to_return = deepcopy(cls.extra_species_list)
        cls.extra_species_list = []
        return to_return

    @classmethod
    def query_str(cls, species, characteristics):
        result = '( '

        for species_string in cls.sbml_species_for_query:
            species_string_split = species_string.split('_dot_')
            species_name = species_string_split[0]
            species_string_split = species_string_split[1:]

            if species.get_name() == species_name:
                if all([char in species_string_split for char in characteristics]):
                    result += species_string + ' + '
                    cls.extra_species_list.append(species_string)

        if result == '( ':
            simlog.error('Query did not return anything \n'
                         'Is ' + species.get_name() + ' on the simulation? \n'
                         ' ')
        else:
            result = result[:-3] + ' )'
            return result

    @classmethod
    def override_get_item(cls, object_to_return, item):
        """
            Due to priority in python the item is stored before the reaction
            So it is stored in the Compiler level and passed to the reaction object in the end
        """
        cls.last_rate = item
        return object_to_return

    @classmethod
    def name_all_involved_species(cls, names=None):
        """
            Name species automatically according to their variable names
        """
        if not names:
            simlog.error('Species must be named' +
                         'Please set ( names == globals() ) in the MobsPy constructor for automatic naming')

        for name, key in names.items():
            # If the symbol '$' is in there, the species has not yet been named
            if isinstance(key, Species) and '$' in key.get_name():
                key.name(name)

    @classmethod
    def compile(cls, species_to_simulate, volume=1, names=None, type_of_model='deterministic', verbose=True,
                default_order=Default):
        # Set in model condition
        Reacting_Species.in_model = False
        Species.in_model = False
        
        # Define dictionaries to return to avoid compatibility problems

        # TODO no events for now
        # If there is only a single species

        """
            Here we construct the species for Copasi using their objects
            
            First we extract all the species from all the given objects and add them to a set
            and for the species dictionary for the Compiler
            
            And construct the Species_from_characteristic dictionary 
            This dictionary returns the species with an specific characteristic using it's key
            
            Also we construct the Species_from_object dictionary
            This dictionary returns the species from an specific object using the object as key
        """

        # We name the species according to variables names for convenience
        cls.name_all_involved_species(names)
        names_used = set()
        for species in species_to_simulate:
            if '$' in species.get_name():
                simlog.error('An error has occurred and one of the species was not named')
            if species.get_name() in names_used:
                simlog.error(f'Names must be unique for all species\n' +
                             f'The name is {species.get_name()}')
            names_used.add(species.get_name())

        # Construct structures necessary for reactions
        cls.ref_characteristics_to_object = mcu.create_orthogonal_vector_structure(species_to_simulate)

        # Start by creating the Mappings for the SBML
        # Convert to user friendly format as well
        mappings_for_sbml = {}
        for spe_object in species_to_simulate:
            mappings_for_sbml[spe_object.get_name()] = []

        # List of Species objects
        for spe_object in species_to_simulate:
            list_of_definitions = []
            for reference in spe_object.get_references():
                list_of_definitions.append(reference.get_characteristics())
            cls.species_string_dict[spe_object] = mcu.create_species_strings(spe_object,
                                                                             list_of_definitions)
        # Set of reactions involved
        cls.reactions_set = set()
        for spe_object in species_to_simulate:
            for reference in spe_object.get_references():
                cls.reactions_set = cls.reactions_set.union(reference.get_reactions())

        # Dimension check. Here we check based on the units the area dimension
        dimension = None

        for reaction in cls.reactions_set:
            if isinstance(reaction.rate, Quantity):
                if uh.extract_length_dimension(str(reaction.rate.dimensionality), dimension):
                    dimension = uh.extract_length_dimension(str(reaction.rate.dimensionality), dimension)

        for spe_object in species_to_simulate:
            for count in spe_object.get_quantities():
                if isinstance(count['quantity'], Quantity):
                    if uh.extract_length_dimension(str(count['quantity'].dimensionality), dimension):
                        dimension = uh.extract_length_dimension(str(count['quantity'].dimensionality), dimension)

        volume = uh.convert_volume(volume, dimension)

        # assign them default order
        for reaction in cls.reactions_set:
            if reaction.order is None:
                reaction.order = default_order

        # Setting Species for SBML and 0 for counts
        species_for_sbml = {}
        for spe_object in species_to_simulate:
            for species_string in cls.species_string_dict[spe_object]:
                species_for_sbml[species_string] = 0
        cls.sbml_species_for_query = species_for_sbml

        # BaseSpecies the mappings for sbml
        for spe_object in species_to_simulate:
            for species_string in cls.species_string_dict[spe_object]:
                mappings_for_sbml[spe_object.get_name()].append(species_string.replace('_dot_', '.'))

        # Set initial counts for SBML
        for spe_object in species_to_simulate:
            for count in spe_object.get_quantities():

                if count['quantity'] == 0:
                    continue

                count_set = \
                    mcu.complete_characteristics_with_first_values(spe_object,
                                                                   count['characteristics'],
                                                                   cls.ref_characteristics_to_object)

                for species_string in species_for_sbml.keys():
                    species_set = mcu.extract_characteristics_from_string(species_string)
                    if species_set == count_set:
                        temp_count = uh.convert_counts(count['quantity'], volume, dimension)
                        if type(temp_count) == float and not type_of_model == 'deterministic':
                            simlog.warning('The stochastic simulation rounds floats to integers')
                            species_for_sbml[species_string] = int(temp_count)
                        else:
                            species_for_sbml[species_string] = temp_count
                        break

        # BaseSpecies reactions for SBML with theirs respective parameters and rates
        # What do I have so far
        # Species_String_Dict and a set of reaction objects in Reactions_Set
        reactions_for_sbml, parameters_for_sbml = rc.create_all_reactions(cls.reactions_set,
                                                                          cls.species_string_dict,
                                                                          cls.ref_characteristics_to_object,
                                                                          type_of_model, volume, dimension)
        parameters_for_sbml['volume'] = (volume, f'dimensionless')

        # O(n^2) reaction check for doubles
        for i, r1 in enumerate(reactions_for_sbml):
            for j, r2 in enumerate(reactions_for_sbml):
                if i == j:
                    continue

                if reactions_for_sbml[r1]['re'] == reactions_for_sbml[r2]['re'] \
                        and reactions_for_sbml[r1]['pr'] == reactions_for_sbml[r2]['pr']:
                    simlog.warning('The following reaction: \n' +
                                   f'{reactions_for_sbml[r1]} \n' +
                                   'Is doubled. Was that intentional? \n')

        if verbose:
            simlog.debug()
            simlog.debug('Species')
            for species in species_for_sbml:
                simlog.debug(species.replace('_dot_', '.') + ',' + str(species_for_sbml[species]))
            simlog.debug()

            simlog.debug('Mappings')
            for mapping in mappings_for_sbml:
                simlog.debug(str(mapping) + ' :')
                for element in mappings_for_sbml[mapping]:
                    simlog.debug(element)
            simlog.debug()

            simlog.debug('Parameters')
            for parameters in parameters_for_sbml:
                simlog.debug(parameters + ',' + str(parameters_for_sbml[parameters]))
            simlog.debug()

            simlog.debug('Reactions')
            for reaction in reactions_for_sbml:
                simlog.debug(reaction + ',' + str(reactions_for_sbml[reaction]).replace('_dot_', '.'))
        
        # Compilation finish, in model not needed anymore:
        Reacting_Species.in_model = False
        Species.in_model = False
        return species_for_sbml, reactions_for_sbml, parameters_for_sbml, mappings_for_sbml


class Reactions:
    """
        This is the Reaction class. It contains the reactants, products, rate and order
        Reactions are created with the >> operator. Each reaction is stored in all involved objects
    """

    @staticmethod
    def __create_reactants_string(list_of_reactants):
        """
            Just a simple way to simlog.debug reactions for debbuging
            Not relevant for simulation
        """
        reaction_string = ''
        for i, r in enumerate(list_of_reactants):
            if r['stoichiometry'] > 1:
                reaction_string += str(r['stoichiometry']) + '*' + str(r['object'])
            else:
                reaction_string += str(r['object'])

            reaction_string += '.' + '.'.join(r['characteristics'])

            if i != len(list_of_reactants) - 1:
                reaction_string += ' '
                reaction_string += '+'
                reaction_string += ' '

        return reaction_string

    def __str__(self):
        return self.__create_reactants_string(self.reactants) + \
               ' -> ' + self.__create_reactants_string(self.products)

    def __getitem__(self, item):
        return Compiler.override_get_item(self, item)

    def __init__(self, reactants, products):
        self.reactants = reactants
        self.products = products

        # Assign default order
        self.order = None

        self.rate = Compiler.last_rate
        if Compiler.last_rate is not None:
            Compiler.last_rate = None

        # Here we extract all involved objects to pact them in a set
        # This is done to find the reactions associated with the species when the Compiler is started
        for reactant in reactants:
            reactant['object'].add_reaction(self)
        for product in products:
            product['object'].add_reaction(self)

    def set_rate(self, rate):
        self.rate = rate


class Reacting_Species:
    """
        If a species is involved in a reaction this object is created
        It contains the species object reference, the characteristics involved in the reaction and finally the stoichiometry
        It adds up by appending elements to a list
    """
    # String querying inside the model
    in_model = False

    def __str__(self):

        if not self.in_model:
            if len(self.list_of_reactants) == 1:
                to_return = str(self.list_of_reactants[0]['object'])
                for cha in self.list_of_reactants[0]['characteristics']:
                    to_return += '.' + cha
                return to_return
            else:
                return self.list_of_reactants
        else:
            if len(self.list_of_reactants) == 1:
                return Compiler.query_str(self.list_of_reactants[0]['object'], self.list_of_reactants[0]['characteristics'])
            else:
                simlog.error('Query not accepted. Only one species to str at a time')
    
    # Labels and value function implementation
    def c(self, item):
        return self.__getattr__(item)

    def label(self, label):
        if len(self.list_of_reactants) == 1:
            self.list_of_reactants[0]['label'] = label
        else:
            simlog.error('Error assigning label to reactant')
        return self

    def __getitem__(self, item):
        return Compiler.override_get_item(self, item)

    def __init__(self, object_reference, characteristics, stoichiometry=1, label=None):
        if object_reference.get_name() == 'S0' and characteristics == set():
            self.list_of_reactants = []
        else:
            self.list_of_reactants = [{'object': object_reference, 'characteristics': characteristics,
                                       'stoichiometry': stoichiometry, 'label': label}]

    def __rmul__(self, stoichiometry):
        if type(stoichiometry) == int:
            self.list_of_reactants[0]['stoichiometry'] = stoichiometry
        else:
            simlog.error(f'Stoichiometry can only be an int - Received {stoichiometry}')
        return self

    def __add__(self, other):
        if isinstance(other, Species):
            other = Reacting_Species(other, set())
        self.list_of_reactants += other.list_of_reactants
        return self

    def __rshift__(self, other):
        if isinstance(other, Species):
            p = Reacting_Species(other, set())
        else:
            p = other

        reaction = Reactions(self.list_of_reactants, p.list_of_reactants)
        return reaction

    def __call__(self, quantity):
        if type(quantity) == int or type(quantity) == float or isinstance(quantity, Quantity):
            if len(self.list_of_reactants) != 1:
                simlog.error('Assignment used incorrectly')
            species_object = self.list_of_reactants[0]['object']
            characteristics = self.list_of_reactants[0]['characteristics']
            species_object.add_quantities(characteristics, quantity)
        else:
            simlog.error('Reactant_species does not support this type of call')
        return self

    def __getattr__(self, characteristic):
        for reactant in self.list_of_reactants:

            species_object = reactant['object']
            characteristics_from_references = mcu.unite_characteristics(species_object.get_references())

            if characteristic not in characteristics_from_references:
                species_object.add_characteristic(characteristic)

            reactant['characteristics'].add(characteristic)

        return self


class ParallelSpecies:
    """
        Class to simulate several species at the same time
        It is basically a list of species
        Species become this using the | operator
    """
    def __init__(self, list_of_species):
        self.list_of_species = list_of_species

    def append(self, species):
        if not isinstance(species, Species):
            simlog.error('Only Species can be appended')
        self.list_of_species.append(species)

    def __getitem__(self, item):
        return self.list_of_species[item]

    def __str__(self):
        to_return = []
        for spe in self.list_of_species:
            to_return.append(spe.get_name())
        return str(to_return)

    def __or__(self, other):
        if isinstance(other, Species):
            self.append(other)
            return self
        elif isinstance(other, ParallelSpecies):
            return ParallelSpecies(self.list_of_species + other.list_of_species)
        else:
            simlog.error('Operator must only be used in Species on ParallelSpecies')

    def __iter__(self):
        for spe in self.list_of_species:
            yield spe


class Species:
    """
        Fundamental class
        Contains the characteristics, the species name, the reactions it is involved in
        and finally the other species it references
        Objects store all the basic information necessary to create an SBML file and construct a model
        So we construct the species through reactions and __getattr__ to form a model
    """
    # String querying inside the model
    in_model = False

    def __str__(self):
        if not self.in_model:
            return self._name
        else:
            return Compiler.query_str(self, set())

    # Def c to get the value
    def c(self, item):
        return self.__getattr__(item)

    # Def labels
    def label(self, label):
        return Reacting_Species(self, set(), label=label)

    # Get data from species for debugging Compiler
    def show_reactions(self):
        simlog.debug(str(self) + '_dot_')
        for reference in self._references:
            for reaction in reference.get_reactions():
                simlog.debug(reaction)

    def show_characteristics(self):
        simlog.debug(str(self) + ' has the following characteristics referenced:')
        for i, reference in enumerate(self.get_references()):
            if reference.get_characteristics():
                simlog.debug(str(reference) + ': ', end='')
                reference.simlog.debug_characteristics()

    def show_references(self):
        simlog.debug(str(self) + '_dot_')
        simlog.debug('{', end='')
        for i, reference in enumerate(self.get_references()):
            if reference.get_characteristics():
                simlog.debug(' ' + str(reference) + ' ', end='')
        simlog.debug('}')

    def show_quantities(self):
        simlog.debug(self._species_counts)

    # Creation of ParallelSpecies For Simulation ##################
    def __or__(self, other):
        if isinstance(other, ParallelSpecies):
            other.append(self)
            return other
        else:
            return ParallelSpecies([self, other])

    def __iter__(self):
        yield self

    # Creation of reactions using entities ########################
    def __getitem__(self, item):
        return Compiler.override_get_item(self, item)

    def __rmul__(self, stoichiometry):
        if type(stoichiometry) == int:
            r = Reacting_Species(self, set(), stoichiometry)
        else:
            simlog.error(f'Stoichiometry can only be an int - Received {stoichiometry}')
        return r

    def __add__(self, other):
        r1 = Reacting_Species(self, set())
        if isinstance(other, Reacting_Species):
            r2 = other
        else:
            r2 = Reacting_Species(other, set())
        return r1 + r2

    def __radd__(self, other):
        Species.__add__(self, other)

    def __rshift__(self, other):
        myself = Reacting_Species(self, set())

        if isinstance(other, Species):
            p = Reacting_Species(other, set())
        elif other == 0:
            exit()
        else:
            p = other

        reaction = Reactions(myself.list_of_reactants, p.list_of_reactants)
        return reaction

    # Adding characteristics to an Entity or Property ####################
    def __getattr__(self, characteristic):
        """
            characteristic: string of characteristic
            We add characteristics using the following manner Species.characteristic
            If it is the first time is called we add it to the set of characteristics
            Second time just creates a Reacting_Species for reaction construction

            IMPORTANT - DON'T USE DEEPCOPY - IT DOES NOT WORK WITH __getattr__
        """

        characteristics_from_references = mcu.unite_characteristics(self.get_references())
        characteristics = {characteristic}

        if characteristic not in characteristics_from_references:

            if len(self.get_characteristics()) == 0:
                self.first_characteristic = characteristic

            self.add_characteristic(characteristic)

        return Reacting_Species(self, characteristics)

    # Adding counts to species
    def __call__(self, quantity):
        # It's called quantity because of the original design
        self.in_model = False
        if type(quantity) == int or type(quantity) == float or isinstance(quantity, Quantity):
            self.add_quantities('std$', quantity)
        elif isinstance(quantity, frc.Specific_Species_Operator):
            for cha in str(quantity).split('_dot_')[1:]:
                if cha in self._characteristics:
                    return cha
            simlog.error(f'{self._name} contains no characteristics in {quantity}')
        elif type(quantity) == str:
            for cha in str(quantity).split('.')[1:]:
                if cha in self._characteristics:
                    return cha
            simlog.error(f'{self._name} contains no characteristics in {quantity}')
        self.in_model = True
        return self

    def add_quantities(self, characteristics, quantity):
        '''
        :param characteristics: characteristics of the species to be set
        :param quantity: counts of that specific species

            We do it like this because Python can't accept sets as dictionary keys
            Set the quantity of an specific string of species
        '''
        already_in = False
        for e in self._species_counts:
            if characteristics == e['characteristics']:
                e['quantity'] = quantity
                already_in = True
        if not already_in:
            self._species_counts.append({'characteristics': characteristics, 'quantity': quantity})

    def reset_quantities(self):
        self._species_counts = []

    def get_quantities(self):
        return self._species_counts

    # Creation of other objects through multiplication ##################
    def __mul__(self, other):
        """
            Multiplications are used to construct more complex species
            We do not concatenate characteristics. This keeps the structure intact for reaction construction
            Instead we combine the sets of references. Every processes references itself
        """
        Compiler.entity_counter += 1
        name = 'N$' + str(Compiler.entity_counter)
        new_entity = Species(name)
        new_entity.set_references(mcu.combine_references(self, other))
        new_entity.add_reference(new_entity)

        mcu.check_orthogonality_between_references(new_entity.get_references())

        return new_entity

    def __init__(self, name):
        self._name = name
        self._characteristics = set()
        self._references = {self}

        # This is necessary for the empty objects generated when we perform multiplication with more than 2 Properties
        self.first_characteristic = None

        # Each object stores the reactions it is involved in
        self._reactions = set()

        # This will store the quantities relating to the species counts
        self._species_counts = []

    def name(self, name):
        self._name = name

    def get_name(self):
        return self._name

    def get_characteristics(self):
        return self._characteristics

    def add_characteristic(self, characteristic):
        self._characteristics.add(characteristic)

    def remove_characteristic(self, characteristic):
        self._characteristics.remove(characteristic)

    def print_characteristics(self):
        simlog.debug(self._characteristics)

    def get_references(self):
        return self._references

    def add_reference(self, reference):
        self._references.add(reference)

    def set_references(self, reference_set):
        self._references = reference_set

    def get_reactions(self):
        return self._reactions

    def set_reactions(self, reactions):
        self._reactions = reactions

    def add_reaction(self, reaction):
        self._reactions.add(reaction)


# Property Call to return several properties as called
def BaseSpecies(number_of_properties=1):
    """
        A function to create multiple species at once
        Reduce the number of lines
    """
    to_return = []
    for i in range(number_of_properties):
        Compiler.entity_counter += 1
        name = 'N$' + str(Compiler.entity_counter)
        to_return.append(Species(name))

    if number_of_properties == 1:
        return to_return[0]
    else:
        return tuple(to_return)


"""
    The Zero is defined here
    Used this for reactions that result in nothing or come from nothing
    The One is used to create new similar instances of a base species
    It's used to implement Copy
"""
__S0, __S1 = BaseSpecies(2)
__S0.name('S0')
__S1.name('S1')

Zero = __S0


def New(species, n=1):
    if type(n) == int:
        to_return = [species*__S1 for _ in range(n)]
        if n == 1:
            return to_return[0]
        else:
            return tuple(to_return)
    elif type(n) == str:
        to_return = species*__S1
        to_return.name(n)
        return to_return
    else:
        simlog.error('New only supports numbers for creating multiple species \n '
                     ' or exclusively string for single species name')


if __name__ == '__main__':
    pass
