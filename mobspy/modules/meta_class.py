"""
    The meta_class.py model is responsible for defining the meta-species, meta-reactions and other "meta-classes"
    It is also responsible for compiling the models into data structures that will be converted into an SBML file later
    Bear in mind that they are not actually Python meta-classes. The first design utilized this feature but know there
    are just regular classes
"""
import mobspy.modules.meta_class_utils as mcu
from mobspy.modules.order_operators import *
from mobspy.modules.function_rate_code import *
import mobspy.modules.reaction_construction_nb as rc
import mobspy.modules.function_rate_code as frc
from pint import Quantity
import mobspy.modules.event_functions as eh
from mobspy.modules.logic_operator_objects import *
import mobspy.modules.species_string_generator as ssg
from copy import deepcopy


# Easter Egg: I finished the first version on a sunday at the BnF in Paris
# If anyone is reading this, I highly recommend you study there, it is quite a nice place
class Compiler:
    """
        The compiler is responsible for constructing the model from the meta-species and meta-reactions given
        and checking if it is a valid MobsPy model

        :param entity_counter: (int) counts the number of defined meta-species (for having a unique placeholder name
        for species before the user or MobsPy name them)
        :param  last_rate: (float, callable, Quantity)  due to python operator priority the __getitem__
        operation executes in the beginning of a reaction definition, so the compiler stores it until it is needed
        :param ref_characteristics_to_object: (dict)  dictionary where the keys are the characteristics and the value
         is the meta-species that the characteristic was added directly to (no inheritance)
    """

    # Basic storage of defined variables for Compiler
    entity_counter = 0
    ref_characteristics_to_object = {}
    last_rate = None

    @classmethod
    def override_get_item(cls, object_to_return, item):
        """
            Due to priority in Python the item is stored before the reaction
            So it is stored in the Compiler level and passed to the reaction object in the end
            Important: the rate is used before the compilation level, it's used when the reaction has been completely
            defined

            :param object_to_return: (Species or Reacting) returns the object that the __getitem__ was performed on
            :param item: (int, float, callable, Quantity) stored reaction rate
        """
        cls.last_rate = item
        return object_to_return

    @classmethod
    def name_all_involved_species(cls, names=None):
        """
            Name species automatically according to their variable names

            :param names: (dict) dictionary with other name options - follows globals() format
        """
        if not names:
            simlog.error('Species must be named' +
                         'Please set ( names == globals() ) in the MobsPy constructor for automatic naming')

        for name, key in names.items():
            # If the symbol '$' is in there, the species has not yet been named
            if isinstance(key, Species) and '$' in key.get_name():
                key.name(name)

    @classmethod
    def compile(cls, meta_species_to_simulate, reactions_set, species_counts,
                volume=1, names=None, type_of_model='deterministic', verbose=True,
                default_order=Default, event_dictionary=None,
                continuous_sim=False, ending_condition=None):
        """
            Here we construct the species for Copasi using their objects

            It handles the construction of species_for_sbml, reactions_for_sbml, parameters_for_sbml,
            mappings_for_sbml, model_str which will be used to write the sbml_string for basiCO
            It also performs checks to see if the model is valid

            :param meta_species_to_simulate: (ParallelSpecies or Species) Species to generate the model
            :param volume: (int, flot) Simulation volume
            :param names: (dict) alternative naming for species - see name_all_involved_species
            :param type_of_model: (str) deterministic or stochastic
            :param verbose: (bool) print or not the model results by generating a model_str
            :param default_order: (Reaction Operator) Operator that decides the order of the
            meta-reaction ( or Default we have the round-robin) and how MobsPy will handle meta-species in
            the products that are not referenced in the reactants - see order_operators.py
            :param event_dictionary: (dict) dictionary with all the events to be added to the model
            :param continuous_sim: (bool) simulation with conditional duration or not
            :param ending_condition: (MetaSpeciesLogicResolver) trigger for the ending contion

            :returns: species_for_sbml (dict) = dictionary where the species strings (in MobsPy format) are the keys and
            the values are their counts inside the model. reactions_for_sbml (dict) = Reaction in dictionary format
            with keys 're', 'pr' and 'rate'. Values are strings parameters_for_sbml (dict) = value with parameter
            name and tuple with quantity and unit. For standard models it only holds the volume mappings_for_sbml
            (dict) = Dictionary that maps a meta-species to the set of species strings that belong to it. model_str
            (str) = str of the variables defined above in a user-readable format. events_for_sbml (dict) = dictionary
            containing the event trigger and count assignments. assigned_species (list) = list of meta-species which
            had counts assigned to them during this model execution (used when composing simulations)
        """
        # We name the species according to variables names for convenience
        cls.name_all_involved_species(names)
        names_used = set()
        for i, species in enumerate(meta_species_to_simulate):
            if '$' in species.get_name():
                simlog.error('An error has occurred and one of the species was not named')
            if species.get_name() in names_used:
                simlog.error(f'Names must be unique for all species\n' +
                             f'The repeated name is {species.get_name()} in position {i}\n' +
                             f'Another possibility could be a repeated meta-species in the model')
            names_used.add(species.get_name())

        # Construct structures necessary for reactions
        cls.ref_characteristics_to_object = mcu.create_orthogonal_vector_structure(meta_species_to_simulate)

        # Order the species references for later usage
        for species in meta_species_to_simulate:
            species.order_references()

        # Start by creating the Mappings for the SBML
        # Convert to user friendly format as well
        mappings_for_sbml = {}
        for spe_object in meta_species_to_simulate:
            mappings_for_sbml[spe_object.get_name()] = []

        # List of Species objects
        species_for_sbml = {}
        mappings_for_sbml = {}
        for spe_object in meta_species_to_simulate:
            species_string_list = ssg.construct_all_combinations(spe_object, 'std$',
                                                                 cls.ref_characteristics_to_object)
            for x in species_string_list:
                x[0] = x[0].get_name()
            mappings_for_sbml[spe_object.get_name()] = ['.'.join(x) for x in species_string_list]
            for species_string in species_string_list:
                species_for_sbml['_dot_'.join(species_string)] = 0

        # Dimension check. Here we check based on the units the area dimension
        dimension = None

        for reaction in reactions_set:
            if isinstance(reaction.rate, Quantity):
                if uh.extract_length_dimension(str(reaction.rate.dimensionality), dimension):
                    dimension = uh.extract_length_dimension(str(reaction.rate.dimensionality), dimension)

        for count in species_counts:
            if isinstance(count['quantity'], Quantity):
                if uh.extract_length_dimension(str(count['quantity'].dimensionality), dimension):
                    dimension = uh.extract_length_dimension(str(count['quantity'].dimensionality), dimension)

        volume = uh.convert_volume(volume, dimension)

        # assign them default order
        for reaction in reactions_set:
            if reaction.order is None:
                reaction.order = default_order

        # Add the flag species used for verifying if the simulation is over
        if continuous_sim:
            species_for_sbml[EndFlagSpecies.get_name()] = 0

        # Set initial counts for SBML
        # Create the list HERE
        assigned_species = []
        for count in species_counts:

            species_string = ssg.construct_species_char_list(count['object'], count['characteristics'],
                                                             cls.ref_characteristics_to_object,
                                                             symbol='_dot_')

            temp_count = uh.convert_counts(count['quantity'], volume, dimension)
            if type(temp_count) == float and not type_of_model == 'deterministic':
                simlog.warning('The stochastic simulation rounds floats to integers')
                species_for_sbml[species_string] = int(temp_count)
                assigned_species.append(species_string)
            else:
                species_for_sbml[species_string] = temp_count
                assigned_species.append(species_string)

        # BaseSpecies reactions for SBML with theirs respective parameters and rates
        # What do I have so far
        # Species_String_Dict and a set of reaction objects in Reactions_Set
        reactions_for_sbml, parameters_for_sbml = rc.create_all_reactions(reactions_set,
                                                                          meta_species_to_simulate,
                                                                          cls.ref_characteristics_to_object,
                                                                          type_of_model, dimension)
        parameters_for_sbml['volume'] = (volume, f'dimensionless')

        # O(n^2) reaction check for doubles
        for i, r1 in enumerate(reactions_for_sbml):
            for j, r2 in enumerate(reactions_for_sbml):
                if i == j:
                    continue

                if reactions_for_sbml[r1]['re'] == reactions_for_sbml[r2]['re'] \
                        and reactions_for_sbml[r1]['pr'] == reactions_for_sbml[r2]['pr'] \
                        and reactions_for_sbml[r1]['kin'] == reactions_for_sbml[r2]['kin']:
                    simlog.warning('The following reaction: \n' +
                                   f'{reactions_for_sbml[r1]} \n' +
                                   'Is doubled. Was that intentional? \n')

        # Event implementation here
        events_for_sbml = eh.format_event_dictionary_for_sbml(species_for_sbml, event_dictionary,
                                                              cls.ref_characteristics_to_object,
                                                              volume, dimension, meta_species_to_simulate)

        if continuous_sim:
            end_event = {'trigger': ending_condition.generate_string(cls.ref_characteristics_to_object),
                         'delay': '0',
                         'assignments': [('End_Flag_MetaSpecies', '1')]}

            reactions_for_sbml['reaction_end'] = {'re': [(10, 'End_Flag_MetaSpecies')], 'pr': [],
                                                  'kin': 'End_Flag_MetaSpecies * 1e-60'}
            events_for_sbml['end_event'] = end_event

        model_str = ''
        if verbose:
            model_str = '\n'
            model_str += 'Species' + '\n'
            species_alpha = list(sorted(species_for_sbml.keys()))
            for spe in species_alpha:
                model_str += spe.replace('_dot_', '.') + ',' + str(species_for_sbml[spe]) + '\n'

            model_str += '\n'
            model_str += 'Mappings' + '\n'
            mappings_alpha = list(sorted(mappings_for_sbml.keys()))
            for map in mappings_alpha:
                model_str += map + ' :' + '\n'
                for element in sorted(mappings_for_sbml[map]):
                    model_str += element + '\n'

            model_str += '\n'
            model_str += 'Parameters' + '\n'
            parameters_alpha = list(sorted(parameters_for_sbml.keys()))
            for par in parameters_alpha:
                model_str += par + ',' + str(parameters_for_sbml[par][0]) + '\n'

            model_str += '\n'
            model_str += 'Reactions' + '\n'
            remove_hack_end_reaction = deepcopy(reactions_for_sbml)
            remove_hack_end_reaction.pop('reaction_end', None)
            reaction_alpha = [str(x[1]).replace('_dot_', '.') for x in
                              list(sorted(remove_hack_end_reaction.items(), key=lambda x: str(x[1])))]

            for i, reac in enumerate(reaction_alpha):
                model_str += 'reaction_' + str(i) + ',' + reac + '\n'

            if events_for_sbml != {}:
                model_str += '\n'
                model_str += 'Events' + '\n'
                list_to_sort = [str(events_for_sbml[key]) for key in events_for_sbml]
                list_to_sort = sorted(list_to_sort)
                for i in range(len(list_to_sort)):
                    model_str += ('event_' + str(i) + ',' + list_to_sort[i] + '\n').replace('_dot_', '.')

        return species_for_sbml, reactions_for_sbml, parameters_for_sbml, mappings_for_sbml, model_str, \
               events_for_sbml, assigned_species


class Reactions:
    """
        This is the Reaction class. It contains the reactants, products, rate and order
        Reactions are created with the >> operator. In each reaction is stored in all involved objects

        :param reactants: (list of Species) list of meta-species used as reactants in the meta-reaction
        :param products: (list of Species) list of meta-species used as products in the meta-reaction
        :param order: (Order Operator) reaction order operator - default in a Round-Robin base
        :param rate: (int, float, callable, Quantity) reaction rate
    """

    @staticmethod
    def __create_reactants_string(list_of_reactants):
        """
            Just a simple way to simlog.debug reactions for debbuging
            Not relevant for simulation
            Important: here reactants are used interchangeably with products, this works for a list_of_products too

            :param list_of_reactants: (list of Species) list of products or reactants in the reaction
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
        """
            Prints the meta-reaction in the format A + B -> C + D
        """
        return self.__create_reactants_string(self.reactants) + \
               ' -> ' + self.__create_reactants_string(self.products)

    def __getitem__(self, item):
        """
            Override of __getitem__ for dealing with reaction rates

            :param item: (int, float, callable, Quantity) = reaction rate
        """
        return Compiler.override_get_item(self, item)

    def __init__(self, reactants, products):
        """
            Constructor of the reaction object. For the object construction only the reactants and products are
            necessary - the order and the rate are assigned later by the compiler

            :param reactants: (list of Species) list of meta-species used as reactants in the meta-reaction
            :param products: (list of Species) list of meta-species used as products in the meta-reaction
        """
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
        """
            Sets the stored reaction rate. See Compiler.override_get_item(self, item)

            :param rate: (int, float, callable, Quantity) = reaction rate
        """
        self.rate = rate


# HERE FABRICIO
class Reacting_Species(ReactingSpeciesComparator):
    """
        This is a intermediary object created when a species is used in a reaction. It is created when a species is
        part of a reaction, so it's constructor is called by the Species class in several of it's methods. It
        transforms a Species object into a list packable object to pass to create the reaction
        object. This list format allows for easy concatenation of reactants by summing

        The >> operator calls the Reaction constructor to finally define a reaction object

        :param  list_of_reactants: (list of dict) used interchangeably for products and reactants.
        A list of dictionaries containing the following {'object': meta-species object reference,
        'characteristics': the characteristics queried, 'stoichiometry': stoichiometry value,
        'label': label if used (None)}
    """

    def __str__(self):
        """
            String representation of the list of reactants
        """
        if len(self.list_of_reactants) == 1:
            to_return = str(self.list_of_reactants[0]['object'])
            for cha in self.list_of_reactants[0]['characteristics']:
                to_return += '.' + cha
            return to_return
        else:
            return self.list_of_reactants

    # Labels and value function implementation
    def c(self, item):
        """
            c query implementation, queries by the value inside a variable instead of the name
            it just calls __getattr__ with the value inside the variable

            Parameters:
                item : value to query over
        """
        item = str(item)
        if item in ['c']:
            simlog.error(f'Characteristic {item} is not allowed. Please pick another name')
        return self.__getattr__(item)

    def label(self, label):
        """
            Label function implementation. This function matches equal meta-species if they have equal labels

            Parameters:
                label (int, float, str) = value for the label for matching
        """
        if len(self.list_of_reactants) == 1:
            self.list_of_reactants[0]['label'] = label
        else:
            simlog.error('Error assigning label to reactant')
        return self

    def __getitem__(self, item):
        """
            Override of __getitem__ for dealing with reaction rates

            Parameters:
                item (int, float, callable, Quantity) = reaction rate
        """
        return Compiler.override_get_item(self, item)

    def __init__(self, object_reference, characteristics, stoichiometry=1, label=None):
        """
            Reacting_Species constructor. It receives the meta-species object reference, the characteristics that
            have been used as a query in the reaction, the stoichiometry of the meta-species in the reaction, and
            finally a label if used.

            Parameters:
                object_reference (Species) = meta-species object reference
                characteristics (str) = characteristics used to query over the meta-species inside this reaction
                stoichiometry (int) = stoichiometry value of the meta-species in the reaction
                label (int, float, str) = value for the label for matching used in this reaction
        """
        super(Reacting_Species, self).__init__()
        if object_reference.get_name() == 'S0' and characteristics == set():
            self.list_of_reactants = []
        else:
            self.list_of_reactants = [{'object': object_reference, 'characteristics': characteristics,
                                       'stoichiometry': stoichiometry, 'label': label}]

    def __rmul__(self, stoichiometry):
        """
            Multiplication by the stoichiometry for reactions

            Parameters:
                stoichiometry (int) = stoichiometry value of the meta-species in the reaction
        """
        if type(stoichiometry) == int:
            self.list_of_reactants[0]['stoichiometry'] = stoichiometry
        else:
            simlog.error(f'Stoichiometry can only be an int - Received {stoichiometry}')
        return self

    def __add__(self, other):
        """
            Addition of meta-species to construct the reaction. Here we can add species ( by transforming into reacting
            species ) or other reacting species. With this we concatenate different reacting species in the
            list_of_reactants to eventually transform into a reaction with the >> operator

            Parameter:
                other (Species or Reacting Species) = other object being added
        """
        if isinstance(other, Species):
            other = Reacting_Species(other, set())
        self.list_of_reactants += other.list_of_reactants
        return self

    def __rshift__(self, other):
        """
            The >> operator for defining reactions. It passes two instances of reacting species to construct the
            list_of_reactants and list_of_products in the reaction object

            Parameter:
                other (Species or Reacting Species) = product side of the reaction being added
        """
        if isinstance(other, Species):
            p = Reacting_Species(other, set())
        else:
            p = other

        reaction = Reactions(self.list_of_reactants, p.list_of_reactants)
        return reaction

    def __call__(self, quantity):
        """
            The call operator here is used to add counts to species non-default state. This stores the characteristics
            that have been called using the dot operator to assign the count after the call operation

            Parameters:
                quantity (int, float, Quantity) = count to be assigned to the species
        """
        if type(quantity) == int or type(quantity) == float or isinstance(quantity, Quantity):
            if len(self.list_of_reactants) != 1:
                simlog.error('Assignment used incorrectly')
            species_object = self.list_of_reactants[0]['object']
            characteristics = self.list_of_reactants[0]['characteristics']
            quantity_dict = species_object.add_quantities(characteristics, quantity)
        else:
            simlog.error('Reactant_species does not support this type of call')

        dummy_object = self.list_of_reactants[0]['object']
        if dummy_object._simulation_context is not None:
            try:
                dummy_object._simulation_context.current_event_count_data.append({'species': dummy_object,
                                                                                  'characteristics': quantity_dict[
                                                                                      'characteristics'],
                                                                                  'quantity': quantity_dict[
                                                                                      'quantity']})
            except:
                simlog.error('Only species count assignments are allowed in a model context')
        else:
            return self

    def __getattr__(self, characteristic):
        """
            This is the implementation of the .dot operation. Adds characteristics to the species and/or perform
            queries

            Parameters:
                characteristic (str) = characteristic for the query
        """
        for reactant in self.list_of_reactants:

            species_object = reactant['object']
            characteristics_from_references = mcu.unite_characteristics(species_object.get_references())

            if characteristic not in characteristics_from_references:
                species_object.add_characteristic(characteristic)

            reactant['characteristics'].add(characteristic)

        return self

    @classmethod
    def is_species(cls):
        return False


class ParallelSpecies:
    """
        This class just stores species after the | operation. It creates a list of species than can be loop through or
        given to the simulator

        Attributes:
            list_of_species (Species) = Meta-species list to store the meta-species
    """

    def __init__(self, list_of_species):
        """
            Constructor not usually used - but a list_of_species is given one can construct the ParallelSpecies from it
            It is advised to use the | operator

            Parameters:
                list_of_species (Species) = Meta-species list to store the meta-species
        """
        self.list_of_species = list_of_species

    def append(self, species):
        """
            Add species to the list_of_species called by the | operator

            Parameters:
                species (Species) = meta-species to be added to the ParallelSpecies
        """
        if not isinstance(species, Species):
            simlog.error('Only Species can be appended')
        self.list_of_species.append(species)

    def __getitem__(self, item):
        """
            Equal to list[item]

            Parameters:
                item (int) = value of index to get the element
        """
        return self.list_of_species[item]

    def __str__(self):
        """
            String representation of the ParallelSpecies just the list of meta-species
        """
        to_return = []
        for spe in self.list_of_species:
            to_return.append(spe.get_name())
        return str(to_return)

    def __or__(self, other):
        """
            Implementation of the | operator. Adds species to the Parallel species

            Parameters:
                other (Species) = meta-species to be added to the ParallelSpecies
        """
        if isinstance(other, Species):
            self.append(other)
            return self
        elif isinstance(other, ParallelSpecies):
            return ParallelSpecies(self.list_of_species + other.list_of_species)
        else:
            simlog.error('Operator must only be used in Species on ParallelSpecies')

    def __iter__(self):
        """
            Iterator through the species inside the ParallelSpecies
        """
        for spe in self.list_of_species:
            yield spe


class Species(SpeciesComparator):
    """
        Fundamental class - The meta-species object
        Contains the characteristics, the species name, the reactions it is involved in
        and finally the other species it references
        Objects store all the basic information necessary to create an SBML file and construct a model
        So we construct the species through reactions and __getattr__ to form a model

        Parameters:
            _name (str) = name of the species - named firstly as N$Counter as a placeholder. During compilation it is
            named by MobsPy to be equal to the variable name. It can also be named by the user using one the methods
            _characteristics (str) = set of characteristics DIRECTLY added to a species, does not contain inherited
            characteristics
            _references (set) = set of meta-species a meta-species has inherited from (used in it's construction).
            Every species contains itself in this set
            first_characteristic (str) = first characteristic added to the species
            _reactions (set) = Every species stores all reactions is involved in. The compiler performs the union of
            the set with all species and species in the _references set to get the reactions back
            _species_counts (list) = counts listed for the species. Stores dictionaries with the characteristics as
            keys and the counts of values

        Methods:
            __str__, c, label, show_reactions, show_characteristics, show_references, show_quantities, __or__,
            __iter__, __getitem__, __rmul__, __radd__, __add__, __rshift__, __getattr__, __call__, add_quantity,
            reset_quantities, get_quantities, __mul__, __init__, get_name, add_characteristic, remove_characteristic,
            print_characteristics, get_references, add_reference, set_references, set_reactions, add_reactions

    """

    def __str__(self):
        """
            String representation, just returns the species name
        """
        return self._name

    # Def c to get the value
    def c(self, item):
        """
            c query implementation, queries by the value inside a variable instead of the name
            it just calls __getattr__ with the value inside the variable

            Parameters:
                item : value to query over
        """
        item = str(item)
        if item in ['c']:
            simlog.error(f'Characteristic {item} is not allowed. Please pick another name')
        return self.__getattr__(item)

    # Def labels
    def label(self, label):
        """
            Label function implementation. This function matches equal meta-species if they have equal labels
            Here it returns a Reacting_Species as labels need to be used in reactions

            Parameters:
                label (int, float, str) = value for the label for matching

            Returns:
                Reacting_Species object created with the label in the constructor
        """
        return Reacting_Species(self, set(), label=label)

    # Get data from species for debugging Compiler
    def show_reactions(self):
        """
            Prints the reactions inside the object
        """
        simlog.debug(str(self) + '_dot_')
        for reference in self._references:
            for reaction in reference.get_reactions():
                simlog.debug(reaction)

    def show_characteristics(self):
        """
            Prints the directly characteristics inside the object
        """
        simlog.debug(str(self) + ' has the following characteristics referenced:')
        for i, reference in enumerate(self.get_references()):
            if reference.get_characteristics():
                simlog.debug(str(reference) + ': ', end='')
                reference.simlog.debug_characteristics()

    def show_references(self):
        """
            Prints the objects this object has inherited from (including self)
        """
        simlog.debug(str(self) + '_dot_')
        simlog.debug('{', end='')
        for i, reference in enumerate(self.get_references()):
            if reference.get_characteristics():
                simlog.debug(' ' + str(reference) + ' ', end='')
        simlog.debug('}')

    def show_quantities(self):
        """
            Shows the species counts stored in this object
        """
        simlog.debug(self._species_counts)

    # Creation of ParallelSpecies For Simulation ##################
    def __or__(self, other):
        """
            Creates an instance of ParallelSpecies using the | operator

            Parameters:
                other (Species or ParallelSpecies) = Species or ParallelSpecies to combine into
        """
        if isinstance(other, ParallelSpecies):
            other.append(self)
            return other
        else:
            return ParallelSpecies([self, other])

    def __iter__(self):
        """
            iter defined to be consistent with ParallelSpecies behavior

            Returns:
                Itself
        """
        yield self

    # Creation of reactions using entities ########################
    def __getitem__(self, item):
        """
            Override of __getitem__ for dealing with reaction rates

            Parameters:
                item (int, float, callable, Quantity) = reaction rate
        """
        return Compiler.override_get_item(self, item)

    def __rmul__(self, stoichiometry):
        """
            Multiplication by the stoichiometry

            Parameters:
                stoichiometry (int) = Stoichiometry of the meta-species in the meta-reaction

            Returns:
                r (Reacting_Species) = Reacting_Species with the stoichiometry added to it
        """
        if type(stoichiometry) == int:
            r = Reacting_Species(self, set(), stoichiometry)
        else:
            simlog.error(f'Stoichiometry can only be an int - Received {stoichiometry}')
        return r

    def __add__(self, other):
        """
            Implementation of addition. In case a Species object is being used in a reaction, so it can then generate
            a Reacting_Species object

            Parameters:
                other (Species or Reacting Species object) = Other objected added to construct a reaction

            Returns:
                r1 + r2 (Reacting Species) = Reacting Species objected created by the sum of the two
        """
        r1 = Reacting_Species(self, set())
        if isinstance(other, Reacting_Species):
            r2 = other
        else:
            r2 = Reacting_Species(other, set())
        return r1 + r2

    def __radd__(self, other):
        """
            Just making addition symmetric check __add__
        """
        Species.__add__(self, other)

    def __rshift__(self, other):
        """
            Reaction definition (>> operator) in case one uses only a single species. Allows it to transform
            into a reacting species before adding it into the reaction

            Parameters:
                other (Species or Reacting_Species) = Reaction products

            Returns
                The reaction
        """
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
            This adds characteristics using the following manner Species.characteristic
            If it is the first time is called we add it to the set of characteristics
            Second time, it just creates a Reacting_Species for reaction construction

            Parameters:
                characteristic (str) = characteristic to be added or to be use as a query in the reaction

            Returns:
                Reacting_Species with the characteristic added for querying
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
        """
            The __call__ operator handles two things for Species objects
            First, it adds characteristics to the default state of the meta-species if the quantity is a real
            Secondly, it returns the name of the characteristic from a species string that belongs to this meta-species

            Parameters:
                quantity (int, float, Quantity) = for count assignment
                        (Specific_Species_Operator object) = for characteristic extraction

            Returns:
                self = to allow for assigning counts mid-reaction
        """
        if type(quantity) == int or type(quantity) == float or isinstance(quantity, Quantity):
            quantity_dict = self.add_quantities('std$', quantity)
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

        if self._simulation_context is not None:
            try:
                self._simulation_context.current_event_count_data.append({'species': self,
                                                                          'characteristics': quantity_dict[
                                                                              'characteristics'],
                                                                          'quantity': quantity_dict['quantity']})
            except:
                simlog.error('Only species count assignments are allowed in a model context')
        else:
            return self

    def add_quantities(self, characteristics, quantity):
        '''
            This function was implemented because Python can't accept sets as dictionary keys
            Set the quantity of an specific string of species

            Parameters
                characteristics (str) = characteristics of the species to be set
                quantity (int, float, Quantity) = counts of that specific species
        '''
        if self._simulation_context is None:
            already_in = False
            for e in self._species_counts:
                if characteristics == e['characteristics']:
                    e['quantity'] = quantity
                    already_in = True
            if not already_in:
                self._species_counts.append({'characteristics': characteristics, 'quantity': quantity})
        else:
            return {'characteristics': characteristics, 'quantity': quantity}

    def reset_quantities(self):
        """
            Just reset the counts inside a species - potentially useful in jupyter environment
        """
        self._species_counts = []

    def get_quantities(self):
        """
            Returns:
                The list of species_counts
        """
        return self._species_counts

    # Creation of other objects through multiplication ##################
    def __mul__(self, other):
        """
            Multiplications are used to construct more complex species
            We do not concatenate characteristics. This keeps the structure intact for reaction construction
            Instead we combine the sets of references, objects can access the characteristics from other through the
            reference set

            Parameters:
                other (Species) = Species for the multiplication and creation of a higher order species

            Returns:
                new_entity (Species) = Higher order species resulted from the multiplication
        """
        Compiler.entity_counter += 1
        name = 'N$' + str(Compiler.entity_counter)
        new_entity = Species(name)
        new_entity.set_references(mcu.combine_references(self, other))
        new_entity.add_reference(new_entity)

        mcu.check_orthogonality_between_references(new_entity.get_references())

        return new_entity

    def __init__(self, name):
        """
            Object constructor - We recommend using BaseSpecies instead

            Parameters:
                name (str) = Name of the species (can be placeholder if named with N$)
        """
        super(Species, self).__init__()
        self._name = name
        self._characteristics = set()
        self._references = {self}
        self._ordered_references = []
        self._simulation_context = None
        self._reference_index_dictionary = {}

        # This is necessary for the empty objects generated when we perform multiplication with more than 2 Properties
        self.first_characteristic = None

        # Each object stores the reactions it is involved in
        self._reactions = set()

        # This will store the quantities relating to the species counts
        self._species_counts = []

    def name(self, name):
        """
            Function to name a species. Allows the user to use a name different from the variable's name

            Parameters:
                name (str) = name of the species
        """
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

    def reset_references(self):
        self._references = {self}

    def get_reactions(self):
        return self._reactions

    def set_reactions(self, reactions):
        self._reactions = reactions

    def reset_reactions(self):
        self._reactions = set()

    def add_reaction(self, reaction):
        self._reactions.add(reaction)

    def reset_counts(self):
        self._species_counts = []

    def set_simulation_context(self, sim):
        if self._simulation_context is None:
            self._simulation_context = sim
        else:
            simlog.error('A different Simulation Object was assigned to a meta-species object under context \n'
                         'Please use only one Simulation Object per context assignment')

    def reset_simulation_context(self):
        self._simulation_context = None

    def order_references(self):
        cleaned_references = [x for x in self.get_references() if x.get_characteristics() != set()]
        self._ordered_references = sorted(cleaned_references, key=lambda x: sorted(list(x.get_characteristics())))
        i = 1
        for reference in self._ordered_references:
            self._reference_index_dictionary[reference] = i
            i = i + 1

    def get_ordered_references(self):
        return self._ordered_references

    def get_index_from_reference_dict(self, reference):
        return self._reference_index_dictionary[reference]

    @classmethod
    def is_species(cls):
        return True


def compile_species_number_line(function, code_line):
    pos = code_line.find('=')
    if pos == -1:
        return 0
    else:
        code_line = code_line[0:pos]
    n = code_line.count(',') + 1

    if n > 1:
        return [function(1) for _ in range(n)]
    else:
        return function(1)


# Property Call to return several properties as called
def BaseSpecies(number_of_properties=None):
    """
        Function that returns a number of base species according to the number_of_properties
        BaseSpecies are just Species with no inheritance attached to it

        Parameters:
            number_of_properties (int) = number of base species wanted

        Returns:
            Species objects according to the number of species requested
    """
    if number_of_properties is None:
        code_line = inspect.stack()[1].code_context[0][:-1]
        return compile_species_number_line(BaseSpecies, code_line)

    to_return = []
    for i in range(number_of_properties):
        Compiler.entity_counter += 1
        name = 'N$' + str(Compiler.entity_counter)
        to_return.append(Species(name))

    if len(to_return) == 1:
        return to_return[0]
    else:
        return tuple(to_return)


__S0, __S1, __SF = BaseSpecies(3)
__S0.name('S0')
__S1.name('S1')
__SF.name('End_Flag_MetaSpecies')
EndFlagSpecies = __SF
Zero = __S0


def New(species, n=None):
    """
        New is just a multiplication in disguise. It contains a __S1 meta-species that multiplies the species used as a
        parameter. The species __S1 contains no characteristics and it's only used to create a higher-dimensional
        species from the species used as an argument
        It can also return several meta-species objects at a time and be used for designing list of species ( check
        tutorial)

        Parameters:
            n (int) = number of species to return
    """
    if n is None:
        code_line = inspect.stack()[1].code_context[0][:-1]
        return compile_species_number_line(lambda i: New(species, i), code_line)

    if type(n) == int:
        to_return = [species * __S1 for _ in range(n)]
        if n == 1:
            return to_return[0]
        else:
            return tuple(to_return)
    elif type(n) == str:
        to_return = species * __S1
        to_return.name(n)
        return to_return
    else:
        simlog.error('New only supports numbers for creating multiple species \n '
                     ' or exclusively string for single species name')


if __name__ == '__main__':
    pass
