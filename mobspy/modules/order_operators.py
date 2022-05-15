"""
    This module deals with order_operators
    Thus it deals with the Round-Robin order assignment and with species in the products that
    are not referenced in the reactants (Born Species)
    It has also the implementation of the Rev operator, since that is a also a reaction operator
"""
from copy import deepcopy
import mobspy.modules.meta_class_utils as mcu
import mobspy.modules.meta_class as mc
import mobspy.simulation_logging.log_scripts as simlog


class __Operator_Base:
    """
        This is the order operator base. It contains some functions that can be useful the reaction operators defined
        Reaction Operators are now responsible for order assignment and how to deal with species in the products that
        are not referenced in the reactants (Born Species)
        It also overrides the getitem method so we can use:
        Operator[Reaction] to assign an order to the reaction
    """

    # Assign order structure
    def __getitem__(self, item):
        item.order = self

    # Transform product function
    @staticmethod
    def transform_species_string(species_string, characteristics_to_transform,
                                 ref_characteristics_to_object):
        """
            This function handles the characteristic change due to a chemical reaction. It receives a string in the
            MobsPy format and the characteristic one has written on the product side of the reaction
            MobsPy them converts the characteristics in the same position as that one (that have been directly
            added to the same meta-species) to the one in the characteristics_to_transform variable

            Parameters:
                species_string (str) = A string associated with a species in the MobsPy format
                characteristics_to_transform (str) = Characteristics received by the species object in the product
                ref_characteristics_to_object (dict) = Dictionary with the characteristics as keys and meta-species
                objects as values
        """
        species_to_return = deepcopy(species_string)
        for characteristic in characteristics_to_transform:
            replaceable_characteristics = ref_characteristics_to_object[characteristic].get_characteristics()

            for rep_cha in replaceable_characteristics:
                species_to_return = species_to_return.replace('_dot_' + rep_cha, '_dot_' + characteristic)

        return species_to_return

    @staticmethod
    def find_all_string_references_to_born_species(species_referenced_by, characteristics, species_string_dictionary):
        """
            Finds all the string species that reference the born-meta species, including inheritors
            Returns them all in a list for the product construction

            Parameters:
                species_referenced_by (meta-species list) = List of species that inherit from the born-meta species
                in the reaction
                characteristics (str) = List of characteristics queried on the born species in the product
                species_string_dictionary (dict) = Dictionary with meta-species as keys and a list of species strings
                as values
        """
        to_return = []
        for species in species_referenced_by:
            for key in species_string_dictionary:
                if species in key.get_references():
                    for species_string in species_string_dictionary[key]:
                        species_string_split = species_string.split('_dot_')[1:]
                        if all([char in species_string_split for char in characteristics]):
                            to_return.append(species_string)

        return to_return

    @staticmethod
    def find_all_default_references_to_born_species(species_referenced_by, characteristics, species_string_dictionary,
                                                    ref_characteristics_to_object):
        """
            Finds only the DEFAULT string species that reference the born-meta species, including inheritors
            Returns them all in a list for the product construction.

            Parameters:
                species_referenced_by (meta-species list) = List of species that inherit from the born-meta species
                in the reaction
                characteristics (str) = List of characteristics queried on the born species in the product
                species_string_dictionary (dict) = Dictionary with meta-species as keys and a list of species strings
                as values
                ref_characteristics_to_object (dict) =  Dictionary with the characteristics as keys and meta-species
                objects that they have been directly added to as values
        """
        to_return = []

        for species in species_referenced_by:
            characteristics_to_find = mcu.complete_characteristics_with_first_values(species, characteristics,
                                                                                     ref_characteristics_to_object)

            characteristics_to_find.remove(species.get_name())
            characteristics_to_find.union(characteristics)

            for key in species_string_dictionary:
                if species in key.get_references():
                    for species_string in species_string_dictionary[key]:
                        species_string_split = species_string.split('_dot_')[1:]
                        if all([char in species_string_split for char in characteristics_to_find]):
                            to_return.append(species_string)

        return to_return


class __Round_Robin_Base(__Operator_Base):
    """
        Here we have the implementation of the round robin order it goes like this:

        2*Ecoli >> 4*Ecoli

        Since an Ecoli is a generic object that refers to all types of Ecoli
        This object refers to multiple reactions and to multiple strings
        Therefore we follow this rule to know which string goes where:
        1 Ecoli reactant >> 1 Ecoli product
        2 Ecoli reactant >> 2 Ecoli product
        1 Ecoli reactant >> 3 Ecoli product
        2 Ecoli reactant >> 4 Ecoli product
        This is a cycle (round robin) between the reactants to be assigned positions in the product
        The products will keep the characteristics of the reactants except if stated otherwise with .
        For completely new species (no reactant of the same species) we use ALL possible combinations
        For only the default option see the code bellow

        Methods:
            __call__
    """

    def __call__(self, order_dictionary, product_species,
                 species_string_dictionary,
                 ref_characteristics_to_object):
        """
            This function is responsible for parrying the products with the reactants according to their meta-species
            It parries them using the order_dictionary and an index-system for the round-robin application
            If it cannot find a parrying it considers it a born species

            Parameters:
                order_dictionary (dict) = Dictionary with the meta-species as keys and the list of meta-species strings
                product_species (dict) = Product species dict with the meta-species object (key: species),
                label (key: label), and characteristics (key: characteristics)
                ref_characteristics_to_object (dict) =  Dictionary with the characteristics as keys and meta-species
                objects that they have been directly added to as values
        """
        round_robin_index = {}
        for species, label in [(e['species'], e['label']) for e in product_species]:
            round_robin_index[(species, label)] = 0

        products = []
        for species, label, characteristics in [(e['species'], e['label'], e['characteristics']) for e in
                                                product_species]:

            # Simple round robin
            try:
                species_to_transform_string = order_dictionary[(species, label)][round_robin_index[(species, label)]]
                round_robin_index[(species, label)] = (round_robin_index[(species, label)] + 1) % len(
                    order_dictionary[(species, label)])

                # Return in list of lists format for combination later
                products.append([self.transform_species_string(species_to_transform_string, characteristics,
                                                               ref_characteristics_to_object)])

            # If the species is not on the reactants - order_dictionary
            except KeyError:

                # Find all the species that reference the one in the reaction
                species_is_referenced_by = []
                for key in species_string_dictionary:
                    if species in key.get_references():
                        species_is_referenced_by.append(key)

                if len(species_is_referenced_by) == 0:
                    simlog.error(f'{species} or inheritors are not in the model. Please add at least one')

                products.append(self.find_all_string_references_to_born_species(species_is_referenced_by,
                                                                                characteristics,
                                                                                species_string_dictionary))

        return products


# Define class to override the operators
All = __Round_Robin_Base()


class __RR_Default_Base(__Operator_Base):
    """
        Only default options for born species (no reference in reactant)
        See __Round_Robin_Base for clarification

        Methods:
            __call__
    """

    # Here is the default order requested by Thomas
    def __call__(self, order_dictionary, product_species,
                 species_string_dictionary,
                 ref_characteristics_to_object):
        """
            This function is responsible for parrying the products with the reactants according to their meta-species
            It parries them using the order_dictionary and an index-system for the round-robin application
            If it cannot find a parrying it considers it a born species

            Parameters:
                order_dictionary (dict) = Dictionary with the meta-species as keys and the list of meta-species strings
                product_species (dict) = Product species dict with the meta-species object (key: species),
                label (key: label), and characteristics (key: characteristics)
                ref_characteristics_to_object (dict) =  Dictionary with the characteristics as keys and meta-species
                objects that they have been directly added to as values
        """
        round_robin_index = {}
        for species, label in [(e['species'], e['label']) for e in product_species]:
            round_robin_index[(species, label)] = 0

        products = []
        for species, label, characteristics in [(e['species'], e['label'], e['characteristics']) for e in
                                                product_species]:

            # Simple round robin
            try:
                species_to_transform_string = order_dictionary[(species, label)][round_robin_index[(species, label)]]
                round_robin_index[(species, label)] = (round_robin_index[(species, label)] + 1) % len(
                    order_dictionary[(species, label)])

                # Return in list of lists format for combination later
                products.append([self.transform_species_string(species_to_transform_string, characteristics,
                                                               ref_characteristics_to_object)])

            # If the species is not on the reactants - order_dictionary
            except KeyError:

                # Find all the species that reference the one in the reaction
                species_is_referenced_by = []
                for key in species_string_dictionary:
                    if species in key.get_references():
                        species_is_referenced_by.append(key)

                if len(species_is_referenced_by) == 0:
                    simlog.error(f'{species} or inheritors are not in the model. Please add at least one')

                products.append(self.find_all_default_references_to_born_species(species_is_referenced_by,
                                                                                 characteristics,
                                                                                 species_string_dictionary,
                                                                                 ref_characteristics_to_object))

        return products


Default = __RR_Default_Base()


class __Set_Reversible_Rate:
    """
        This class is responsible for dealing with reversible reactions
        Since reversible reactions have too rates
        The process is in __getitem__
    """
    def __getitem__(self, both_rates):
        """
            This function extracts both reaction rates from a tuple

            Parameters:
                both_rates (rate functions) = rate function from the direct and reverse reaction
        """
        if len(both_rates) != 2:
            simlog.error('The reversible reaction must receive 2 rates')

        self.reaction_direct.rate = both_rates[0]
        self.reaction_reverse.rate = both_rates[1]

    def __init__(self):
        """
            This a dummy constructor, the class object Rev is the one with the __getitem__ overridden
        """
        self.reaction_direct = None
        self.reaction_reverse = None

    def set_reactions(self, reaction_direct, reaction_reverse):
        self.reaction_direct = reaction_direct
        self.reaction_reverse = reaction_reverse


# Reversible reaction operator
class __Reversible_Base:

    def __getitem__(self, reaction):
        """
            This __getitem__ uses the rate setter which is a instance of __Set_Reversible_Rate to set the reversible
            reaction rates. It also creates the reverse reaction with the constructor

            Parameters:
                reaction (Reaction object) = object from the reaction class
        """
        reaction_direct = reaction
        reaction_reverse = mc.Reactions(reaction_direct.products, reaction_direct.reactants)
        self.rate_setter.set_reactions(reaction_direct, reaction_reverse)
        return self.rate_setter

    def __init__(self, rate_setter):
        """
            This a dummy constructor, the class object Rev is the one with the __getitem__ overridden
        """
        self.rate_setter = rate_setter


Rev = __Reversible_Base(__Set_Reversible_Rate())
