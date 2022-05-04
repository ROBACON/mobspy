from copy import deepcopy
import mobspy.modules.meta_class_utils as mcu
import mobspy.modules.meta_class as mc
import mobspy.simulation_logging.log_scripts as simlog


class __Operator_Base:
    """
        This is the order operator base. It contains some functions that can be useful for all operators
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
            species_string : A string associated with a species in the form of Species.c1.c2.c3
            characteristics_to_transform : Characteristics received by the species object in the product
            refs : orthogonal structure - get the object referring to the characteristic and the characteristics
            associated with the object

            For the products we take the referenced characteristic and every string that does not have that characteristic
            will have it transformed to become a reaction

            We only transform the characteristics from the same base_object, respecting the structure
            described on the paper
        """
        species_to_return = deepcopy(species_string)
        for characteristic in characteristics_to_transform:
            replaceable_characteristics = ref_characteristics_to_object[characteristic].get_characteristics()

            for rep_cha in replaceable_characteristics:
                species_to_return = species_to_return.replace('_dot_' + rep_cha, '_dot_' + characteristic)

        return species_to_return

    @staticmethod
    def find_all_string_references_to_born_species(species, characteristics, species_string_dictionary):
        to_return = []
        for key in species_string_dictionary:
            if species in key.get_references():
                for species_string in species_string_dictionary[key]:
                    species_string_split = species_string.split('_dot_')[1:]
                    if all([char in species_string_split for char in characteristics]):
                        to_return.append(species_string)

        return to_return

    @staticmethod
    def find_all_default_references_to_born_species(species, characteristics, species_string_dictionary,
                                                    ref_characteristics_to_object):
        to_return = []

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
    """

    def __call__(self, order_dictionary, product_species,
                 species_string_dictionary,
                 ref_characteristics_to_object):

        round_robin_index = {}
        for species, label in [(e['species'], e['label']) for e in product_species]:
            round_robin_index[(species, label)] = 0

        products = []
        for species, label, characteristics in [(e['species'], e['label'], e['characteristics']) for e in product_species]:

            # species.in_model = False
            # print(species, label, characteristics)

            # Simple round robin
            try:
                species_to_transform_string = order_dictionary[(species, label)][round_robin_index[(species, label)]]
                round_robin_index[(species, label)] = (round_robin_index[(species, label)] + 1) % len(order_dictionary[(species, label)])

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

                for referenciator_species in species_is_referenced_by:
                    products.append(self.find_all_string_references_to_born_species(referenciator_species,
                                                                                    characteristics,
                                                                                    species_string_dictionary))

        return products


# Define class to override the operators
All = __Round_Robin_Base()


class __RR_Default_Base(__Operator_Base):
    """
        Only default options for born species (no reference in reactant)
        See comment above for clarification
    """

    # Here is the default order requested by Thomas
    def __call__(self, order_dictionary, product_species,
                 species_string_dictionary,
                 ref_characteristics_to_object):

        round_robin_index = {}
        for species, label in [(e['species'], e['label']) for e in product_species]:
            round_robin_index[(species, label)] = 0

        products = []
        for species, label, characteristics in [(e['species'], e['label'], e['characteristics']) for e in product_species]:

            # Simple round robin
            try:
                species_to_transform_string = order_dictionary[(species, label)][round_robin_index[(species, label)]]
                round_robin_index[(species, label)] = (round_robin_index[(species, label)] + 1) % len(order_dictionary[(species, label)])

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

                for referenciator_species in species_is_referenced_by:
                    products.append(self.find_all_default_references_to_born_species(referenciator_species,
                                                                                     characteristics,
                                                                                     species_string_dictionary,
                                                                                     ref_characteristics_to_object))

        return products


Default = __RR_Default_Base()


class __Set_Reversible_Rate:
    def __getitem__(self, both_rates):
        if len(both_rates) != 2:
            simlog.error('The reversible reaction must receive 2 rates')

        self.reaction_direct.rate = both_rates[0]
        self.reaction_reverse.rate = both_rates[1]

    def __init__(self):
        self.reaction_direct = None
        self.reaction_reverse = None

    def set_reactions(self, reaction_direct, reaction_reverse):
        self.reaction_direct = reaction_direct
        self.reaction_reverse = reaction_reverse


# Reversible reaction operator
class __Reversible_Base:

    def __getitem__(self, reaction):
        # The reactions are added by the constructor
        reaction_direct = reaction
        reaction_reverse = mc.Reactions(reaction_direct.products, reaction_direct.reactants)
        self.rate_setter.set_reactions(reaction_direct, reaction_reverse)
        return self.rate_setter

    def __init__(self, rate_setter):
        self.rate_setter = rate_setter


Rev = __Reversible_Base(__Set_Reversible_Rate())
