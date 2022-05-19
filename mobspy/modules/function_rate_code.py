"""
    mobspy.modules.process_result_data.py

    Handles rate_assignment to the reactions:
        The passing of arguments to be executed by the user supplied rate function
        The construction of mass action kinetics
        Handles the different types of rates supplied by the user
"""

import mobspy.modules.meta_class as mc
import mobspy.modules.reaction_construction_nb as rc
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.modules.meta_class as mc
import mobspy.modules.unit_handler as uh
from pint import Quantity


class Bool_Override:
    """
        Just a base class for implementing the . operation in the rate function arguments
        through boolean overriding

        Attributes:
            _stocked_characteristics (str) = stocks the characteristics of the queries performed by the user
            species_string (str) = string value from an individual species in MobsPY format

        Methods:
            __bool__

        related classes:
            Specific_Species_Operator - inherits from it for the __bool__ method
    """

    def __bool__(self):
        """
            The implementation of the .dot operation for rate function arguments
            Returns true when the argument possesses the characteristics
            Bool is called after the .dot operations

            Parameters:
                self

            Returns:
                True if the object contains all the characteristics queried
                False otherwise
        """
        if self.species_string == '$Null':
            return False

        species_string_split = self.species_string.split('_dot_')[1:]
        if all([char in species_string_split for char in self._stocked_characteristics]):
            to_return_boolean = True
        else:
            to_return_boolean = False

        self._stocked_characteristics = set()
        return to_return_boolean


class Specific_Species_Operator(Bool_Override):
    """
        This class creates objects from the species strings from the meta-species to pass them to the rate functions
        as arguments. It uses the Bool_Override class to return true or false to the .dot operation inside the rate
        functions

        Attributes:
            species_string (str) = A string from MobsPy meta-species format
            species_object (Species) = Meta-species set for which the species_string is contained in
            _stocked_characteristics (str) = stocks the characteristics of the queries performed by the user

        Methods:
            __getattr__, __str__, is_a, add

    """

    def __init__(self, species_string, species_object):
        """
            Constructs the object from the species strings from the meta-species to pass them to the rate functions
            as arguments.

            Parameters:
                species_string (str) = A string from MobsPy meta-species format
                species_object (Species) = Meta-species set for which the species_string is contained in

        """
        self.species_string = species_string
        self._stocked_characteristics = set()
        self._species_object = species_object

    def __getattr__(self, characteristic):
        """
            Stores the characteristics for the boolean query inside the rate function by adding them to the set

            Parameters:
                characteristic (str) = characteristic being queried
        """
        self._stocked_characteristics.add(characteristic)
        return self

    def __str__(self):
        """
            Returns the species_string from the MobsPy meta-species used in the object construction
        """
        return self.species_string

    def is_a(self, reference):
        """
            This function checks to see if the meta-species the species_string belong to has inherited from the
            parameter reference (reminder: every meta-species inherits from itself)

            Parameters:
                reference (Species) = Meta-species object

            Returns:
                True if the meta-species in Specific_Species_Operator has inherited from the reference
                False otherwise
        """
        if not self._stocked_characteristics:
            if reference in self._species_object.get_references():
                return True
            else:
                return False
        else:
            simlog.error('Concatenation of is_a and dot operator still not supported. Please use them separately')

    def add(self, characteristic):
        """
            Adds characteristics to the set of characteristics

            Parameters:
                characteristic (str) = characteristic to add to the set
        """
        self._stocked_characteristics.add(characteristic)


def extract_reaction_rate(combination_of_reactant_species, reactant_string_list
                          , reaction_rate_function, type_of_model, dimension):
    """
        This function is responsible for returning the reaction rate string for the model construction. To do this it
        does a different action depending on the type of the reaction_rate_function (we consider constants as functions)
        It passes the rate as a string expression in MobsPy standard units

        Parameters:
            combination_of_reactant_species (list of Species) = Meta-species currently being used in this reaction
            reactant_string_list (list of strings) = list of species strings in order they apear in the reaction
            reaction_rate_function (float, callable, Quantity) = rate stored in the reaction object
            dimension (int) = system dimension (for the rate conversion)

        Returns:
            reaction_rate_string (str) = the reaction kinetics as a string for SBML
            extra_species (list) = species that appear in the rate but are not reactants (COPASI will throw an error)
    """
    extra_species = []
    if type(reaction_rate_function) == int or type(reaction_rate_function) == float or isinstance(reaction_rate_function, Quantity):
        reaction_rate_function = uh.convert_rate(reaction_rate_function, len(reactant_string_list), dimension)
        reaction_rate_string = basic_kinetics_string(reactant_string_list,
                                                     reaction_rate_function, type_of_model)

    elif callable(reaction_rate_function):
        arguments = prepare_arguments_for_callable(combination_of_reactant_species,
                                                   reactant_string_list, reaction_rate_function.__code__.co_varnames)
        rate = reaction_rate_function(**arguments)
        rate = uh.convert_rate(rate, len(reactant_string_list), dimension)

        if type(rate) == int or type(rate) == float:
            reaction_rate_string = basic_kinetics_string(reactant_string_list,
                                                         rate, type_of_model)
        elif type(rate) == str:
            reaction_rate_string = rate
        elif rate is None:
            simlog.error('There is a reaction rate missing for the following reactants: \n'
                         + str(reactant_string_list))
        else:
            simlog.error('The function return a non-valid value')
    elif reaction_rate_function is None:
        simlog.error('There is a reaction rate missing for the following reactants: \n'
                     + str(reactant_string_list))
    else:
        simlog.debug(type(reaction_rate_function))
        simlog.error('The rate type is not supported')

    return reaction_rate_string


def basic_kinetics_string(reactants, reaction_rate, type_of_model):
    """
        This constructs the bases for mass-action kinetics. Both for stochastic and deterministic depending on the type
        of model

        Parameters:
            reactants (list of str) = list of reactants in MobsPy str format
            reaction_rate (float) = reaction constant
            type_of_model (str) = stochastic or deterministic - rate expressions are different depending on each case

        Returns:
            kinetics_string (str) = mass action kinetics expression for the reaction
    """
    counts = rc.count_string_dictionary(reactants)

    if type_of_model.lower() not in ['stochastic', 'deterministic']:
        simlog.error('Type of model not supported. Only stochastic or deterministic')

    kinetics_string = ""
    for name, number in counts.items():
        if type_of_model.lower() == 'stochastic':
            kinetics_string += stochastic_string(name, number)
        elif type_of_model.lower() == 'deterministic':
            kinetics_string += deterministic_string(name, number)
        kinetics_string += ' * '

    kinetics_string += str(reaction_rate)

    n = kinetics_string.count('*') - 1
    if n > 0:
        kinetics_string += f' * volume^{-n}'

    return kinetics_string


def stochastic_string(reactant_name, number):
    """
        This function returns the stochastic string expression for mass action kinetics
        For instance the reaction 2A -> 3A would imply A*(A-1)/2
        It only does so for one reactant, so it must be called for all reactants in the reaction

        Parameters:
            reactant_name (str) = species string involved in the reaction
            number (int) = stoichiometry (number of times it appears)

        Returns:
            to_return_string (str) = the mass action kinetics string expression for only that species
    """
    to_return_string = ''
    for i in range(number):
        if i == 0:
            to_return_string = reactant_name
        else:
            to_return_string += f' * ({reactant_name} - {i})/{i + 1}'

    return to_return_string


def deterministic_string(reactant_name, number):
    """
        This function returns the deterministic string expression for mass action kinetics
        For instance the reaction 2A -> 3A would imply A*A
        It only does so for one reactant, so it must be called for all reactants in the reaction

        Parameters:
            reactant_name (str) = species string involved in the reaction
            number (int) = stoichiometry (number of times it appears)

        Returns:
            to_return_string (str) = the mass action kinetics string expression for only that species
    """
    to_return_string = ''
    for i in range(number):
        if i == 0:
            to_return_string = reactant_name
        else:
            to_return_string += f' * {reactant_name}'
    return to_return_string


def prepare_arguments_for_callable(combination_of_reactant_species, reactant_string_list, rate_function_arguments):
    """
        This function prepares the requested arguments to the rate function for a given reaction by creating objects
        of the Specific_Species_Operator class

        Parameters:
            combination_of_reactant_species : meta-species involved in the reaction
            reactant_string_list : species strings involved in the reaction
            rate_function_arguments : arguments received by the rate function
    """
    argument_dict = {}
    i = 0
    for i, (species, reactant_string) in enumerate(zip(combination_of_reactant_species, reactant_string_list)):
        try:
            species = species['object']
            argument_dict[rate_function_arguments[i]] = Specific_Species_Operator(species_string=reactant_string,
                                                                                  species_object=species)
        except IndexError:
            continue

    while len(argument_dict) < len(rate_function_arguments):
        argument_dict[rate_function_arguments[i]] = Specific_Species_Operator('$Null', mc.Zero)
        i += 1

    return argument_dict


if __name__ == '__main__':
    pass
