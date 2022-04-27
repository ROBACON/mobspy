import mobspy.modules.meta_class as mc
import mobspy.modules.reaction_construction_nb as rc
import mobspy.simulation_logging.log_scripts as simlog
import mobspy.modules.meta_class as mc
import mobspy.modules.unit_handler as uh
from pint import Quantity


class Bool_Override:
    """
        This class override the boolean operator of all function rate objects
    """

    def __bool__(self):
        # Return false if the species does not exist in reaction
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
        Here we create a new species object. It is called specific because it refers to specific
        species strings

        The previous object was generic - Ecoli would refer to ALL Ecoli
        Here we deal with specific species like Ecoli.blue.happy.alive

        We override the operators to allow us to write is Ecoli.alive on function rate definition

        Operators use the constructor and what is inside it to define themselves
        After that the boolean method is called on the constructed object
        Thus allowing to set rates based on the specific characteristics of the string created in the beginning
    """

    def __init__(self, species_string, species_object):
        self.species_string = species_string
        self._stocked_characteristics = set()
        self._species_object = species_object

    def __getattr__(self, characteristic):
        self._stocked_characteristics.add(characteristic)
        return self

    def __str__(self):
        return self.species_string

    def is_a(self, reference):
        if not self._stocked_characteristics:
            if reference in self._species_object.get_references():
                return True
            else:
                return False
        else:
            simlog.error('Concatenation of is_a and dot operator still not supported. Please use them separately')

    def add(self, characteristic):
        self._stocked_characteristics.add(characteristic)


def extract_reaction_rate(combination_of_reactant_species, reactant_string_list
                          , reaction_rate_function, type_of_model, volume, dimension):
    '''
        The order of the reactants appears in species_string_dict appears equality
        to the order they appear on the reaction

        If an int or a flot is returned we apply basic kinetics from CRN
        If a str is returned we use that as it is.
        Remember to set parameters_for_sbml if needed

        remember 2-Ecoli means Ecoli + Ecoli so stoichiometry must be taken into consideration

        species_string_list is a dictionary with {'reactants' :[ list_of_reactants ], 'products':[ list_of_products ]}
        reaction_rate_function is the rate function passed by the user
        parameters_for_sbml are the parameters used for the SBML construction

        returns: the reaction kinetics as a string for SBML
    '''
    extra_species = []
    if type(reaction_rate_function) == int or type(reaction_rate_function) == float or isinstance(reaction_rate_function, Quantity):
        reaction_rate_function = uh.convert_rate(reaction_rate_function, len(reactant_string_list), dimension)
        reaction_rate_string = basic_kinetics_string(reactant_string_list,
                                                     reaction_rate_function, type_of_model, volume)

    elif callable(reaction_rate_function):
        arguments = prepare_arguments_for_callable(combination_of_reactant_species,
                                                   reactant_string_list, reaction_rate_function.__code__.co_varnames)
        rate = reaction_rate_function(**arguments)
        rate = uh.convert_rate(rate, len(reactant_string_list), dimension)

        if type(rate) == int or type(rate) == float:
            reaction_rate_string = basic_kinetics_string(reactant_string_list,
                                                         rate, type_of_model,
                                                         volume)
        elif type(rate) == str:
            reaction_rate_string = rate
            extra_species = mc.Compiler.get_extra_species_list()
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

    return reaction_rate_string, extra_species


def basic_kinetics_string(reactants, reaction_rate, type_of_model, volume):
    """
        Just assign basic kinetics string based on the received reactans and rate
        parameters_for_sbml is for the construction of the model later
        Type of model is stochastic or deterministic - Used to determine which expression to use
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
        In the form S * (S - 1)/2 * .....
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
        Just the reactants strings multiplied
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
        combination_of_reactant_species : reactant species OBJECTS involved in the reaction
        reactant_string_list : reactant species STRINGS involved in the reaction
        rate_function_arguments : arguments received by the rate function

        This function extracts the requested arguments by the rate function for a given reaction
        So then it can be given to it using **kwargs
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
