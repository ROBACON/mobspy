import modules.meta_class as mc
import modules.reaction_construction_nb as rc
import simulation_logging.log_scripts as simlog


class Bool_Override:
    """
        This class override the boolean operator of all function rate objects
    """

    def __bool__(self):
        # Return false if the species does not exist in reaction
        if self.species_string == '$Null':
            return False

        species_string_split = self.species_string.split('&')
        if all(char in species_string_split for char in self._stocked_characteristics):
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

    def add(self, characteristic):
        self._stocked_characteristics.add(characteristic)


class IsInstance(Bool_Override):
    """
        Check if a specific species references another.
        Example:
            Something = Age*Live
            IsInstance(Something, Age) = TRUE
            IsInstance(Something, Live) = TRUE
            IsInstance(Something, What) = FALSE
        It does this by overriding the Boolean method after an object has been constructed
    """

    @staticmethod
    def contains_reference(sso, reference):

        if reference in sso._species_object.get_references():
            return True
        else:
            return False

    def __init__(self, sso, reference):
        self._operator = None
        self._reference = reference
        if isinstance(sso, Specific_Species_Operator):
            self._concatenate_sso_object = [sso]
        else:
            simlog.error('Specific_Species_Operator')

    def __bool__(self):
        if not self._operator:
            if len(self._concatenate_sso_object) > 1:
                simlog.error('IsInstance is only used for single species')
            if self.contains_reference(self._concatenate_sso_object[0], self._reference):
                return True


def extract_reaction_rate(combination_of_reactant_species, reactant_string_list
                          , reaction_rate_function, parameters_for_sbml, type_of_model):
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

    if type(reaction_rate_function) == int or type(reaction_rate_function) == float:
        reaction_rate_string = basic_kinetics_string(reactant_string_list,
                                                     reaction_rate_function, parameters_for_sbml, type_of_model)

    elif callable(reaction_rate_function):
        arguments = prepare_arguments_for_callable(combination_of_reactant_species,
                                                   reactant_string_list, reaction_rate_function.__code__.co_varnames)
        rate = reaction_rate_function(**arguments)

        if type(rate) == int or type(rate) == float:
            reaction_rate_string = basic_kinetics_string(reactant_string_list,
                                                         rate, parameters_for_sbml, type_of_model)
        elif type(rate) == str:
            return rate
        elif rate is None:
            simlog.error('There is a reaction rate missing for the following reactants: \n'
                         + str(reactant_string_list))
        else:
            simlog.error('The function return a non-valid value')

    elif type(reaction_rate_function) == str:
        return reaction_rate_function
    elif reaction_rate_function is None:
        simlog.error('There is a reaction rate missing for the following reactants: \n'
                     + str(reactant_string_list))
    else:
        simlog.debug(type(reaction_rate_function))
        simlog.error('The rate type is not supported')

    return reaction_rate_string


def basic_kinetics_string(reactants, reaction_rate, parameters_for_sbml, type_of_model):
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

    rate_str = 'rate_' + str(len(parameters_for_sbml))
    parameters_for_sbml[rate_str] = reaction_rate

    kinetics_string += rate_str

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
