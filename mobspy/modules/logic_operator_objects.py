import mobspy.simulation_logging.log_scripts as simlog


class SpeciesComparator:

    def __init__(self):
        self._simulation_context = None

    def start_check(self):
        if self._simulation_context is not None:
            self._simulation_context.event_context_add()

    def add_operation_and_number(self, symbol, number):
        operation = [{'object': self, 'characteristics': set()}, symbol, number]
        logic_re_object = MetaSpeciesLogicResolver(operation, self._simulation_context)

        if self._simulation_context is not None:
            self._simulation_context.number_of_context_comparisons += 1
            self._simulation_context.trigger_list.append(logic_re_object)

        return logic_re_object

    def __lt__(self, number):
        self.start_check()
        return self.add_operation_and_number('<', number)

    def __le__(self, number):
        self.start_check()
        return self.add_operation_and_number('<=', number)

    def __gt__(self, number):
        self.start_check()
        return self.add_operation_and_number('>', number)

    def __ge__(self, number):
        self.start_check()
        return self.add_operation_and_number('>=', number)

    def __eq__(self, other):
        if self._simulation_context is not None:
            simlog.error('Equality assignment not allowed for event condition in MobsPy.\n'
                         'Please if necessary use ( >= ) & ( =< )')
        else:
            return id(self) == id(other)

    def __neg__(self, other):
        return id(self) != id(other)

    def __hash__(self):
        return hash(id(self))


class ReactingSpeciesComparator(SpeciesComparator):

    def __init__(self):
        super(ReactingSpeciesComparator, self).__init__()
        self._simulation_context = None

    def start_check(self):
        self.check_context()
        if self._simulation_context is not None:
            self._simulation_context.event_context_add()

    def check_context(self):
        for react_dict in self.list_of_reactants:
            if react_dict['object']._simulation_context is not None:
                self._simulation_context = react_dict['object']._simulation_context
                break
        else:
            self._simulation_context = None

    def add_operation_and_number(self, symbol, number):
        operation = []

        for i, react_dict in enumerate(self.list_of_reactants):
            if i > 0:
                dl = ['+']
            else:
                dl = []
            dl = dl + [react_dict['stoichiometry'], '*']
            dl = dl + [{'object': react_dict['object'], 'characteristics': react_dict['characteristics']}]
            operation = operation + dl
        operation = operation + [symbol, number]
        logic_re_object = MetaSpeciesLogicResolver(operation, self._simulation_context)

        self.check_context()
        if self._simulation_context is not None:
            self._simulation_context.number_of_context_comparisons += 1
            self._simulation_context.trigger_list.append(logic_re_object)

        return logic_re_object


class MetaSpeciesLogicResolver:

    def __init__(self, operation, model_context=None):
        self.operation = operation
        self.model_context = model_context
        if self.model_context is not None:
            model_context.trigger_list.append(self)
        self.order = 0

    def __bool__(self):
        if self.model_context is not None:
            self.model_context.bool_number_call += 1
        return True

    def __and__(self, other):
        if not isinstance(other, MetaSpeciesLogicResolver):
            print('ERROR')
            exit()

        new_operation = ['('] + ['('] + self.operation + [')'] + ['&&'] + ['('] + other.operation + [')'] + [')']
        self.operation = new_operation
        self.order += 1
        if self.model_context is not None:
            self.model_context.trigger_list.append(self)
        return self

    def __or__(self, other):
        if not isinstance(other, MetaSpeciesLogicResolver):
            print('ERROR')
            exit()

        new_operation = ['('] + ['('] + self.operation + + [')'] + ['||'] + ['(']  + other.operation + [')'] + [')']
        self.operation = new_operation
        self.order += 1
        if self.model_context is not None:
            self.model_context.trigger_list.append(self.operation)
        return self

    @classmethod
    def find_all_species_strings(cls, species, characteristics, species_for_sbml):
        reference_set = set(characteristics)
        reference_set.add(str(species))
        return [x for x in species_for_sbml.keys() if reference_set.issubset(set(x.split('_dot_')))]

    def generate_string_from_vec_space(self, species_for_sbml):
        copasi_str = ''
        for i, e in enumerate(self.operation):
            if type(e) == int or type(e) == float or type(e) == str:
                copasi_str = copasi_str + str(e) + ' '
            else:
                copasi_str = copasi_str + '('
                ite = self.find_all_species_strings(e['object'], e['characteristics'], species_for_sbml)
                for j, species_str in enumerate(ite):
                    if j == 0:
                        copasi_str = copasi_str + f'{species_str}'
                    else:
                        copasi_str = copasi_str + f' + {species_str}'
                copasi_str = copasi_str + ')' + ' '

        return copasi_str
