"""
    The class_of_meta_specie_named_any.py model is responsible for defining the Context_specie_named_any class.
"""
from mobspy.modules.meta_class import Species
from inspect import stack as inspect_stack
from mobspy.simulation_logging.log_scripts import error as simlog_error

class Context_specie_named_any(Species):
    """
        Class which inherits from Species. It only has one object, Any, which is defined at the end of this script.
        It is used to simplify the syntax of reactions and count setting inside the body of a
        "with Any.example_characteristic1 :" statement.
        Is is compatible with event_condition and event_time methods of the Simulation class. 
        It can be nested and used several time in the same with statement.
    """
    # This is the set of characteristics which are currently under the active Any context.
    _set_of_characteristics_currently_under_the_any_context = set()

    # This is the list of all the nested Any contexts which are currently active.
    _list_of_nested_any_contexts = []

    def __getattr__(self, item):
        """
            This method is called when an attribute is called on the Any specie.
            It is used to add the characteristic to the currently active Any context.

            :param item: (str) characteristic to be added to the currently active Any context.
            :raise simlog.error: if the Any specie is given a characteristic outside of a context.
            :return self: to allow for assigning multiple characteristics to the Any specie in the same line.
        """
        code_line = inspect_stack()[1].code_context[0][:-1]
        code_line = code_line.split(' ')
        is_with = [x for x in code_line if x != ''][0]
        if is_with == 'with' or item == "__sphinx_mock__":
            pass
        else:
            simlog_error('Characteristics cannot be added to the Any specie outside of a context')
        self._set_of_characteristics_currently_under_the_any_context.add(item)
        return self

    def __enter__(self):
        """
            Context manager for Any's characteristics.
            Called in "with Any.example_characteristic :" format, when entering.
        """
        self.context_initiator_for_meta_specie_named_any()
        return 0

    def __exit__(self, *args):
        """
            Context manager for Any's characteristics.
            Called in "with Any.example_characteristic :" format, when exiting.
        """
        self.context_finish_for_meta_specie_named_any()

    def context_initiator_for_meta_specie_named_any(self):
        """
            This adds the current context in _list_of_nested_any_contexts,
            and then updates the Any context in all meta-species.
        """
        self._list_of_nested_any_contexts.append(set(self._set_of_characteristics_currently_under_the_any_context))
        Species.update_meta_specie_named_any_context(Species.meta_specie_named_any_context.union(self._set_of_characteristics_currently_under_the_any_context))

    def context_finish_for_meta_specie_named_any(self):
        """
            This removes the context which is ending from _list_of_nested_any_contexts and updates the current Any context.
            Then, it updates the Any context in all meta-species.
        """
        self._previous_set_of_characteristics_under_the_any_context = self._list_of_nested_any_contexts.pop()
        if len(self._list_of_nested_any_contexts) > 0:
            self._set_of_characteristics_currently_under_the_any_context = self._list_of_nested_any_contexts[-1]
        else:
            self._set_of_characteristics_currently_under_the_any_context = set()
        Species.update_meta_specie_named_any_context(Species.meta_specie_named_any_context - self._previous_set_of_characteristics_under_the_any_context)

    def __call__(self, quantity):
        """
            The call operator is overloaded as the Any specie cannot be called, as it is not supposed to
            be used this way. Thus, it raises an error when called.
        """
        simlog_error('The Any specie cannot be called')
    
    def __add__(self, other):
        """
            The add operator is overloaded as the Any specie cannot be added. Thus, it raises an error when added.
        """
        simlog_error('The Any specie cannot be added')
    
    def __radd__(self, other):
        """
            The add operator is overloaded as the Any specie cannot be added. Thus it raises an error when added.
        """
        simlog_error('The Any specie cannot be added')
    
    def __rmul__(self, other):
        """
            The multiplication operator is overloaded as the Any specie cannot be multiplied.
            Thus it raises an error when multiplied.
        """
        simlog_error('The Any specie cannot be multiplied')

    def __rshift__(self, other):
        """
            The >> operator is overloaded so that  it raises an error when used.
        """
        simlog_error('The >> operator cannot be used on the Any specie')

   
# Any is the only object of the Any_specie class that will be used. It is defined here.
__SAny = Context_specie_named_any('Context_Any_MetaSpecies')
Any = __SAny
