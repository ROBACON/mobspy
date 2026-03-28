"""Reaction-related classes: Reactions, Reacting_Species, and helpers.

Contains _Last_rate_storage, Assignment_Opp_Imp, Reactions, and
Reacting_Species. These are tightly coupled with Species but separated
here for organizational clarity.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Self

from numpy import floating as np_float_
from numpy import integer as np_int_
from pint import Quantity

from mobspy.exceptions import ReactionError, ValidationError
from mobspy.modules.assignments_implementation import (
    Asg as asgi_Asg,
)
from mobspy.modules.assignments_implementation import (
    Assign as asgi_Assign,
)
from mobspy.modules.logic_operator_objects import (
    ReactingSpeciesComparator as lop_ReactingSpeciesComparator,
)
from mobspy.modules.meta_class_utils import (
    unite_characteristics as mcu_unite_characteristics,
)
from mobspy.modules.mobspy_expressions import (
    ExpressionDefiner as me_ExpressionDefiner,
)
from mobspy.modules.mobspy_expressions import (
    OverrideQuantity as me_OverrideQuantity,
)
from mobspy.modules.mobspy_parameters import (
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)

if TYPE_CHECKING:
    from mobspy.modules.species import Species


class _Last_rate_storage:
    last_rate: Any = None
    entity_counter: int = 0

    @classmethod
    def override_get_item(cls, object_to_return: Any, item: Any) -> Any:
        """Store the rate before the reaction is fully defined.

        Due to priority in Python the item is stored before the
        reaction. It is stored at the Compiler level and passed
        to the reaction object in the end.

        :param object_to_return: returns the object that
            __getitem__ was performed on
        :param item: stored reaction rate
        """
        cls.last_rate = item
        return object_to_return

    @classmethod
    def process_rate(cls, rate: Any) -> Any:
        from mobspy.modules.species import Species

        if isinstance(rate, (np_int_, np_float_)):
            rate = float(rate)

        if isinstance(rate, (Species, Reacting_Species, Reactions)):
            raise ReactionError(
                "Reaction rate of type "
                + str(type(_Last_rate_storage.last_rate))
                + " not valid"
            )

        if not (
            isinstance(rate, (int, float, str))
            or callable(rate)
            or isinstance(
                rate,
                (
                    me_OverrideQuantity,
                    Quantity,
                    mp_Mobspy_Parameter,
                    me_ExpressionDefiner,
                ),
            )
            or rate is None
        ):
            raise ReactionError(
                "Reaction rate of type " + str(type(rate)) + " not valid"
            )

        return rate


class Reactions:
    """Reaction class storing reactants, products, rate and order.

    Reactions are created with the ``>>`` operator and stored
    in all involved objects.

    :param reactants: list of meta-species used as reactants
    :param products: list of meta-species used as products
    :param order: reaction order operator (Round-Robin default)
    :param rate: reaction rate
    """

    def __init__(
        self,
        reactants: list[dict[str, Any]],
        products: list[dict[str, Any]],
        rate: Any = None,
    ) -> None:
        """Construct a reaction from reactants and products.

        The order and the rate are assigned later by the compiler.

        :param reactants: list of meta-species reactants
        :param products: list of meta-species products
        :raise simlog.error: if defined under simulation context
        :raise simlog.error: if called from Zero to Zero
        """
        from mobspy.modules.species import Species

        for p in products:
            if p["object"].get_name() == "Context_MetaSpecies":
                raise ReactionError("The Any meta-species cannot be used in reactions")

        if asgi_Assign.check_context():
            asgi_Assign.reset_context()
            raise ReactionError(
                "A MobsPy context error has happen. "
                "A reaction was defined with the assignment context activated."
                "The assignment context was deactivated. "
                "Please try to redefine the model"
            )

        # Setup rate - consider reversible reactions
        try:
            flag_len = len(_Last_rate_storage.last_rate) == 2
        except (TypeError, AttributeError):
            flag_len = False

        if flag_len and rate is None:
            store_last_rate = _Last_rate_storage.last_rate[0]
            Reactions(
                reactants=products,
                products=reactants,
                rate=_Last_rate_storage.last_rate[1],
            )
            rate = _Last_rate_storage.process_rate(store_last_rate)

        elif rate is not None and not flag_len:
            rate = _Last_rate_storage.process_rate(rate)

        elif rate is None and not flag_len:
            rate = _Last_rate_storage.process_rate(_Last_rate_storage.last_rate)

        elif rate is not None and flag_len:
            rate = _Last_rate_storage.process_rate(rate)

        else:
            raise ReactionError("To many rate provided")
        self.rate = rate
        if _Last_rate_storage.last_rate is not None:
            _Last_rate_storage.last_rate = None

        if len(Species.meta_specie_named_any_context) != 0:
            for j in Species.meta_specie_named_any_context:
                for r in reactants:
                    r["object"].c(j)
                    r["characteristics"].add(j)
                for p in products:
                    p["object"].c(j)
                    p["characteristics"].add(j)

        try:
            _ = reactants[0]["object"]
        except IndexError:
            try:
                _ = products[0]["object"]
            except IndexError as e:
                raise ReactionError("No Meta-Species detected in the reaction") from e

        if Species.get_simulation_context() is not None:
            raise ReactionError(
                "Reactions cannot be defined under event context. Only species counts"
            )

        self.reactants = reactants
        self.products = products

        self.order = None

        for reactant in reactants:
            reactant["object"].add_reaction(self)
        for product in products:
            product["object"].add_reaction(self)

    @staticmethod
    def __create_reactants_string(list_of_reactants: list[dict[str, Any]]) -> str:
        """Format a list of reactants/products as a string."""
        reaction_string = ""
        for i, r in enumerate(list_of_reactants):
            if r["stoichiometry"] > 1:
                reaction_string += str(r["stoichiometry"]) + "*" + str(r["object"])
            else:
                reaction_string += str(r["object"])

            if len(r["characteristics"]) > 0:
                reaction_string += "." + ".".join(r["characteristics"])

            if i != len(list_of_reactants) - 1:
                reaction_string += " "
                reaction_string += "+"
                reaction_string += " "

        return reaction_string

    def __str__(self) -> str:
        """Print the meta-reaction in the format A + B -> C + D."""
        return (
            self.__create_reactants_string(self.reactants)
            + " -> "
            + self.__create_reactants_string(self.products)
        )

    def __getitem__(self, item: Any) -> Self:
        """Override of __getitem__ for dealing with reaction rates.

        :param item: (int, float, callable, Quantity) = reaction rate
        """
        return _Last_rate_storage.override_get_item(self, item)  # type: ignore[no-any-return]

    def set_rate(self, rate: Any) -> None:
        """Set the stored reaction rate.

        :param rate: (int, float, callable, Quantity) = reaction rate
        """
        self.rate = rate


class Assignment_Opp_Imp:
    def __add__(self, other: Any) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.add(self, other)
        else:
            raise ValidationError(
                "Addition not implemented for meta-species in this context"
            )

    def __radd__(self, other: Any) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.add(other, self)
        else:
            raise ValidationError(
                "Addition not implemented for meta-species in this context"
            )

    def __sub__(self, other: Any) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.sub(self, other)
        else:
            raise ValidationError(
                "Subtraction not implemented for meta-species in this context"
            )

    def __rsub__(self, other: Any) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.sub(other, self)
        else:
            raise ValidationError(
                "Subtraction not implemented for meta-species in this context"
            )

    def __truediv__(self, other: Any) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.div(self, other)
        else:
            raise ValidationError(
                "Division not implemented for meta-species in this context"
            )

    def __rtruediv__(self, other: Any) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.div(other, self)
        else:
            raise ValidationError(
                "Division not implemented for meta-species in this context"
            )

    def __pow__(self, other: Any) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.pow(self, other)
        else:
            raise ValidationError(
                "Division not implemented for meta-species in this context"
            )

    def __rpow__(self, other: Any) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.pow(other, self)
        else:
            raise ValidationError(
                "Division not implemented for meta-species in this context"
            )

    def __mul__(self, other: Any) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.mul(self, other)
        else:
            raise ValidationError(
                "Multiplication not implemented for meta-species in this sense"
            )

    def __rmul__(self, other: Any) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.mul(other, self)
        else:
            raise ValidationError(
                "Multiplication not implemented for meta-species in this sense"
            )


class Reacting_Species(lop_ReactingSpeciesComparator, Assignment_Opp_Imp):
    """Intermediary object created when a species enters a reaction.

    Transforms a Species object into a list-compatible format
    for the reaction object. The ``>>`` operator calls the
    Reaction constructor.

    Attributes:
        list_of_reactants: Represents reactants or products.
            Each dict contains 'object', 'characteristics',
            'stoichiometry', and optionally 'label'.
    """

    def __init__(
        self,
        object_reference: Species,
        characteristics: set[str],
        stoichiometry: int | float = 1,
        label: int | float | str | None = None,
    ) -> None:
        """Construct a Reacting_Species.

        :param object_reference: meta-species object reference
        :param characteristics: characteristics used to query
        :param stoichiometry: stoichiometry in the reaction
        :param label: label value for matching
        """
        super().__init__()
        if object_reference.get_name() == "_S0" and characteristics == set():
            self.list_of_reactants: list[dict[str, Any]] = []
        else:
            self.list_of_reactants = [
                {
                    "object": object_reference,
                    "characteristics": characteristics,
                    "stoichiometry": stoichiometry,
                    "label": label,
                }
            ]

    def __enter__(self) -> int:
        """Enter context manager for characteristics."""
        self.context_initiator_for_reacting_specie()
        return 0

    def __exit__(self, *args: Any) -> None:
        """Exit context manager for characteristics."""
        self.context_finish_for_reacting_specie()

    def __str__(self) -> str:
        """String representation of the list of reactants."""
        from mobspy.modules.species import Species

        species_object = self.list_of_reactants[0]["object"]
        characteristics = self.list_of_reactants[0]["characteristics"]
        if len(self.list_of_reactants) == 1:
            if Species.get_simulation_context() is None:
                to_return = str(species_object)
                for cha in self.list_of_reactants[0]["characteristics"]:
                    to_return += "." + cha
                return to_return
            else:
                return Species.str_under_context(species_object, characteristics)
        else:
            if Species.get_simulation_context() is not None:
                raise ReactionError(
                    "Please separate the species when using "
                    "string based assignments under event "
                    "context. Ex: str(A) + str(B)"
                )
            return str(self.list_of_reactants)

    def c(self, item: Any) -> Reacting_Species:
        """Query by value instead of name.

        Calls __getattr__ with the value inside the variable.

        :param item: value to query over
        """
        from mobspy.modules.species import Species

        item = str(item)
        Species.check_if_valid_characteristic(self, item)
        return self.__getattr__(item)  # type: ignore[no-any-return]

    def label(self, label: int | float | str) -> Self:
        """Assign a label to a meta-species for compiler matching.

        :param label: value for the label for matching
        """
        if len(self.list_of_reactants) == 1:
            self.list_of_reactants[0]["label"] = label
        else:
            raise ReactionError(
                "Labels cannot be assigned to multiple "
                "reacting species at the same time."
            )
        return self

    def __getitem__(self, item: Any) -> Self:
        """Override of __getitem__ for dealing with reaction rates.

        :param item: (int, float, callable, Quantity) = reaction rate
        """
        return _Last_rate_storage.override_get_item(self, item)  # type: ignore[no-any-return]

    def get_spe_object(self) -> Species:
        if len(self.list_of_reactants) != 1:
            raise ReactionError(
                "The internal method get_queried_characteristics can only be used for "
                "Reacting_Species with a single "
            )
        return self.list_of_reactants[0]["object"]  # type: ignore[no-any-return]

    def get_query_characteristics(self) -> Any:
        if len(self.list_of_reactants) != 1:
            raise ReactionError(
                "The internal method get_queried_characteristics can only be used for "
                "Reacting_Species with a single "
            )
        return self.list_of_reactants[0]["characteristics"]

    def __rmul__(self, stoichiometry: Any) -> Self | Any:
        """Multiply by stoichiometry for reactions.

        :param stoichiometry: stoichiometry value
        """
        if not asgi_Assign.check_context():
            if isinstance(stoichiometry, (int, float)):
                self.list_of_reactants[0]["stoichiometry"] = stoichiometry
            else:
                raise ReactionError(
                    "Stoichiometry can only be an int or "
                    f"float - Received {stoichiometry}"
                )
            return self
        else:
            return asgi_Assign.mul(stoichiometry, self)

    def __add__(self, other: Species | Reacting_Species | Any) -> Self | Any:
        """Addition of meta-species to construct the reaction.

        :param other: (Species or Reacting Species) other object being added
        """
        from mobspy.modules.species import Species

        if not asgi_Assign.check_context():
            if isinstance(other, Species):
                other = Reacting_Species(other, set())
            try:
                self.list_of_reactants += other.list_of_reactants
            except AttributeError as e:
                raise ReactionError(
                    "Addition between meta-species and "
                    f"types {type(other)} is not supported"
                ) from e
            return self
        else:
            return asgi_Assign.add(self, other)

    def __radd__(self, other: Any) -> Self | Any:
        if not asgi_Assign.check_context():
            return Reacting_Species.__add__(self, other)
        else:
            return asgi_Assign.add(other, self)

    def __invert__(self) -> Reacting_Species:
        return self.c("not$")

    def __neg__(self) -> Any:
        if asgi_Assign.check_context():
            return asgi_Assign.mul(-1, self)
        else:
            raise ValidationError(
                "The negative operator was applied to a "
                "Reacting Species in the wrong context"
            )

    def __rshift__(self, other: Species | Reacting_Species) -> Reactions:
        """The ``>>`` operator for defining reactions.

        :param other: (Species or Reacting Species)
            product side of the reaction being added
        """
        import sys

        from mobspy.modules.species import Species, _get_multiline_code_context

        frame = sys._getframe(1)
        code_line = _get_multiline_code_context(
            frame.f_code.co_filename, frame.f_lineno
        )
        line_number = frame.f_lineno
        Species._compile_defined_reaction(code_line, line_number)

        if isinstance(other, Species):  # noqa: SIM108
            p = Reacting_Species(other, set())
        else:
            p = other

        reaction = Reactions(self.list_of_reactants, p.list_of_reactants)
        return reaction

    def __call__(self, quantity: Any) -> Self | None:  # type: ignore[return]
        """Assign counts to species non-default state.

        :param quantity: (int, float, Quantity) count to be assigned
        """
        from mobspy.modules.species import Species

        if isinstance(quantity, (np_int_, np_float_)):
            quantity = float(quantity)

        if len(Species.meta_specie_named_any_context) > 0:
            for i in Species.meta_specie_named_any_context:
                self = self.c(i)  # type: ignore[assignment]

        species_object = self.list_of_reactants[0]["object"]
        characteristics = self.list_of_reactants[0]["characteristics"]
        simulation_under_context = self.list_of_reactants[0][
            "object"
        ].get_simulation_context()
        if (
            isinstance(quantity, (int, float, Quantity, mp_Mobspy_Parameter))
        ) and not asgi_Assign.check_context():
            if len(self.list_of_reactants) != 1:
                raise ReactionError(
                    "Assignment used incorrectly. Only one species at a time"
                )
            quantity_dict = species_object.add_quantities(characteristics, quantity)
        elif asgi_Assign.check_context():
            dummy_rsp = species_object
            dummy_rsp = [dummy_rsp.c(char) for char in characteristics][-1]
            dummy_rsp.assign(quantity)
        elif simulation_under_context is None:
            raise ReactionError(
                "Reactant_Species count assignment does "
                f"not support the type {type(quantity)}"
            )

        if simulation_under_context is not None:
            try:
                if isinstance(quantity, str):
                    quantity_dict = species_object.add_quantities(
                        characteristics, quantity
                    )
                simulation_under_context.current_event_count_data.append(
                    {
                        "species": species_object,
                        "characteristics": quantity_dict["characteristics"],
                        "quantity": quantity_dict["quantity"],
                    }
                )
            except Exception as e:
                raise ReactionError(
                    str(e)
                    + "\n Only species count assignments are allowed in a model context"
                ) from e
        else:
            return self  # pyright: ignore[reportReturnType]

    def __getattr__(self, characteristic: str) -> Any:
        """Implementation of the .dot operation.

        :param characteristic: (str) characteristic for the query
        """
        from mobspy.modules.species import Species

        if characteristic == "_ipython_canary_method_should_not_exist_":
            return 0

        Species.check_if_valid_characteristic(self, characteristic)

        if characteristic == "assign":
            return asgi_Asg(self, species_or_reacting=False)

        for reactant in self.list_of_reactants:
            species_object = reactant["object"]
            characteristics_from_references = mcu_unite_characteristics(
                species_object.get_references()
            )

            if (
                characteristic not in characteristics_from_references
                and "$" not in characteristic
            ):
                if len(species_object.get_characteristics()) == 0:
                    species_object.first_characteristic = characteristic

                species_object.add_characteristic(characteristic)

            reactant["characteristics"].add(characteristic)

        return self

    @classmethod
    def is_species(cls) -> bool:
        return False

    @classmethod
    def is_spe_or_reac(cls) -> bool:
        return True

    old_context: set[str] = set()

    def context_initiator_for_reacting_specie(self) -> None:
        """Add the current context and update the Cts context in all meta-species."""
        from mobspy.modules.species import Species

        if len(self.list_of_reactants) == 1:
            self.old_context = Species.meta_specie_named_any_context
            new_context = Species.meta_specie_named_any_context.union(
                self.list_of_reactants[0]["characteristics"]
            )
            Species.update_meta_specie_named_any_context(new_context)
        else:
            raise ReactionError(
                "Contexts can only be used on basic Reacting meta species"
            )

    def context_finish_for_reacting_specie(self) -> None:
        """Remove the ending context and update."""
        from mobspy.modules.species import Species

        Species.update_meta_specie_named_any_context(self.old_context)


_methods_Reacting_Species = set(dir(Reacting_Species))
