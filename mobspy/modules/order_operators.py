"""This module deals with order_operators.

Thus it deals with the Round-Robin order assignment and with
species in the products that are not referenced in the
reactants (Born Species).
It has also the implementation of the Rev operator, since
that is a also a reaction operator.
"""

from __future__ import annotations

import warnings
from copy import deepcopy
from typing import Any, Protocol, runtime_checkable

from mobspy.constants import ALL_CHAR, DOT_SEPARATOR, SET_SPECIES_PREFIX
from mobspy.exceptions import ReactionError
from mobspy.modules.reactions import Reactions
from mobspy.modules.species import Species
from mobspy.modules.species_string_generator import (
    construct_all_combinations as ssg_construct_all_combinations,
)
from mobspy.modules.species_string_generator import (
    construct_species_char_list as ssg_construct_species_char_list,
)

_REVERSIBLE_RATE_PAIR_LEN = 2


@runtime_checkable
class OrderOperator(Protocol):
    """Protocol for reaction order operators (All, Default).

    Allows type-annotating parameters that accept order operators::

        def my_func(order: OrderOperator) -> None: ...
        my_func(All)   # ok
        my_func(Default)  # ok
    """

    def __getitem__(self, item: Any) -> Any: ...

    def __call__(
        self,
        order_dictionary: dict[tuple[Any, Any], list[Any]],
        product_species: list[dict[str, Any]],
        model: list[Any] | set[Any],
        ref_characteristics_to_object: dict[str, Any],
        all_reactions: bool = False,
    ) -> list[list[tuple[Any, Any]]]: ...


class __Operator_Base:  # noqa: N801
    """Order operator base with shared utilities.

    Contains functions useful for the reaction operators.
    Reaction Operators are responsible for order assignment
    and how to deal with species in the products that are
    not referenced in the reactants (Born Species).
    It also overrides the getitem method so we can use:
    Operator[Reaction] to assign an order to the reaction.
    """

    # Assign order structure
    def __getitem__(self, item: Any) -> Any:
        try:
            if isinstance(item, str):
                return item + "." + ALL_CHAR

            if item.is_species():
                return item.c(ALL_CHAR)
            if not item.is_species():
                for reactant in item.list_of_reactants:
                    reactant["characteristics"].add(ALL_CHAR)
                return item
        except AttributeError as e:
            raise ReactionError(
                "All can only be used on species, reacting"
                " species and strings under set_count"
            ) from e
        return None

    # Transform product function
    @staticmethod
    def transform_species_string(
        species_string: list[Any],
        characteristics_to_transform: set[str] | str,
        ref_characteristics_to_object: dict[str, Any],
    ) -> str | tuple[str, float]:
        """Handle characteristic change from a chemical reaction.

        Receives a string in the MobsPy format and the
        characteristic written on the product side. MobsPy
        converts the characteristics in the same position as
        that one (directly added to the same meta-species) to
        the one in characteristics_to_transform.

        Args:
            species_string: Species string in MobsPy format.
            characteristics_to_transform: Characteristics received by the species object
                in the product.
            ref_characteristics_to_object: Characteristics as keys and meta-species
                objects as values.
        """
        species_object = species_string[0]
        species_to_return_list: list[Any] = (
            [species_object.get_name(), *deepcopy(species_string[1:])]
            if len(species_string) > 1
            else [species_object.get_name()]
        )

        for characteristic in characteristics_to_transform:
            obj = ref_characteristics_to_object[characteristic]
            i = species_object.get_index_from_reference_dict(obj)
            species_to_return_list[i] = characteristic

        species_to_return: str | tuple[str, float]
        if isinstance(species_to_return_list[-1], str):
            species_to_return = DOT_SEPARATOR.join(species_to_return_list)
        elif isinstance(species_to_return_list[-1], float):
            species_to_return = (
                DOT_SEPARATOR.join(species_to_return_list[:-1]),
                species_to_return_list[-1],
            )
        else:
            species_to_return = DOT_SEPARATOR.join(species_to_return_list)

        return species_to_return

    @staticmethod
    def find_all_string_references_to_born_species(
        species_referenced_by: list[Any],
        characteristics: set[str] | str,
        ref_characteristics_to_object: dict[str, Any],
    ) -> list[Any]:
        """Find all string species referencing the born-meta species.

        Includes inheritors. Returns them all in a list
        for the product construction.

        Args:
            species_referenced_by: (meta-species list) Species that inherit from the
                born-meta species in the reaction.
            characteristics: Characteristics queried on the born species in the product.
            meta_species_in_model: Meta-species in model as values.
        """
        to_return: list[Any] = []
        for species in species_referenced_by:
            to_return += ssg_construct_all_combinations(
                species,
                characteristics,
                ref_characteristics_to_object,
                symbol=DOT_SEPARATOR,
            )

        return to_return

    @staticmethod
    def find_all_default_references_to_born_species(
        species_referenced_by: list[Any],
        characteristics: set[str] | str,
        ref_characteristics_to_object: dict[str, Any],
    ) -> list[Any]:
        """Find DEFAULT string species referencing born-meta species.

        Includes inheritors. Returns them all in a list
        for the product construction.

        Args:
            species_referenced_by: (meta-species list) Species that inherit from the
                born-meta species in the reaction.
            characteristics: Characteristics queried on the born species in the product.
            meta_species_in_model: Meta-species in model.
            ref_characteristics_to_object: Characteristics as keys and meta-species
                objects directly added to as values.
        """
        to_return: list[Any] = []
        for species in species_referenced_by:
            to_return += [
                ssg_construct_species_char_list(
                    species,
                    characteristics,
                    ref_characteristics_to_object,
                    symbol=DOT_SEPARATOR,
                )
            ]

        return to_return

    def __call__(
        self,
        order_dictionary: dict[tuple[Any, Any], list[Any]],
        product_species: list[dict[str, Any]],
        model: list[Any] | set[Any],
        ref_characteristics_to_object: dict[str, Any],
        all_reactions: bool = False,
    ) -> list[list[tuple[Any, Any]]]:
        """Parry products with reactants by meta-species.

        Uses the order_dictionary and an index-system for
        round-robin application. If it cannot find a
        parrying it considers it a born species.

        Args:
            order_dictionary: Meta-species as keys and list of meta-species strings.
            product_species: Product species dict with meta-species object (key:
                species), label (key: label), and characteristics (key:
                characteristics).
            ref_characteristics_to_object: Characteristics as keys and meta-species
                objects directly added to as values.
        """
        round_robin_index: dict[tuple[Any, Any], int] = {}
        for species, label in [(e["species"], e["label"]) for e in product_species]:
            round_robin_index[(species, label)] = 0

        products: list[list[tuple[Any, Any]]] = []
        for species, label, characteristics, stoichiometry in [
            (e["species"], e["label"], e["characteristics"], e["stoichiometry"])
            for e in product_species
        ]:
            if ALL_CHAR in characteristics:
                species_is_referenced_by: list[Any] = [
                    spe_obe for spe_obe in model if species in spe_obe.get_references()
                ]
                all_strings = self.find_all_string_references_to_born_species(
                    species_is_referenced_by,
                    characteristics,
                    ref_characteristics_to_object,
                )
                products.append([(stoichiometry, s) for s in all_strings])
                continue

            # Simple round robin
            try:
                species_to_transform_string = order_dictionary[(species, label)][
                    round_robin_index[(species, label)]
                ]
                round_robin_index[(species, label)] = (
                    round_robin_index[(species, label)] + 1
                ) % len(order_dictionary[(species, label)])

                # Return in list of lists format for combination later
                products.append(
                    [
                        (
                            stoichiometry,
                            self.transform_species_string(
                                species_to_transform_string,
                                characteristics,
                                ref_characteristics_to_object,
                            ),
                        )
                    ]
                )

            # If the species is not on the reactants - order_dictionary
            except KeyError:
                species_is_referenced_by = []
                for spe_obe in model:
                    if species in spe_obe.get_references():
                        species_is_referenced_by.append(spe_obe)

                if len(species_is_referenced_by) == 0:
                    raise ReactionError(
                        f"Species {species} was used in a "
                        "reaction but itself or any inheritors"
                        " are not in the model."
                        " Please add at least one"
                    ) from None

                if all_reactions:
                    # Find all the species that reference the one in the reaction
                    products.append(
                        [
                            (stoichiometry, s)
                            for s in self.find_all_string_references_to_born_species(
                                species_is_referenced_by,
                                characteristics,
                                ref_characteristics_to_object,
                            )
                        ]
                    )
                else:
                    # Find only default state of the species
                    products.append(
                        [
                            (stoichiometry, s)
                            for s in self.find_all_default_references_to_born_species(
                                species_is_referenced_by,
                                characteristics,
                                ref_characteristics_to_object,
                            )
                        ]
                    )

        return products


class __Round_Robin_Base(__Operator_Base):  # noqa: N801
    """Here we have the implementation of the round robin order it goes like this:

    2*Ecoli >> 4*Ecoli

    Since an Ecoli is a generic object that refers to all types of Ecoli
    This object refers to multiple reactions and to multiple strings
    Therefore we follow this rule to know which string goes where:
    1 Ecoli reactant >> 1 Ecoli product
    2 Ecoli reactant >> 2 Ecoli product
    1 Ecoli reactant >> 3 Ecoli product
    2 Ecoli reactant >> 4 Ecoli product
    This is a cycle (round robin) between the reactants
    to be assigned positions in the product.
    The products will keep the characteristics of the
    reactants except if stated otherwise with dot notation.
    For completely new species (no reactant of the same
    species) we use ALL possible combinations.
    For only the default option see the code below.
    """

    def __call__(
        self,
        order_dictionary: dict[tuple[Any, Any], list[Any]],
        product_species: list[dict[str, Any]],
        meta_species_in_model: list[Any] | set[Any],
        ref_characteristics_to_object: dict[str, Any],
        all_reactions: bool = True,
    ) -> list[list[tuple[Any, Any]]]:
        """Parry products with reactants by meta-species.

        Uses the order_dictionary and an index-system for
        round-robin application. If it cannot find a
        parrying it considers it a born species.

        Args:
            order_dictionary: Meta-species as keys and list of meta-species strings.
            product_species: Product species dict with meta-species object (key:
                species), label (key: label), and characteristics (key:
                characteristics).
            ref_characteristics_to_object: Characteristics as keys and meta-species
                objects directly added to as values.
        """
        return super().__call__(
            order_dictionary,
            product_species,
            meta_species_in_model,
            ref_characteristics_to_object,
            all_reactions,
        )


# Define class to override the operators
All = __Round_Robin_Base()


class __RR_Default_Base(__Operator_Base):  # noqa: N801
    """Only default options for born species (no reference in reactant)
    See __Round_Robin_Base for clarification
    """

    # Here is the default order requested by Thomas
    def __call__(
        self,
        order_dictionary: dict[tuple[Any, Any], list[Any]],
        product_species: list[dict[str, Any]],
        meta_species_in_model: list[Any] | set[Any],
        ref_characteristics_to_object: dict[str, Any],
        all_reactions: bool = False,
    ) -> list[list[tuple[Any, Any]]]:
        """Parry products with reactants by meta-species.

        Uses the order_dictionary and an index-system for
        round-robin application. If it cannot find a
        parrying it considers it a born species.

        Args:
            order_dictionary: Meta-species as keys and list of meta-species strings.
            product_species: Product species dict with meta-species object (key:
                species), label (key: label), and characteristics (key:
                characteristics).
            ref_characteristics_to_object: Characteristics as keys and meta-species
                objects directly added to as values.
        """
        return super().__call__(
            order_dictionary,
            product_species,
            meta_species_in_model,
            ref_characteristics_to_object,
            all_reactions,
        )


Default = __RR_Default_Base()


class __Set_Reversible_Rate:  # noqa: N801
    """This class is responsible for dealing with reversible reactions.
    Since reversible reactions have two rates, the process is in __getitem__.
    """

    def __getitem__(self, both_rates: tuple[Any, Any]) -> None:
        """This function extracts both reaction rates from a tuple

        Args:
            both_rates: Rate function from the direct and reverse reaction.
        """
        try:
            if len(both_rates) != _REVERSIBLE_RATE_PAIR_LEN:
                raise ReactionError("The reversible reaction must receive 2 rates")
        except TypeError as e:
            raise ReactionError("The reversible reaction must receive 2 rates") from e

        self.reaction_direct.rate = both_rates[0]
        self.reaction_reverse.rate = both_rates[1]

    def __init__(self) -> None:
        """Dummy constructor for rate setter."""
        self.reaction_direct: Any = None
        self.reaction_reverse: Any = None

    def set_reactions(self, reaction_direct: Any, reaction_reverse: Any) -> None:
        """Store the forward and reverse reaction pair."""
        self.reaction_direct = reaction_direct
        self.reaction_reverse = reaction_reverse


# Reversible reaction operator
class __Reversible_Base:  # noqa: N801
    def __getitem__(self, reaction: Reactions) -> __Set_Reversible_Rate:
        """Set reversible reaction rates via rate setter.

        Uses __Set_Reversible_Rate instance to set rates.
        Also creates the reverse reaction.

        .. deprecated::
            Use ``A >> B @ (k_fwd, k_rev)`` instead.

        Args:
            reaction: Object from the reaction class.
        """
        warnings.warn(
            "Rev[reaction][k_fwd, k_rev] is deprecated. "
            "Use 'A >> B @ (k_fwd, k_rev)' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        reaction_direct = reaction
        reaction_reverse = Reactions(
            reaction_direct.products, reaction_direct.reactants
        )
        self.rate_setter.set_reactions(reaction_direct, reaction_reverse)
        return self.rate_setter

    def __init__(self, rate_setter: __Set_Reversible_Rate) -> None:
        """Dummy constructor for reversible base."""
        self.rate_setter = rate_setter


Rev = __Reversible_Base(__Set_Reversible_Rate())
"""Create reversible reactions: ``Rev[A >> B][k_fwd, k_rev]``."""


class _Set_Reaction_User_Base:  # noqa: N801
    def __getitem__(self, reaction: Reactions) -> _Set_Reaction_Method:
        """Assign a reaction to a Set operator.

        Args:
            reaction: Object from the reaction class.
        """
        return _Set_Reaction_Method(reaction)


Set = _Set_Reaction_User_Base()
"""Define ordered reaction sets: ``Set[A >> B][rate1, rate2, ...]``."""


class _Set_Reaction_Method:  # noqa: N801
    _number_of_calls: int = 0

    def __init__(self, reaction: Reactions) -> None:
        if reaction.rate is None:
            raise ReactionError(
                "A reaction rate was not found in the reaction used in the Set Operator"
            )

        self.reaction = reaction

    def at(self, time: int | float) -> _Set_Reaction_Method:  # noqa: ARG002
        """Schedule this reaction to activate at the given time."""
        self._number_of_calls += 1
        self._set_species = Species("dummy")
        self._set_species._bypass_name(SET_SPECIES_PREFIX + str(self._number_of_calls))

        dict_insert_species: dict[str, Any] = {
            "object": self._set_species,
            "characteristics": set(),
            "stoichiometry": 1,
            "label": None,
        }
        self._set_species(1)
        self.reaction.reactants.append(dict_insert_species)
        self.reaction.products.append(dict_insert_species)

        for spe_obj in self.reaction.reactants + self.reaction.products:
            spe_obj = spe_obj["object"]  # noqa: PLW2901
            spe_obj.link_a_species(self._set_species)

        return self
