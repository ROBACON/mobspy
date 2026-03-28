"""Construction of individual reactions from a meta-reaction."""

from __future__ import annotations

from collections.abc import Generator, Sequence
from inspect import signature as inspect_signature
from itertools import product as itertools_product
from typing import TYPE_CHECKING, Any

from mobspy.exceptions import CompilationError
from mobspy.types import ReactionData

if TYPE_CHECKING:
    from mobspy.types import ReactionsForSbml

# from mobspy.modules.mobspy_parameters import *
from mobspy.modules.context_related_scripts import (
    Unit_Context_Setter as crs_Unit_Context_Setter,
)
from mobspy.modules.function_rate_code import (
    extract_reaction_rate as fr_extract_reaction_rate,
)
from mobspy.modules.meta_class import Reactions
from mobspy.modules.meta_class_utils import (
    count_stoichiometry as mcu_count_string_dictionary,
)
from mobspy.modules.order_operators import Default
from mobspy.modules.species_string_generator import (
    construct_all_combinations as ssg_construct_all_combinations,
)


def iterator_for_combinations(
    list_of_lists: list[list[Any]],
) -> Generator[tuple[Any, ...], None, None]:
    """Iterate through all the combinations of a list of lists
    [[A,B,C], [D, F, E], [G]] means:
    ADG, AFG, AEG, BDG, BFG, BEG, CDG, CFG, CEG .....

    Parameter:
        list_of_lists (list of lists) = list of all lists
            to iterate through all combinations
    """
    for i in itertools_product(*list_of_lists):  # noqa: UP028
        yield i


def check_for_invalid_reactions(
    reactions: set[Any],
    ref_characteristics_to_object: dict[str, Any],
) -> None:
    """Check for queries referencing two independent
    characteristics from the same meta-species, which
    would result in an empty set. Like Ecoli.live.dead
    (assuming live and dead belong to the same
    meta-species characteristics set). If that ever
    happens inside a meta-reaction we just pop an error.

    :param reactions: (set of meta-reactions) set of
        meta-reactions inside the model
    :param ref_characteristics_to_object: (dict)
        dictionary with the characteristics as keys
        and their respective object as values
    """
    for reaction in reactions:
        for reactant in reaction.reactants:
            check_for_duplicates: dict[Any, str] = {}
            for cha in reactant["characteristics"]:
                if cha == "all$":
                    continue

                try:
                    check_for_duplicates[ref_characteristics_to_object[cha]]
                    raise CompilationError(
                        f"Illegal reaction: {reaction}. \n"
                        "There is a query with two "
                        "characteristics "
                        f"{ref_characteristics_to_object[cha].get_characteristics()}"
                        " from the same"
                        " vector axis resulting in an"
                        " impossible query. \n"
                        "As all those characteristics"
                        " have been directly added to"
                        f" {ref_characteristics_to_object[cha]},"
                        " they are located in the same"
                        " vector axis. \n"
                        "To solve this either assign the"
                        " characteristics to new"
                        " meta-species and use "
                        "inheritance or create a new"
                        " reaction for each desired"
                        " query."
                    )
                except KeyError:
                    try:
                        check_for_duplicates[ref_characteristics_to_object[cha]] = cha
                    except KeyError as e:
                        raise CompilationError(
                            "A base object for"
                            f" characteristic {cha} was"
                            " not found in the species"
                            " supplied to the "
                            "simulator \n"
                            "Perhaps a species is missing ? "
                        ) from e

        for product in reaction.products:
            check_for_duplicates = {}
            for cha in product["characteristics"]:
                try:
                    check_for_duplicates[ref_characteristics_to_object[cha]]
                    raise CompilationError(
                        f"Illegal reaction: {reaction}. \n"
                        "There is a transformation"
                        " with two characteristics "
                        f"{ref_characteristics_to_object[cha].get_characteristics()}"
                        " from the same"
                        " vector axis resulting in an"
                        " undefinable"
                        " transformation. \n"
                        "As all those characteristics"
                        " have been directly added"
                        " to "
                        f"{ref_characteristics_to_object[cha]},"
                        " they are located in the"
                        " same vector axis. \n"
                        "To solve this either assign"
                        " the characteristics to new"
                        " meta-species and use "
                        "inheritance or create a new"
                        " reaction for each desired"
                        " query."
                    )
                except KeyError:
                    if "$" not in cha:
                        check_for_duplicates[ref_characteristics_to_object[cha]] = cha


def construct_reactant_structures(
    reactant_species: tuple[Any, ...], ref_characteristics_to_object: dict[str, Any]
) -> list[list[Any]]:
    """Find corresponding strings for the reactant
    meta-species with or without a query, pack them in
    a list and return.

    :param reactant_species: (meta-species object)
        species objects of the involved species
    :param ref_characteristics_to_object: (dict)
        dictionary with characteristics as keys and
        species objects as values
    """
    species_string_combinations: list[list[Any]] = []

    for reactant in reactant_species:
        species_string_combinations.append(
            ssg_construct_all_combinations(
                reactant["object"],
                reactant["characteristics"],
                ref_characteristics_to_object,
            )
        )

    return species_string_combinations


def construct_order_structure(
    species_order_list: list[tuple[Any, Any]],
    current_species_string_list: tuple[Any, ...],
) -> dict[Any, list[Any]]:
    """Order structure for reaction order operations.
    Returns the cyclic_dictionary to be used by the
    order operator. The meta-species objects are the
    keys of this dictionary and a lists of species
    strings currently being used in the reaction are
    the values, allowing the product to find its
    corresponding species-string in a future step.

    :param species_order_list: (list of meta-species
        objects) list of meta-species objects as they
        appear in the meta-reaction
    :param current_species_string_list: (list of
        strings) list of strings in MobsPy format of
        the species currently in this specific reaction

    :return: cyclic_dict (dict) Dictionary where the
        keys are meta-species objects and the values
        are lists of species
    """
    cyclic_dict: dict[Any, list[Any]] = {}
    for species_object, species_string in zip(  # noqa: B905
        species_order_list, current_species_string_list
    ):
        try:
            cyclic_dict[species_object].append(species_string)
        except KeyError:
            cyclic_dict[species_object] = [species_string]

    return cyclic_dict


def construct_product_structure(reaction: Reactions) -> list[dict[str, Any]]:
    """This function unpacks the products in a meta-reaction

    :param: reaction meta-reaction currently being analysed

    :return: product_list = A list of dictionaries for
        each product with the meta-species object,
        the label and the characteristics
    """
    product_list: list[dict[str, Any]] = []
    for product in reaction.products:
        if isinstance(product["stoichiometry"], float):
            product_list.append(
                {
                    "species": product["object"],
                    "label": product["label"],
                    "characteristics": product["characteristics"],
                    "stoichiometry": product["stoichiometry"],
                }
            )
        else:
            for _ in range(product["stoichiometry"]):
                product_list.append(
                    {
                        "species": product["object"],
                        "label": product["label"],
                        "characteristics": product["characteristics"],
                        "stoichiometry": 1,
                    }
                )

    return product_list


def construct_single_reaction_for_sbml(
    reactant_species_string_list: list[str],
    product_species_string_list: Sequence[tuple[float, str]],
    reaction_rate: str,
) -> ReactionData:
    """Construct the reactions for SBML for the
    conversion by the model builder script.
    It follows the following structure 're':[('stoichmetry', reactantant_string) ....
    The reaction rate must be a string containing the reaction kinetics
    This returns a single reaction to be appended by the reactions_for_sbml dictionary

    :param reactant_species_string_list: (list of
        strings) list of reactants in MobsPy format
    :param product_species_string_list: (list of
        strings) list of products in MobsPy format
    :param reaction_rate: (str) reaction rate expression as a string

    :return: to_return (dict) = dictionary that packs the reactants products and rate
    """
    to_return = ReactionData(reactants=[], products=[], kinetics=reaction_rate)
    reactant_count_dict = mcu_count_string_dictionary(reactant_species_string_list)
    product_count_dict = mcu_count_string_dictionary(product_species_string_list)

    for key in reactant_count_dict:
        to_return.reactants.append((reactant_count_dict[key], key))

    for key in product_count_dict:
        to_return.products.append((product_count_dict[key], key))

    return to_return


def get_involved_species(
    reaction: Reactions, meta_species_in_model: list[Any]
) -> tuple[list[tuple[Any, Any]], list[list[dict[str, Any]]]]:
    """Extract all involved meta-species inside a
    reaction. This function implements the inheritance
    mechanism by finding within each meta-species
    references set if they reference the meta-species
    in the reaction.

    :param reaction: (meta-reaction object)
    :param meta_species_in_model: (list) list of
        meta-species used in the model

    :return: base_species_order (list of meta-species
        objects) = order that the meta-species appear
        in the meta-reaction,
        reactant_species_combination_list (list of
        lists of meta-species) = list of lists of all
        meta-species that have inherited from the
        meta-species in the meta-reaction
    """
    reactant_species_combination_list: list[list[dict[str, Any]]] = []
    base_species_order: list[tuple[Any, Any]] = []

    for reactant in reaction.reactants:
        flag_absent_reactant = False
        for _ in range(reactant["stoichiometry"]):
            species_for_reactant: list[dict[str, Any]] = []
            base_species_order.append((reactant["object"], reactant["label"]))

            for species in meta_species_in_model:
                if reactant["object"] in species.get_references():
                    species_for_reactant.append(
                        {
                            "object": species,
                            "characteristics": reactant["characteristics"],
                            "stoichiometry": reactant["stoichiometry"],
                        }
                    )
                    flag_absent_reactant = True

            if not flag_absent_reactant:
                raise CompilationError(
                    f"Species {reactant['object']} or any"
                    " inheritors were not found"
                    " in model \n"
                    f"For reaction {reaction} \n"
                    f"Please add the species or remove the reaction"
                )

            reactant_species_combination_list.append(species_for_reactant)

    return base_species_order, reactant_species_combination_list


def construct_rate_function_arguments(
    rate_function: Any,
    reaction: Reactions,
) -> list[str]:
    rate_function_arguments = str(inspect_signature(rate_function))

    black_list = ["*", "="]
    if any(i in rate_function_arguments for i in black_list):
        raise CompilationError(
            f"Rate arguments must not contain = or *. \n"
            f"Error in reaction {reaction}. \n"
            "Error in rate function"
            f" {rate_function} in signature"
            f" {inspect_signature(rate_function)!s}"
        )

    rate_function_arguments = str(rate_function_arguments).replace("(", "")
    rate_function_arguments = str(rate_function_arguments).replace(")", "")
    rate_function_arguments = str(rate_function_arguments).replace(" ", "")
    rate_function_arguments = rate_function_arguments.split(",")
    return rate_function_arguments


def create_all_reactions(
    reactions: set[Any],
    meta_species_in_model: list[Any],
    ref_characteristics_to_object: dict[str, Any],
    type_of_model: str,
    dimension: int,
    parameter_exist: dict[str, Any],
    parameters_in_reaction: dict[str, Any],
    skip_check: bool,
) -> tuple[ReactionsForSbml, dict[str, Any]]:
    """This function creates all reactions
    Returns the reactions_for_sbml and parameters_for_sbml dictionary
    Those will be used by another module to create the SBML file

    :param reactions: (meta-reaction objects) reactions
        objects constructed by the meta_class module
    :param meta_species_in_model: (list) list of
        meta-species in model
    :param ref_characteristics_to_object: (dict)
        Characteristics as keys objects as values
    :param type_of_model: (str) stochastic or deterministic
    :param dimension: (int) model dimension 1D, 2D, 3D, .....

    :returns: reactions_for_sbml (dict) = dictionary
        with all reactions that will be added to the
        sbml model file, parameters_for_sbml (dict) =
        parameters for the sbml model file
    """
    reactions_for_sbml: ReactionsForSbml = {}

    check_for_invalid_reactions(reactions, ref_characteristics_to_object)

    # Initiate expressions
    with crs_Unit_Context_Setter():
        for reaction in reactions:
            base_species_order, reactant_species_combination_list = (
                get_involved_species(reaction, meta_species_in_model)
            )

            for combination_of_reactant_species in iterator_for_combinations(
                reactant_species_combination_list
            ):
                reactant_species_string_combination_list = (
                    construct_reactant_structures(
                        combination_of_reactant_species, ref_characteristics_to_object
                    )
                )

                for reactant_string_list in iterator_for_combinations(
                    reactant_species_string_combination_list
                ):
                    product_object_list = construct_product_structure(reaction)
                    order_structure = construct_order_structure(
                        base_species_order, reactant_string_list
                    )

                    if reaction.order is None:
                        product_species_species_string_combination_list = Default(
                            order_structure,
                            product_object_list,
                            meta_species_in_model,
                            ref_characteristics_to_object,
                        )
                    else:
                        product_species_species_string_combination_list = (
                            reaction.order(
                                order_structure,
                                product_object_list,
                                meta_species_in_model,
                                ref_characteristics_to_object,
                            )
                        )

                    for product_string_list in iterator_for_combinations(
                        product_species_species_string_combination_list
                    ):
                        reaction_rate_arguments = None
                        if callable(reaction.rate):
                            reaction_rate_arguments = construct_rate_function_arguments(
                                reaction.rate, reaction
                            )

                        reactant_strings = [
                            "_dot_".join([reactant[0].get_name()] + reactant[1:])
                            if len(reactant) > 1
                            else reactant[0].get_name()
                            for reactant in reactant_string_list
                        ]

                        try:
                            rate_string, parameters_in_reaction = (
                                fr_extract_reaction_rate(
                                    combination_of_reactant_species,
                                    reactant_strings,
                                    reaction.rate,
                                    type_of_model,
                                    dimension,
                                    reaction_rate_arguments,
                                    parameter_exist,
                                    parameters_in_reaction,
                                    skip_check,
                                )
                            )
                        except TypeError as e:
                            raise CompilationError(
                                f"On reaction {reaction} \n" + str(e)
                            ) from e

                        if rate_string == 0:
                            continue

                        reactions_for_sbml[
                            "reaction_" + str(len(reactions_for_sbml))
                        ] = construct_single_reaction_for_sbml(
                            reactant_strings, product_string_list, rate_string
                        )

    return reactions_for_sbml, parameters_in_reaction


if __name__ == "__main__":
    pass
