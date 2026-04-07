"""Construction of individual reactions from a meta-reaction."""

from __future__ import annotations

from inspect import signature as inspect_signature
from itertools import product as itertools_product
from typing import TYPE_CHECKING, Any

from mobspy.constants import ALL_CHAR
from mobspy.exceptions import CompilationError
from mobspy.types import CompilationContext, ConcreteSpeciesId, ReactionData

if TYPE_CHECKING:
    from collections.abc import Generator, Sequence

    from mobspy.modules.reactions import Reactions
    from mobspy.types import ReactionsForSbml

from mobspy.modules.expression_context import (
    Unit_Context_Setter as crs_Unit_Context_Setter,
)
from mobspy.modules.order_operators import Default
from mobspy.modules.rate_functions import (
    extract_reaction_rate as fr_extract_reaction_rate,
)
from mobspy.modules.species_string_generator import (
    construct_all_combinations as ssg_construct_all_combinations,
)
from mobspy.modules.species_utils import (
    count_stoichiometry as mcu_count_string_dictionary,
)


def iterator_for_combinations(
    list_of_lists: list[list[Any]],
) -> Generator[tuple[Any, ...], None, None]:
    """Iterate through all the combinations of a list of lists
    [[A,B,C], [D, F, E], [G]] means:
    ADG, AFG, AEG, BDG, BFG, BEG, CDG, CFG, CEG .....

    Args:
        list_of_lists: List of all lists to iterate through all combinations.
    """
    yield from itertools_product(*list_of_lists)


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

    Args:
        reactions: (set of meta-reactions) set of meta-reactions inside the model.
        ref_characteristics_to_object: Dictionary with the characteristics as keys and
            their respective object as values.
    """
    for reaction in reactions:
        for reactant in reaction.reactants:
            check_for_duplicates: dict[Any, str] = {}
            for cha in reactant["characteristics"]:
                if cha == ALL_CHAR:
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

    Args:
        reactant_species: (meta-species object) species objects of the involved species.
        ref_characteristics_to_object: Dictionary with characteristics as keys and
            species objects as values.
    """
    return [
        ssg_construct_all_combinations(
            reactant["object"],
            reactant["characteristics"],
            ref_characteristics_to_object,
        )
        for reactant in reactant_species
    ]


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

    Args:
        species_order_list: (list of meta-species objects) list of meta-species objects
            as they appear in the meta-reaction.
        current_species_string_list: List of strings in MobsPy format of the species
            currently in this specific reaction.


    Returns:
        Dictionary where the keys are meta-species objects and the values are lists of
        species.
    """
    cyclic_dict: dict[Any, list[Any]] = {}
    for species_object, species_string in zip(
        species_order_list, current_species_string_list, strict=False
    ):
        try:
            cyclic_dict[species_object].append(species_string)
        except KeyError:
            cyclic_dict[species_object] = [species_string]

    return cyclic_dict


def construct_product_structure(reaction: Reactions) -> list[dict[str, Any]]:
    """This function unpacks the products in a meta-reaction

    :param: reaction meta-reaction currently being analysed


    Returns:
        A list of dictionaries for each product with the meta-species object, the label
        and the characteristics.
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
            product_list.extend(
                {
                    "species": product["object"],
                    "label": product["label"],
                    "characteristics": product["characteristics"],
                    "stoichiometry": 1,
                }
                for _ in range(product["stoichiometry"])
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

    Args:
        reactant_species_string_list: List of reactants in MobsPy format.
        product_species_string_list: List of products in MobsPy format.
        reaction_rate: Reaction rate expression as a string.


    Returns:
        Dictionary that packs the reactants products and rate.
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

    Args:
        reaction: (meta-reaction object).
        meta_species_in_model: List of meta-species used in the model.


    Returns:
        Base_species_order (list of meta-species objects) = order that the meta-species
        appear in the meta-reaction, reactant_species_combination_list (list of lists of
        meta-species) = list of lists of all meta-species that have inherited from the
        meta-species in the meta-reaction.
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
                    "Please add the species or remove the reaction"
                )

            reactant_species_combination_list.append(species_for_reactant)

    return base_species_order, reactant_species_combination_list


def construct_rate_function_arguments(
    rate_function: Any,
    reaction: Reactions,
) -> list[str]:
    """Extract and validate parameter names from a rate function signature."""
    import warnings  # noqa: PLC0415

    sig = inspect_signature(rate_function)

    for param in sig.parameters.values():
        if isinstance(param.annotation, str):
            warnings.warn(
                f"Rate function '{rate_function.__qualname__}' has "
                "stringified type annotations (possibly from "
                "'from __future__ import annotations'). MobsPy "
                "extracts parameter names only, so this is harmless, "
                "but if you encounter issues, remove the __future__ "
                "import from the file defining your rate functions.",
                stacklevel=4,
            )
            break

    for param in sig.parameters.values():
        if param.kind in (param.VAR_POSITIONAL, param.VAR_KEYWORD):
            raise CompilationError(
                "Rate arguments must not contain * or **. \n"
                f"Error in reaction {reaction}. \n"
                f"Error in rate function {rate_function} "
                f"in signature {sig!s}"
            )
        if param.default is not param.empty:
            raise CompilationError(
                "Rate arguments must not have default values. \n"
                f"Error in reaction {reaction}. \n"
                f"Error in rate function {rate_function} "
                f"in signature {sig!s}"
            )

    return list(sig.parameters.keys())


def create_all_reactions(
    reactions: set[Any],
    meta_species_in_model: Any,
    ref_characteristics_to_object: dict[str, Any],
    parameters_in_reaction: Any,
    ctx: CompilationContext,
) -> tuple[ReactionsForSbml, Any]:
    """This function creates all reactions
    Returns the reactions_for_sbml and parameters_for_sbml dictionary
    Those will be used by another module to create the SBML file

    Args:
        reactions: (meta-reaction objects) reactions objects constructed by the
            meta_class module.
        meta_species_in_model: List of meta-species in model.
        ref_characteristics_to_object: Characteristics as keys objects as values.
        type_of_model: Stochastic or deterministic.
        dimension: Model dimension 1D, 2D, 3D, .....


    Returns:
        Dictionary with all reactions that will be added to the sbml model file,
        parameters_for_sbml (dict) = parameters for the sbml model file.
    """
    reactions_for_sbml: ReactionsForSbml = {}

    check_for_invalid_reactions(reactions, ref_characteristics_to_object)

    # Initiate expressions
    with crs_Unit_Context_Setter(model_context=ctx.model_context):
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
                            ConcreteSpeciesId(
                                base=reactant[0].get_name(),
                                characteristics=tuple(reactant[1:]),
                            ).to_sbml_id()
                            if len(reactant) > 1
                            else reactant[0].get_name()
                            for reactant in reactant_string_list
                        ]

                        try:
                            rate_string, parameters_in_reaction = (
                                fr_extract_reaction_rate(
                                    list(combination_of_reactant_species),
                                    reactant_strings,
                                    reaction.rate,
                                    reaction_rate_arguments,
                                    parameters_in_reaction,
                                    ctx,
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
                            reactant_strings, product_string_list, str(rate_string)
                        )

    return reactions_for_sbml, parameters_in_reaction
