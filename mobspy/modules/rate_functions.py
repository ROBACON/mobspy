"""Rate assignment and mass-action kinetics construction for reactions.

Handles rate_assignment to the reactions:
    The passing of arguments to be executed by the user
    supplied rate function. The construction of mass action
    kinetics. Handles the different types of rates supplied
    by the user.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Any

from pint import Quantity

from mobspy.constants import NULL_SPECIES
from mobspy.exceptions import CompilationError
from mobspy.mobspy_logging import get_logger
from mobspy.modules.meta_class import Zero as mc_Zero
from mobspy.modules.mobspy_expressions import (
    ExpressionDefiner as mbe_ExpressionDefiner,
)
from mobspy.modules.mobspy_expressions import (
    MobsPyExpression as mbe_MobsPyExpression,
)
from mobspy.modules.mobspy_expressions import (
    Specific_Species_Operator as mbe_Specific_Species_Operator,
)
from mobspy.modules.mobspy_parameters import (
    Internal_Parameter_Constructor as mp_Mobspy_Parameter,
)
from mobspy.modules.species_utils import (
    count_stoichiometry as mcu_count_string_dictionary,
)
from mobspy.modules.unit_handler import convert_rate as uh_convert_rate
from mobspy.types import CompilationContext

if TYPE_CHECKING:
    from mobspy.modules.model_unit_context import ModelUnitContext

_logger = get_logger(__name__)

_RateFunction = (
    int | float | Quantity | str | Callable[..., Any] | mbe_ExpressionDefiner | None
)


def extract_reaction_rate(  # noqa: PLR0913
    combination_of_reactant_species: list[dict[str, Any]],
    reactant_string_list: list[str],
    reaction_rate_function: _RateFunction,
    function_rate_arguments: list[str] | None,
    parameters_in_reaction: set[mp_Mobspy_Parameter],
    ctx: CompilationContext,
) -> tuple[str | int, set[mp_Mobspy_Parameter]]:
    """Return the reaction rate string for model construction.

    Does a different action depending on the type of the
    reaction_rate_function (we consider constants as
    functions). It passes the rate as a string expression
    in MobsPy standard units.

    Args:
        function_rate_arguments: List of strings of the function rate argument ex:['r1',
            'r2', ....].
        combination_of_reactant_species: Meta-species currently being used in this
            reaction.
        reactant_string_list: List of species strings in order they appear in the
            reaction.
        reaction_rate_function: Rate stored in the reaction object.
        dimension: System dimension (for the rate conversion).


    Returns:
        The reaction kinetics as a string for SBML.
    """
    dimension = ctx.dimension
    type_of_model = ctx.type_of_model
    model_context = ctx.model_context
    parameter_exist = ctx.parameter_exist

    is_count = False
    if isinstance(reaction_rate_function, (int, float, Quantity)):
        reaction_rate_function, dimension, is_count = uh_convert_rate(
            reaction_rate_function,
            len(reactant_string_list),
            dimension,
            model_context=model_context,
        )

        if reaction_rate_function == 0:
            return 0, parameters_in_reaction
        reaction_rate_string = basic_kinetics_string(
            reactant_string_list,
            reaction_rate_function,  # pyright: ignore[reportArgumentType]
            type_of_model,
            is_count,
        )

    elif isinstance(reaction_rate_function, mbe_ExpressionDefiner) and parameter_exist:
        parameters_in_reaction = parameters_in_reaction.union(
            reaction_rate_function._parameter_set
        )
        reaction_rate_string = basic_kinetics_string(
            reactant_string_list,
            reaction_rate_function,
            type_of_model,
            is_count,
        )

    elif function_rate_arguments is not None:
        return _process_callable_rate(
            combination_of_reactant_species,
            reactant_string_list,
            reaction_rate_function,
            function_rate_arguments,
            parameters_in_reaction,
            ctx,
        )
    elif reaction_rate_function is None:
        raise CompilationError(
            "There is a reaction rate missing for the "
            "following reactants: \n" + str(reactant_string_list)
        )
    elif isinstance(reaction_rate_function, str):
        reaction_rate_string = reaction_rate_function
    else:
        _logger.debug(str(type(reaction_rate_function)))
        raise CompilationError(
            f"The {type(reaction_rate_function)},"
            f" from {reaction_rate_function} is not supported"
        )

    return reaction_rate_string, parameters_in_reaction  # type: ignore[possibly-undefined]


def _process_callable_rate(  # noqa: PLR0913
    combination_of_reactant_species: list[dict[str, Any]],
    reactant_string_list: list[str],
    reaction_rate_function: _RateFunction,
    function_rate_arguments: list[str],
    parameters_in_reaction: set[mp_Mobspy_Parameter],
    ctx: CompilationContext,
) -> tuple[str | int, set[mp_Mobspy_Parameter]]:
    """Process a callable rate function and return (rate_string, parameters)."""
    dimension = ctx.dimension
    type_of_model = ctx.type_of_model
    model_context = ctx.model_context
    skip_check = ctx.skip_expression_check

    if function_rate_arguments != [""]:
        arguments = prepare_arguments_for_callable(
            combination_of_reactant_species,
            reactant_string_list,
            function_rate_arguments,
            dimension,
            model_context=model_context,
        )
        rate = reaction_rate_function(**arguments)  # type: ignore[operator, misc]
    else:
        rate = reaction_rate_function()  # type: ignore[operator, misc]

    rate, dimension, is_count = uh_convert_rate(
        rate,
        len(reactant_string_list),
        dimension,
        model_context=model_context,
    )

    if rate == 0:
        return 0, parameters_in_reaction

    if isinstance(rate, (int, float)):
        reaction_rate_string = basic_kinetics_string(
            reactant_string_list,
            rate,
            type_of_model,
            is_count,
        )
    elif isinstance(rate, str):
        reaction_rate_string = rate.replace("$", "")
    elif isinstance(rate, mbe_MobsPyExpression):
        reaction_rate_string, parameters_in_reaction = _process_expression_rate(
            rate,
            reactant_string_list,
            parameters_in_reaction,
            type_of_model,
            dimension,
            skip_check,
        )
    elif isinstance(rate, mp_Mobspy_Parameter):
        parameters_in_reaction.add(rate)
        reaction_rate_string = basic_kinetics_string(
            reactant_string_list,
            str(rate),
            type_of_model,
        )
    elif rate is None:
        raise CompilationError(
            "There is a reaction rate missing for the "
            "following reactants: \n" + str(reactant_string_list)
        )
    else:
        raise CompilationError(
            f"The rate function {reaction_rate_function},"
            " returned a non-valid value. \n"
            "Only int, floats and str are accepted"
        )

    return reaction_rate_string, parameters_in_reaction


def _process_expression_rate(  # noqa: PLR0913
    rate: mbe_MobsPyExpression,
    reactant_string_list: list[str],
    parameters_in_reaction: set[mp_Mobspy_Parameter],
    type_of_model: str,
    dimension: int | None,
    skip_check: bool,
) -> tuple[str, set[mp_Mobspy_Parameter]]:
    """Convert a MobsPyExpression rate to a string."""
    if len(rate._expression_variables) > 0:
        reaction_rate_string, _ = rate.generate_string_operation(
            skip_check=skip_check,
            dimension=dimension,
        )
        parameters_in_reaction = parameters_in_reaction.union(rate._parameter_set)
    else:
        rate_for_mass_action, is_count = rate.generate_string_operation(
            reaction_order=len(reactant_string_list),
            dimension=dimension,
            skip_check=skip_check,
        )
        parameters_in_reaction = parameters_in_reaction.union(rate._parameter_set)
        reaction_rate_string = basic_kinetics_string(
            reactant_string_list,
            rate_for_mass_action,
            type_of_model,
            is_count,
        )
    return reaction_rate_string, parameters_in_reaction


def basic_kinetics_string(
    reactants: list[str],
    reaction_rate: int | float | str | mbe_ExpressionDefiner,
    type_of_model: str,
    is_count: bool = False,
) -> str:
    """Construct mass-action kinetics string.

    Both for stochastic and deterministic depending on the
    type of model.

    Args:
        reactants: List of reactants in MobsPy str format.
        reaction_rate: Reaction constant.
        type_of_model: Stochastic or deterministic, rate
            expressions differ depending on each case.

    Returns:
        Mass action kinetics expression for the reaction.
    """
    counts = mcu_count_string_dictionary(reactants)

    kinetics_string = ""
    for name, number in counts.items():
        if type_of_model.lower() == "stochastic":
            kinetics_string += stochastic_string(name, int(number))
        elif type_of_model.lower() == "deterministic":
            kinetics_string += deterministic_string(name, int(number))
        kinetics_string += " * "

    kinetics_string += str(reaction_rate)

    n: int = 0
    for item in counts.values():
        n += int(item)
    n = n - 1

    if n > 0 and not is_count:
        kinetics_string += f" * volume^{-n}"
    elif n < 0 and not is_count:
        kinetics_string += " * volume"
    else:
        pass

    return kinetics_string


def stochastic_string(reactant_name: str, number: int) -> str:
    """Return stochastic string for mass action kinetics.

    For instance the reaction 2A -> 3A would imply
    A*(A-1)/2. It only does so for one reactant, so it
    must be called for all reactants in the reaction.

    Args:
        reactant_name: Species string involved in the reaction.
        number: Stoichiometry (number of times it appears).

    Returns:
        The mass action kinetics string expression for only that species.
    """
    to_return_string: str = ""
    for i in range(number):
        if i == 0:
            to_return_string = reactant_name
        else:
            to_return_string += f" * ({reactant_name} - {i})/{i + 1}"

    return to_return_string


def deterministic_string(reactant_name: str, number: int) -> str:
    """Return deterministic string for mass action kinetics.

    For instance the reaction 2A -> 3A would imply A*A.
    It only does so for one reactant, so it must be called
    for all reactants in the reaction.

    Args:
        reactant_name: Species string involved in the reaction.
        number: Stoichiometry (number of times it appears).

    Returns:
        The mass action kinetics string expression for only that species.
    """
    to_return_string = ""
    for i in range(number):
        if i == 0:
            to_return_string = reactant_name
        else:
            to_return_string += f" * {reactant_name}"
    return to_return_string


def prepare_arguments_for_callable(
    combination_of_reactant_species: list[dict[str, Any]],
    reactant_string_list: list[str],
    rate_function_arguments: list[str],
    dimension: int | None,
    model_context: ModelUnitContext | None = None,  # noqa: ARG001  # reserved for future unit support
) -> dict[str, mbe_MobsPyExpression | mbe_Specific_Species_Operator]:
    """Prepare arguments for the rate function.

    Creates objects of the Specific_Species_Operator class
    for a given reaction.

    Args:
        combination_of_reactant_species: Meta-species involved
            in the reaction.
        reactant_string_list: Species strings involved in the
            reaction.
        rate_function_arguments: Arguments received by the rate
            function.

    Returns:
        Dictionary with arguments for a rate function.
    """
    argument_dict: dict[
        str,
        mbe_MobsPyExpression | mbe_Specific_Species_Operator,
    ] = {}

    if rate_function_arguments is not None:
        i = 0
        for i, (species, reactant_string) in enumerate(
            zip(
                combination_of_reactant_species,
                reactant_string_list,
                strict=False,
            )
        ):
            try:
                species = species["object"]  # noqa: PLW2901
                argument_dict[rate_function_arguments[i]] = mbe_MobsPyExpression(
                    reactant_string,
                    species,
                    dimension=dimension,
                    count_in_model=True,
                    concentration_in_model=False,
                    count_in_expression=False,
                    concentration_in_expression=False,
                )
            except IndexError:
                continue

        if i != 0:
            while len(argument_dict) < len(rate_function_arguments):
                i += 1
                argument_dict[rate_function_arguments[i]] = (
                    mbe_Specific_Species_Operator(NULL_SPECIES, mc_Zero)
                )
        elif i == 0:
            while len(argument_dict) < len(rate_function_arguments):
                argument_dict[rate_function_arguments[i]] = (
                    mbe_Specific_Species_Operator(NULL_SPECIES, mc_Zero)
                )
                i += 1

    return argument_dict
