"""Run COPASI parameter estimation tasks and collect fitted parameter results."""

from __future__ import annotations

from typing import Any

from pandas import DataFrame

from mobspy.exceptions import ParameterError
from mobspy.import_manager.lazy_import_class import (
    LazyImporter as ipm_LazyImporter,
)
from mobspy.mobspy_logging import get_logger

_logger = get_logger(__name__)

basico = ipm_LazyImporter("basico")


def basiCO_parameter_estimation(
    simulation_object: Any,
    parameters_to_estimate: list[Any] | set[Any] | tuple[Any, ...],
    experimental_data: Any = None,
    bound: dict[str, list[float]] | list[float] | tuple[float, ...] | None = None,
    method: str = "Evolution Strategy (SRES)",
    verbose: bool = True,
    change_parameter_values: bool = True,
) -> dict[str, Any]:
    """
    This function fits a MobsPy parameter in a model to experimental data.
    Bound specifies the search range.
    If the bound is None it takes as bounds /1000 and *1000
    the original value of the parameter
    The function updates the value of the parameter when it is done

    Args:
        simulation_object: - MobsPy Simulation Object.
        parameters_to_estimate: List of MobsPy parameters to estimate.
        experimental_data: - experimental data to fit the parameter with.
        bound: - list of two elements with a lower and upper bound of the parameters
            values, or dictionary with the parameter name and a two element for that
            specific parameter.
        method: - method to be used by basiCO optimisation - The options are Random
            Search, Simulated Annealing, Differential Evolution, Scatter Search, Genetic
            Algorithm, Evolutionary Programming, Genetic Algorithm SR, Evolution
            Strategy (SRES), Particle Swarm.
        verbose: Print the results after finishing or not.
        change_parameter_values: Change/convert parameter values when possible.
    """

    # Check inputs in order ############################
    # Parameters ok?
    if not isinstance(parameters_to_estimate, (list, set, tuple)):
        raise ParameterError(
            "The parameter that will be estimated must be inside a list, set or tuple"
        )

    converted_parameters: list[str] = []
    # If bound is None, we set the auto-bound when checking the parameters
    if bound is None:
        bound = {}
        flag_auto_set = True
    else:
        flag_auto_set = False

    # Conversion of MobsPy parameters to str, construction of bound
    original_parameters = parameters_to_estimate
    for par in parameters_to_estimate:
        if flag_auto_set:
            if not isinstance(bound, dict):
                raise ParameterError(
                    "bound must be a dict when auto-setting parameter bounds"
                )
            bound[str(par)] = [par.value / 1000, par.value * 1000]

        converted_parameters.append(str(par))
    parameters_to_estimate = converted_parameters

    # Simulation Object ok?
    if simulation_object.experimental_data is not None:
        experimental_data = simulation_object.experimental_data.return_pandas()[0]
    if experimental_data is None:
        raise ParameterError(
            "No experimental data found in the simulation object or as argument of the "
            "basiCO_parameter_estimation function"
        )

    # Experimental data ok?
    if not isinstance(experimental_data, (list, set, tuple, DataFrame)):
        raise ParameterError(
            "Experimental for basiCO estimation must be a "
            "list of pandas dataframes or a pandas dataframe"
        )

    # Bound ok?
    new_bound: dict[str, Any] = {}
    if isinstance(bound, dict):
        for key in bound:
            new_bound[str(key)] = bound[key]
        bound = new_bound

        for par in parameters_to_estimate:
            if par not in bound:
                raise ParameterError(
                    "If a dictionary is used for the bounds,"
                    " all parameters range for estimation "
                    "must be specified. Make sure the "
                    "dictionary keys contains all "
                    "parameters and a list with upper and "
                    "lower bound value for the each "
                    "parameter is given as the items "
                )
    else:
        try:
            if len(bound) != 2:
                raise ParameterError(
                    "The bound argument must be a list with "
                    "the lower and upper bound of all "
                    "parameters"
                )
        except TypeError as e:
            raise ParameterError(
                "The bound argument must be a list with "
                "the lower and upper bound of all "
                "parameters"
            ) from e
    ########################################################################################################

    sbml_list = simulation_object.generate_sbml()
    if len(sbml_list) > 1:
        raise ParameterError(
            "BasiCO optimization does not support composite simulation optimization"
        )

    sbml_str = sbml_list[0]
    model = basico.model_io.load_model_from_string(sbml_str)

    fit_list: list[dict[str, Any]] = []
    for par in parameters_to_estimate:
        basico_reaction_dict = find_parameters_in_basico_dataframe(
            basico.get_reaction_parameters(), par
        ).to_dict()
        try:
            basico_parameter_name = next(iter(basico_reaction_dict["reaction"].keys()))
        except IndexError as e:
            raise ParameterError(
                f"Parameter {par} was not found in the "
                "Simulation model. \n "
                "Please make sure that any of the "
                "meta-species used to construct the "
                "simulator use the parameter in one "
                "of their reactions."
            ) from e

        if bound is None:
            pass

        if isinstance(bound, dict):
            fit_dictionary = {
                "name": basico_parameter_name,
                "lower": bound[par][0],
                "upper": bound[par][1],
            }
        else:
            if not isinstance(bound, (list, tuple)):
                raise ParameterError(
                    f"bound must be a list or tuple, got {type(bound).__name__}"
                )
            fit_dictionary = {
                "name": basico_parameter_name,
                "lower": bound[0],
                "upper": bound[1],
            }
        fit_list.append(fit_dictionary)

    if isinstance(experimental_data, (list, set, tuple)):
        for i, exp in enumerate(experimental_data):
            basico.add_experiment("exp" + str(i), exp, model=model)
    else:
        basico.add_experiment("exp1", experimental_data, model=model)

    basico.set_fit_parameters(fit_list, model=model)
    basico_results = basico.run_parameter_estimation(model=model, method=method)

    # WHY IS EVERYTHING PANDAS IN BASICO??????????? - I don't get paid enough for this
    results: dict[str, Any] = {}
    for key, sol in basico_results.to_dict()["sol"].items():
        parameter_name = str(key).replace("Values[", "")
        parameter_name = str(parameter_name).replace("]", "")
        results[parameter_name] = sol

    if change_parameter_values:
        for p in original_parameters:
            if p.has_units():
                p.set_value(results[str(p)]).convert_to_original_unit()
                results[str(p)] = p.value
            else:
                p.set_value(results[str(p)])

    if verbose:
        _logger.info("Parameter estimation complete. The results follow: ")
        for key in results:
            _logger.info(f"{key} {results[key]}")

    return results


def find_parameters_in_basico_dataframe(
    basico_reactions_df: Any, mobspy_parameter_name: str
) -> Any:
    """
    Finds the corresponding mobspy parameter in the basico
    reactions parameters dataframe.

    Args:
        basico_reactions_df: BasiCO reactions df.
        mobspy_parameter_name: Name of the mobspy parameter.
    """

    def has_mobspy_parameter_in_name(
        basico_reaction_name: str, mobspy_parameter_name: str
    ) -> bool:
        """Check if a BasiCO reaction name contains the parameter."""
        parameter_name = basico_reaction_name.split(".")[1]
        return True if parameter_name == mobspy_parameter_name else False  # noqa: SIM210

    def find_common_substrings(df: str) -> bool:
        """Filter predicate for matching the parameter name."""
        return has_mobspy_parameter_in_name(df, mobspy_parameter_name)

    return basico_reactions_df[
        basico_reactions_df.index.to_frame()["name"].apply(find_common_substrings)
    ]


def python_parameter_estimation() -> None:
    """Run parameter estimation (not yet implemented)."""
