"""
Plotting module: plot generation from MobsPy simulation results.

Provides the PlottingMixin class that adds plotting capabilities
(deterministic, stochastic, parametric, raw) to the Simulation class.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from mobspy.exceptions import SimulationError, ValidationError
from mobspy.modules.reactions import Reacting_Species
from mobspy.modules.species import Species
from mobspy.plot_scripts.default_plots import (
    deterministic_plot as dp_deterministic_plot,
)
from mobspy.plot_scripts.default_plots import (
    parametric_plot as dp_parametric_plot,
)
from mobspy.plot_scripts.default_plots import (
    raw_plot as dp_raw_plot,
)
from mobspy.plot_scripts.default_plots import (
    stochastic_plot as dp_stochastic_plot,
)

if TYPE_CHECKING:
    from mobspy.types import CompiledModelDict


if TYPE_CHECKING:
    from mobspy.data_handler.time_series_object import SimulationResults


def plot_results(
    results: Any,
    plot_parameters: dict[str, Any],
    species_strings: set[str] | None = None,
    method: str = "deterministic",
) -> Any:
    """Plot simulation results without a Simulation instance.

    Standalone function for composition-friendly workflows.

    Args:
        results: SimulationResults or similar result object.
        plot_parameters: Plot configuration dictionary.
        species_strings: Species to plot (None for all).
        method: ``"deterministic"`` or ``"stochastic"``.

    Returns:
        The plot object.
    """
    strings = species_strings or set()
    if method == "stochastic":
        return dp_stochastic_plot(strings, results, plot_parameters)
    return dp_deterministic_plot(strings, results, plot_parameters)


class PlottingMixin:
    """Mixin providing plotting capabilities for Simulation."""

    # These attributes are provided by the Simulation class
    _list_of_models: list[CompiledModelDict]
    results: SimulationResults | dict[str, Any]
    plot_parameters: dict[str, Any]

    def extract_plot_essentials(
        self, *species: str | Species | Reacting_Species
    ) -> tuple[set[str], Any, dict[str, Any]]:
        """Extract essential information for plotting.

        Args:
            *species: Meta-species objects or strings to plot.

        Returns:
            Tuple of (species_strings, results, plot_parameters).
        """
        if not species:
            species_strings: set[str] = set()
            for model in self._list_of_models:
                species_strings = species_strings.union(model.mappings)
        else:
            species_strings = set()

        for spe in species:
            if isinstance(spe, (Species, Reacting_Species)):
                species_strings.add(str(spe))
            elif isinstance(spe, str):
                species_strings.add(spe)
            else:
                raise ValidationError(
                    "Only species objects or strings are accepted as plotting arguments"
                )

        return species_strings, self.results, self.plot_parameters

    def plot_stochastic(self, *species: str | Species | Reacting_Species) -> Any:
        """Generate a stochastic plot of the simulation results.

        Args:
            *species: Variable number of species to plot.

        Raises:
            SimulationError: If no results are available for plotting.
        """
        if not hasattr(self, "results") or not self.results:
            raise SimulationError(
                "No simulation results available for plotting. "
                "Call .run() on the Simulation object first."
            )

        spe_strings, results, params = self.extract_plot_essentials(*species)
        return plot_results(results, params, spe_strings, method="stochastic")

    def plot_deterministic(self, *species: str | Species | Reacting_Species) -> Any:
        """Generate a deterministic plot of the simulation results.

        Args:
            *species: Variable number of species to plot.

        Raises:
            SimulationError: If no results are available for plotting.
        """
        if not hasattr(self, "results") or not self.results:
            raise SimulationError(
                "No simulation results available for plotting. "
                "Call .run() on the Simulation object first."
            )

        spe_strings, results, params = self.extract_plot_essentials(*species)
        return plot_results(results, params, spe_strings, method="deterministic")

    def plot_parametric(self, *species: str | Species | Reacting_Species) -> Any:
        """Generate a parametric plot of the simulation results.

        Args:
            *species: Variable number of species to plot.

        Raises:
            SimulationError: If no results are available for plotting.
        """
        if not hasattr(self, "results") or not self.results:
            raise SimulationError(
                "No simulation results available for plotting. "
                "Call .run() on the Simulation object first."
            )

        spe_strings, results, params = self.extract_plot_essentials(*species)
        return dp_parametric_plot(spe_strings, results, params)

    def plot(self, *species: str | Species | Reacting_Species) -> Any:
        """Plot simulation results, auto-detecting the method.

        Dispatches to ``plot_stochastic`` or ``plot_deterministic``
        based on the simulation method used.
        """
        method = self.plot_parameters.get("simulation_method", "deterministic")
        if method == "stochastic":
            return self.plot_stochastic(*species)
        return self.plot_deterministic(*species)

    def plot_raw(
        self, parameters_or_file: str | dict[str, Any], return_fig: bool = False
    ) -> Any:
        """Generate a raw plot with custom parameters.

        Args:
            parameters_or_file: JSON file name or dict with plot config.
            return_fig: If True, return the figure object instead of displaying it.

        Raises:
            SimulationError: If no results are available for plotting.
        """
        if not hasattr(self, "results") or not self.results:
            raise SimulationError(
                "No simulation results available for plotting. "
                "Call .run() on the Simulation object first."
            )

        return dp_raw_plot(self.results, parameters_or_file, return_fig=return_fig)
