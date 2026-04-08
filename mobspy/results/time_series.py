"""
This module implements the class that stores the results from a MobsPy simulation
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pandas as pd

from mobspy.dsl.reactions import Reacting_Species
from mobspy.dsl.species import Species
from mobspy.exceptions import ValidationError

_TUPLE_PAIR_LEN = 2

if TYPE_CHECKING:
    from collections.abc import Iterator

    from mobspy.dsl.mobspy_parameters import (
        Internal_Parameter_Constructor,
    )
    from mobspy.types import TimeSeriesDataDict


class MobsPyTimeSeries:
    """Single time-series result from one simulation run."""

    def __init__(
        self,
        data_dict: TimeSeriesDataDict,
        model_parameters: dict[str, Any] | None = None,
    ) -> None:
        """Creates the MobsPy timeseries object

        Args:
            data_dict: Resulting dictionary from simulation.

        {'data': ...., 'params':....., 'models':.......}
        """
        self.ts_data: dict[str, list[float]] = data_dict["data"]
        self.ts_parameters: dict[str, Any] = data_dict["params"]
        self.ts_models: list[str] = data_dict["models"]
        if model_parameters is None:
            self.ts_model_parameters: dict[str, Any] = {}
        else:
            self.ts_model_parameters = model_parameters


class SimulationResults:
    """Collection of time-series results across multiple simulation runs.

    Supports indexing by run number, species name, or meta-species object.
    """

    def check_parameters_for_deepcopy(self) -> None:
        """Cast non-primitive parameter values to strings for safe deepcopy."""
        for i, d in enumerate(self.ts_parameters):
            for par, val in d.items():
                if not isinstance(
                    val,
                    (
                        int,
                        float,
                        str,
                        bool,
                    ),
                ):
                    # cast to str
                    self.ts_parameters[i][par] = str(val)

    def __init__(
        self,
        list_of_mspy_ts: list[MobsPyTimeSeries],
        model_parameter_objects: (
            dict[str, Internal_Parameter_Constructor] | None
        ) = None,
        fres: bool = False,
    ) -> None:
        """Creates the MobsPy timeseries object

        Args:
            data_dict: Resulting dictionary from simulation.

        {'data': ...., 'params':....., 'models':.......}
        """

        self.ts_data: list[dict[str, Any]] = [dict(x.ts_data) for x in list_of_mspy_ts]
        self.ts_parameters: list[dict[str, Any]] = [
            dict(x.ts_parameters) for x in list_of_mspy_ts
        ]

        # This lines are here to allow the object to be deepcopiable
        for ts_par in self.ts_parameters:
            if ts_par["unit_x"] is not None:
                ts_par["unit_x"] = str(ts_par["unit_x"])
            if ts_par["unit_y"] is not None:
                ts_par["unit_y"] = str(ts_par["unit_y"])

        self.ts_models: list[list[str]] = [list(x.ts_models) for x in list_of_mspy_ts]
        self.ts_model_parameters: list[dict[str, Any]] = [
            dict(x.ts_model_parameters) for x in list_of_mspy_ts
        ]

        self.check_parameters_for_deepcopy()

        if model_parameter_objects is not None:
            need_conversion_dict: set[str] = set()
            for par, par_object in model_parameter_objects.items():
                if par_object._has_units:
                    need_conversion_dict.add(par)

            for i, par_comb in enumerate(self.ts_model_parameters):
                for par_name in par_comb:
                    if par_name not in need_conversion_dict:
                        continue

                    cf = model_parameter_objects[par_name].conversion_factor
                    self.ts_model_parameters[i][par_name] = (
                        self.ts_model_parameters[i][par_name] / cf
                    )

        self.fres: bool = fres

    def to_dict(self) -> dict[str, Any]:
        """

        Returns:
            Data in dict format {'data': ...., 'params':....., 'models':.......}.
        """
        return {
            "data": self.ts_data,
            "params": self.ts_parameters,
            "models": self.ts_models,
        }

    def __len__(self) -> int:
        """
        Length of MobsPy Time Series is equal to the number of runs from simulation
        """
        return len(self.ts_data)

    def __str__(self) -> str:
        tr: str = ""
        for i, (data, params) in enumerate(
            zip(self.ts_data, self.ts_model_parameters, strict=True)
        ):
            if self.ts_model_parameters:
                tr += f"Model Parameters {params} \n"
            tr += f"Time Series {i}: \n"
            tr += f"{data}\n \n"
        return tr

    def add_ts_to_data(self, time_series: dict[str, list[float]]) -> None:
        """
        Add a new time series to the TS data. Used for stochastic plotting
        the average and standard deviation

        Args:
            time_series: Dictionary with species strings as keys and run as value.
        """
        if isinstance(time_series, dict):
            self.ts_data += [time_series]

    _GetItemKey = (
        int
        | str
        | tuple[str | Species | Reacting_Species, int]
        | Species
        | Reacting_Species
    )

    def __getitem__(self, item: _GetItemKey) -> Any:
        """Implement run retrieval using a meta-species object.

        Returns one run if there is only one time-series and
        returns multiple runs if there are multiple time series.
        """
        if isinstance(item, str) and item == "runs":
            raise ValidationError(
                "As of version 2.0.1 MobsPy has changed the data output format. \n"
                "Now data can be accessed through the following syntax: \n"
                "S.results[Meta-Species Object] or "
                "S.results[Meta-Species string name] \n"
                "Both can perform queries"
            )

        series_index = None
        if isinstance(item, tuple):
            item, series_index = self._unpack_tuple_key(item)

        if isinstance(item, int):
            return self.ts_data[item]

        to_return = self._resolve_item(item, series_index)

        if not self.fres:
            return to_return
        return to_return[0]

    def _unpack_tuple_key(
        self, item: tuple[str | Species | Reacting_Species, int]
    ) -> tuple[str | Species | Reacting_Species, int]:
        """Validate and unpack a tuple-based key into (item, series_index)."""
        if len(item) != _TUPLE_PAIR_LEN or not isinstance(item[1], int):
            raise ValidationError(
                "Only len 2 and ints allowed in tuple-based assignments"
            )
        return item[0], item[1]

    def _resolve_item(
        self,
        item: str | Species | Reacting_Species,
        series_index: int | None,
    ) -> list[Any]:
        """Dispatch retrieval by item type."""
        if isinstance(item, str):
            return self._resolve_str_item(item, series_index)
        if isinstance(item, Species):
            return self._resolve_species_item(item, series_index)
        if isinstance(item, Reacting_Species):
            return self._resolve_reacting_species_item(item, series_index)
        return []

    def _resolve_str_item(self, item: str, series_index: int | None) -> list[Any]:
        """Retrieve data for a string key."""
        try:
            return [ts[item] for ts in self.ts_data]
        except KeyError:
            if series_index is None:
                return [
                    self._sum_reacting_species_data(item, i) for i in range(len(self))
                ]
            return [self._sum_reacting_species_data(item, series_index)]

    def _resolve_species_item(
        self, item: Species, series_index: int | None
    ) -> list[Any]:
        """Retrieve data for a Species key."""
        if series_index is None:
            return [ts[item.get_name()] for ts in self.ts_data]
        return [self._sum_reacting_species_data(item.get_name(), series_index)]

    def _resolve_reacting_species_item(
        self, item: Reacting_Species, series_index: int | None
    ) -> list[Any]:
        """Retrieve data for a Reacting_Species key."""
        if series_index is None:
            return [self._sum_reacting_species_data(item, i) for i in range(len(self))]
        return [self._sum_reacting_species_data(item, series_index)]

    def _sum_reacting_species_data(
        self, item: str | Species | Reacting_Species, ts_index: int
    ) -> list[float]:
        """
        Maps meta-species according to characteristics.
        Ex: A.a1 = A.a1.b1 + A.a1.b2 + A.a1.b3

        Args:
            item: Meta-species object or string to be retrieved.

        :ts_index: Index of the time series to perform the sum
        """

        def _sum_element_by_element(l1: list[float], l2: list[float]) -> list[float]:
            rt: list[float] = []
            for e1, e2 in zip(l1, l2, strict=True):
                rt.append(e1 + e2)
            return rt

        time_series = self.ts_data[ts_index]
        to_return: list[float] = [0.0 for _ in range(len(time_series["Time"]))]

        found_flag = False
        if isinstance(item, Reacting_Species):
            for reactant in item.list_of_reactants:
                for key in time_series:
                    reactant_name = reactant["object"].get_name()
                    reactant_full_name = reactant["characteristics"].union(
                        [reactant_name]
                    )
                    if reactant_full_name.issubset(set(key.split("."))):
                        to_return = _sum_element_by_element(to_return, time_series[key])
                        found_flag = True
        elif isinstance(item, str):
            for key in time_series:
                if set(item.split(".")).issubset(set(key.split("."))):
                    to_return = _sum_element_by_element(to_return, time_series[key])
                    found_flag = True

        if not found_flag:
            raise ValidationError(f"{item} was not found in data")

        return to_return

    def __iter__(self) -> Iterator[Any]:
        if not self.fres:
            yield from self.ts_data
        else:
            yield from self.ts_data[0]

    def get_max_time_for_species(self, species: str | Species) -> list[float]:
        """
        Returns the maximum time in all the time-series stored for a given species
        """
        if not isinstance(species, str):
            species = species.get_name()
        max_length: int = 0
        max_ts: dict[str, Any] | None = None
        for ts in self.ts_data:
            if species in ts and len(ts["Time"]) > max_length:
                max_length = len(ts["Time"])
                max_ts = ts
        if max_ts is None:
            msg = "Could not find maximal time series."
            raise ValueError(msg)
        return max_ts["Time"]  # type: ignore[no-any-return]  # dynamic dispatch

    def return_pandas(self) -> list[pd.DataFrame]:
        """Convert each run's time-series data to a pandas DataFrame."""
        return [pd.DataFrame.from_dict(ts) for ts in self.ts_data]
