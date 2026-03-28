"""
This module implements the class that stores the results from a MobsPy simulation
"""

from __future__ import annotations

import inspect
from typing import TYPE_CHECKING, Any

import pandas as pd

from mobspy.exceptions import ValidationError
from mobspy.modules.meta_class import Reacting_Species, Species

if TYPE_CHECKING:
    from collections.abc import Iterator

    from mobspy.modules.mobspy_parameters import (
        Internal_Parameter_Constructor,
    )
    from mobspy.types import TimeSeriesDataDict



class MobsPyTimeSeries:
    def __init__(
        self,
        data_dict: TimeSeriesDataDict,
        model_parameters: dict[str, Any] | None = None,
    ) -> None:
        """Creates the MobsPy timeseries object

        :param data_dict: (dict) resulting dictionary from simulation
        {'data': ...., 'params':....., 'models':.......}
        """
        self.ts_data: dict[str, list[float]] = data_dict["data"]
        self.ts_parameters: dict[str, Any] = data_dict["params"]
        self.ts_models: list[str] = data_dict["models"]
        if model_parameters is None:
            self.ts_model_parameters: dict[str, Any] = {}
        else:
            self.ts_model_parameters = model_parameters


class MobsPyList_of_TS:
    def check_parameters_for_deepcopy(self) -> None:
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

        :param data_dict: (dict) resulting dictionary from simulation
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
                if par_object._has_units == "T":
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
        :return: data in dict format {'data': ...., 'params':....., 'models':.......}
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
            zip(self.ts_data, self.ts_model_parameters, strict=False)
        ):
            if self.ts_model_parameters != {}:
                tr += f"Model Parameters {params} \n"
            tr += f"Time Series {i}: \n"
            tr += f"{data}\n \n"
        return tr

    def add_ts_to_data(self, time_series: dict[str, list[float]]) -> None:
        """
        Add a new time series to the TS data. Used for stochastic plotting
        the average and standard deviation

        :param time_series: (dict) dictionary with species strings
            as keys and run as value
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
        """
        Implements run retrieval using a meta-species object.
        Returns one run if there is only one time-series and
        returns multiple runs if there are multiple time series
        """
        code_line = inspect.stack()[1].code_context[0][:-1]
        if "['runs']" in code_line:
            raise ValidationError(
                "As of version 2.0.1 MobsPy has changed the data output format. \n"
                "Now data can be accessed through the following syntax: \n"
                "S.results[Meta-Species Object] or "
                "S.results[Meta-Species string name] \n"
                "Both can perform queries")

        series_index = None
        if isinstance(item, tuple):
            if len(item) != 2 or not isinstance(item[1], int):
                raise ValidationError(
                    "Only len 2 and ints allowed in tuple-based assignments")
            series_index = item[1]
            item = item[0]

        to_return = []
        if isinstance(item, int):
            return self.ts_data[item]
        elif isinstance(item, str):
            try:
                for ts in self.ts_data:
                    to_return.append(ts[item])
            except KeyError:
                if series_index is None:
                    for i in range(len(self)):
                        to_return.append(self._sum_reacting_species_data(item, i))
                else:
                    to_return.append(
                        self._sum_reacting_species_data(item, series_index)
                    )
        elif isinstance(item, Species):
            if series_index is None:
                for ts in self.ts_data:
                    to_return.append(ts[item.get_name()])
            else:
                to_return.append(self._sum_reacting_species_data(item, series_index))
        elif isinstance(item, Reacting_Species):
            if series_index is None:
                for i in range(len(self)):
                    to_return.append(self._sum_reacting_species_data(item, i))
            else:
                to_return.append(self._sum_reacting_species_data(item, series_index))

        if not self.fres:
            return to_return
        else:
            return to_return[0]

    def _sum_reacting_species_data(
        self, item: str | Reacting_Species, ts_index: int
    ) -> list[float]:
        """
        Maps meta-species according to characteristics.
        Ex: A.a1 = A.a1.b1 + A.a1.b2 + A.a1.b3

        :param item: Meta-species object or string to be retrieved
        :ts_index: Index of the time series to perform the sum
        """

        def _sum_element_by_element(l1: list[float], l2: list[float]) -> list[float]:
            rt: list[float] = []
            for e1, e2 in zip(l1, l2, strict=False):
                rt.append(e1 + e2)
            return rt

        time_series = self.ts_data[ts_index]
        to_return = [0 for _ in range(len(time_series["Time"]))]

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
            if species in ts and len(ts) > max_length:
                max_length = len(ts)
                max_ts = ts
        if max_ts is None:
            msg = "Could not find maximal time series."
            raise ValueError(msg)
        return max_ts["Time"]

    def return_pandas(self) -> list[pd.DataFrame]:
        to_return: list[pd.DataFrame] = []

        for ts in self.ts_data:
            to_return.append(pd.DataFrame.from_dict(ts))

        return to_return
