"""Compute averages and standard deviations across stochastic simulation runs."""

from __future__ import annotations

from typing import Any

import numpy as np

from mobspy.mobspy_logging import get_logger

_logger = get_logger(__name__)


def time_series_average(species_string: str, mobspy_ts: Any) -> list[float]:
    """
    Badly named function - Average between all RUNS inside a single time-series

    Args:
        species_string: Species string to perform the average upon.
        mobspy_ts: MobsPy time series object.


    Returns:
        A list with the average values from all runs.
    """

    list_series: list[Any] = [
        series for series in mobspy_ts if species_string in series
    ]

    war_1, war_2 = (True, True)
    for s1 in list_series:
        for s2 in list_series:
            if s1 == s2:
                continue

            if len(s1) != len(s2) and war_1:
                _logger.warning(
                    "Time Series length is different. \n"
                    "MobsPy disregards time-series that have"
                    " already finished during calculations"
                )
                war_1 = False

            for t1, t2 in zip(s1["Time"], s2["Time"], strict=False):
                if t1 != t2 and war_2:
                    _logger.warning(
                        "Times in Time Series Objects are different. \n"
                        "MobsPy calculates the average by index"
                        " position. Please be careful."
                    )
                war_2 = False

    average_series: list[float] = []
    for j in range(len(mobspy_ts.get_max_time_for_species(species_string))):
        add: float = 0
        size = 0
        for series in list_series:
            try:
                add = add + series[species_string][j]
                size += 1
            except IndexError:
                pass
            except KeyError:
                pass

        average_series.append(add / size)

    return average_series


def standard_deviation(
    species_string: str, mobspy_ts: Any, average_series: list[float] | None = None
) -> list[Any]:
    """
    Standard deviation between all RUNS inside a single time-series

    Args:
        species_string: Species string.
        mobspy_ts: MobsPy time series object.
        average_series: If the average series is given it is not recalculated.


    Returns:
        A list with the standard deviation values from all runs.
    """
    if average_series is None:
        average_series = time_series_average(species_string, mobspy_ts)
    deviation_series: list[Any] = []

    for j in range(len(mobspy_ts.get_max_time_for_species(species_string))):
        add: float = 0
        size = 0
        for series in mobspy_ts:
            try:
                add = add + (average_series[j] - series[species_string][j]) ** 2
                size += 1
            except IndexError:
                pass
            except KeyError:
                pass

        deviation_series.append(np.sqrt(add / size))

    return deviation_series


def average_plus_standard_deviation(
    species_string: str,
    mobspy_ts: Any,
    average_series: list[float] | None = None,
    deviation_series: list[Any] | None = None,
) -> tuple[list[float], list[float], list[float]]:
    """
    Standard deviation between all RUNS inside a single time-series

    Args:
        species_string: Species string.
        mobspy_ts: MobsPy time series object.
        average_series: If the average series is given it is not recalculated.
        deviation_series: If the deviation is given it is not recalculated.


    Returns:
        Average value of the run, plus (list) = average + deviation, minus (list) =
        average - deviation.
    """
    series_average = (
        time_series_average(species_string, mobspy_ts)
        if average_series is None
        else average_series
    )
    series_deviation = (
        standard_deviation(species_string, mobspy_ts)
        if deviation_series is None
        else deviation_series
    )

    plus: list[float] = []
    minus: list[float] = []
    for average, deviation in zip(series_average, series_deviation, strict=False):
        plus.append(average + deviation)
        minus.append(average - deviation)

    return series_average, plus, minus
