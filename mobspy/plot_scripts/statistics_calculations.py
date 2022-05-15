import numpy as np


def time_series_average(list_series):
    """
        Badly named function - Average between all RUNS inside a single time-series

        Parameters:
            list_series (list) = a list of all runs from a time-series to calculate the average

        Returns:
            average_series (list) = a list with the average values from all runs
    """
    average_series = []
    for j in range(len(list_series[0])):

        add = 0
        try:
            for series in list_series:
                add = add + series[j]

        except KeyError:
            raise ValueError('The time_series_average function was programed for same size time series')

        average_series.append(add/len(list_series))

    return average_series


def standard_deviation(list_series):
    """
       Standard deviation between all RUNS inside a single time-series

       Parameters:
           list_series (list) = a list of all runs from a time-series to calculate the average

       Returns:
           deviation_series (list) = a list with the standard deviation values from all runs
   """
    average_series = time_series_average(list_series)
    deviation_series = []

    for j in range(len(list_series[0])):

        add = 0
        try:
            for series in list_series:
                add = add + (average_series[j] - series[j])**2

        except KeyError:
            raise ValueError('The time_series_average function was programed for same size time series')

        deviation_series.append(np.sqrt(add/len(list_series)))

    return deviation_series


def average_plus_standard_deviation(list_series):
    """
        Combines the average with the standard deviation

        Parameters:
           list_series (list) = a list of all runs from a time-series to calculate the average

       Returns:
           series_average (list) = a list with the average values from all runs
           plus (list) = average list + std deviation list (element by element)
           minus (list) = average list - std deviation list (element by element)
    """
    series_average = time_series_average(list_series)
    series_deviation = standard_deviation(list_series)

    plus = []
    minus = []
    for average, deviation in zip(series_average, series_deviation):
        plus.append(average + deviation)
        minus.append(average - deviation)

    return series_average, plus, minus
