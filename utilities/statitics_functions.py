import numpy
import logging

logging.basicConfig(level=logging.WARNING)
module_logging = logging.getLogger(__name__)
logging_array = [module_logging]


def time_series_average(list_series):

    '''
    :param list_series: a list of all time series to calculate the average
    :return: a single series with average value for each time
    '''

    average_series = []
    for j in range(len(list_series[0])):

        add = 0
        try:
            for series in list_series:
                add = add + series[j]

        except KeyError:
            logging_array[0].error('The time_series_average function was programed for same size time series')
            exit()

        average_series.append(add/len(list_series))

    return average_series


def standard_deviation(list_series):

    '''
    :param list_series: a list of all time series to calculate the standard deviation
    :return: the standard deviation per time value
    '''

    average_series = time_series_average(list_series)
    deviation_series = []

    for j in range(len(list_series[0])):

        add = 0
        try:
            for series in list_series:
                add = add + (average_series[j] - series[j])**2

        except KeyError:
            logging_array[0].error('The time_series_average function was programed for same size time series')
            exit()

        deviation_series.append(add/len(list_series))

    return deviation_series
