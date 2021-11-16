import logging
import os, inspect, sys
from pathlib import Path


# TODO: Rework logging? It appears that it sucks
# TODO: Fucking hate this module
class ProgramLogger:

    '''
        This class is responsible for keeping track of the logs through the simulation
        It writes and saves logs
        It also prints the logs to the screen
        It is a simple encapsulation of logging to better suit our software needs
        To unify the logging and the headers
    '''

    levels = {
        'critical': logging.CRITICAL,
        'error': logging.ERROR,
        'warn': logging.WARNING,
        'warning': logging.WARNING,
        'info': logging.INFO,
        'debug': logging.DEBUG
    }

    _instance = False
    _logger = None
    _module_logger = None
    _log_file = None
    _str_formater = "%(asctime)s [%(levelname)-4.7s]  %(message)s"
    _consolehandler = None
    _reference = __name__

    @classmethod
    def set_up_log(cls, level='warning'):

        # Get root
        cls._module_logger = logging.getLogger(cls._reference)

        # Set up console log
        consolehandler = logging.StreamHandler()
        logformatter = logging.Formatter(cls._str_formater)
        consolehandler.setFormatter(logformatter)
        consolehandler.setLevel(level=cls.levels[level])
        cls._consolehandler = consolehandler
        cls._module_logger.addHandler(consolehandler)


    @classmethod
    def set_log_file(cls, log_file_name, absolute_path=False, log_level='info'):

        cls._module_logger = logging.getLogger(cls._reference)

        # If an absolute path is given just use that
        if not absolute_path:
            local_dir = Path(inspect.getfile(inspect.currentframe()))
            sim_dir = local_dir.absolute().parents[1]

            if log_file_name[0] == '/':
                log_file_name = log_file_name.replace('/', '', 1)

            log_file_name = os.path.join(sim_dir, log_file_name)

        os.makedirs(os.path.dirname(log_file_name), exist_ok=True)
        if os.path.isfile(log_file_name):
            os.remove(log_file_name)

        cls._log_file = log_file_name

        filehandler = logging.FileHandler(filename=cls._log_file)
        logformatter = logging.Formatter(cls._str_formater)
        filehandler.setFormatter(logformatter)
        filehandler.setLevel(cls.levels[log_level])
        cls._module_logger.addHandler(filehandler)

    @classmethod
    def set_log_level(cls, log_level):
        if cls._consolehandler is not None:
            cls._consolehandler.setLevel(level=cls.levels[log_level])
        else:
            raise ValueError('Please set up logger before adjusting level')

    @classmethod
    def set_log_reference(cls, reference):
        cls._reference = reference



if __name__ == '__main__':

    ProgramLogger.set_log_file('logs/last_logs.txt')
    ProgramLogger.write_log('I love refrigerators', level='warning')

