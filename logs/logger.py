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
    _log_level = levels['warning']
    _module_logger = None
    _log_file = None

    @classmethod
    def write_log(cls, log_text, level='info'):

        # The first time we call the logging class in the program, so set up the variables
        if not cls._instance:

            # Setup log format
            logformatter = logging.Formatter(f"%(asctime)s [%(levelname)-4.7s]  %(message)s ")

            # Setup log file
            if cls._log_file is not None:
                filehandler = logging.FileHandler(filename=cls._log_file)
                filehandler.setFormatter(logformatter)
                filehandler.setLevel(cls.levels['info'])
                cls._module_logger.addHandler(filehandler)

            # Set up console
            consolehandler = logging.StreamHandler()
            consolehandler.setFormatter(logformatter)
            consolehandler.setLevel(level=cls._log_level)
            cls._module_logger.addHandler(consolehandler)

            cls._instance = True

        # According to the level we call the correct logging function
        level = str(level.lower())
        if level == 'error':
            cls._module_logger.error(log_text)
        elif level == 'info':
            cls._module_logger.info(log_text)
        elif level == 'warning':
            cls._module_logger.warning(log_text)
        elif level == 'debug':
            cls._module_logger.debug(log_text)
        elif level == 'critical':
            cls._module_logger.critical(log_text)


    @classmethod
    def set_log_file(cls, log_file_name, absolute_path=False):

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

    @classmethod
    def set_log_level(cls, log_level):
        cls._log_level = cls.levels[log_level.lower()]


if __name__ == '__main__':

    ProgramLogger.set_log_file('logs/last_logs.txt')
    ProgramLogger.write_log('I love refrigerators', level='warning')

