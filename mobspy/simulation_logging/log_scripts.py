from copy import deepcopy

import sys

global_simlog_level = 3


def error(message):
    print('\033[91m' + 'ERROR: ' + message + '\033[0m', file=sys.stderr)
    exit(1)


def debug(message=''):
    global global_simlog_level
    if global_simlog_level >= 2:
        message_copy = deepcopy(message)
        print(message_copy.replace('_dot_', '.'), file=sys.stderr)
    pass


def warning(message):
    global global_simlog_level
    if global_simlog_level >= 1:
        message_copy = deepcopy(message)
        message.replace('_dot_', '.')
        print('\033[91m' + 'WARNING: ' + message_copy + '\033[0m', file=sys.stderr)


if __name__ == '__main__':
    a = 5
    error(f'Did not found your {a}')
