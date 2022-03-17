from copy import deepcopy

import sys


def error(message):
    print('\033[91m' + 'ERROR: ' + message + '\033[0m', file=sys.stderr)
    exit(1)


def debug(message=''):
    message_copy = deepcopy(message)
    print(message_copy.replace('_dot_', '.'), file=sys.stderr)
    pass


def warning(message):
    message_copy = deepcopy(message)
    message.replace('_dot_', '.')
    print('\033[91m' + 'WARNING: ' + message_copy + '\033[0m', file=sys.stderr)


if __name__ == '__main__':
    a = 5
    error(f'Did not found your {a}')
