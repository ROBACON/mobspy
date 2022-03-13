from copy import deepcopy


def error(message):
    print('\033[91m' + 'ERROR: ' + message + '\033[0m')
    exit(1)


def debug(message=''):
    message_copy = deepcopy(message)
    print(message_copy.replace('_dot_', '.'))
    pass


def warning(message):
    message_copy = deepcopy(message)
    message.replace('_dot_', '.')
    print('\033[91m' + 'WARNING: ' + message_copy + '\033[0m')


if __name__ == '__main__':
    a = 5
    error(f'Did not found your {a}')