

def error(message):
    print('\033[91m' + 'ERROR: ' + message + '\033[0m')
    exit(1)


def debug(message=''):
    # print(message)
    pass


def warning(message):
    print('\033[91m' + 'WARNING: ' + message + '\033[0m')


if __name__ == '__main__':
    a = 5
    error(f'Did not found your {a}')