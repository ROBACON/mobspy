import inspect


def Apply(function, *args):

    for arg in args:
        if type(arg) != tuple:
            raise TypeError('Only tuples are accepted for assignment with Apply')

    for arg in args:
        function(*arg)


if __name__ == '__main__':
    pass

