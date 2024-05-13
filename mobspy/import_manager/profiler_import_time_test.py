import cProfile
from importlib import import_module as implib_import_module

def profile_imports():
    # Import the modules you want to profile
    implib_import_module('mobspy')


if __name__ == "__main__":

    cProfile.run("profile_imports()", sort='cumulative')