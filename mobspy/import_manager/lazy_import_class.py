from typing import Any as type_Any
from importlib import import_module as implib_import_module

class LazyImporter:

    def __init__(self, module_name, *attrs):
        self.module_name = module_name
        self.attrs = attrs
        self._module: type_Any | None = None

    def __getattr__(self, attr):
        if self._module is None:
            self._module = implib_import_module(self.module_name)
        return getattr(self._module, attr)


if __name__ == '__main__':

    basico = LazyImporter('basico')

    a = basico.model_io.load_model_from_string('')

    print(a)

    # basico_model_io.load_model_from_string(sbml_str)


