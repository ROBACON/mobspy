"""Provide a proxy object that defers module import until first attribute access."""

from __future__ import annotations

from importlib import import_module as implib_import_module
from types import ModuleType
from typing import Any


class LazyImporter:
    def __init__(self, module_name: str, *attrs: str) -> None:
        self.module_name = module_name
        self.attrs = attrs
        self._module: ModuleType | None = None

    def __getattr__(self, attr: str) -> Any:
        if self._module is None:
            self._module = implib_import_module(self.module_name)
        return getattr(self._module, attr)
