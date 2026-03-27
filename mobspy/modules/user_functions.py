from __future__ import annotations

from collections.abc import Callable
from typing import Any


def Apply(function: Callable[..., Any], *args: tuple[Any, ...]) -> None:
    """
    Applies one function to a series of tuples with function
    arguments. It is useful when applying the same function over
    and over to several arguments
    """
    for arg in args:
        if type(arg) != tuple:  # noqa: E721
            raise TypeError("Only tuples are accepted for assignment with Apply")

    for arg in args:
        function(*arg)


if __name__ == "__main__":
    pass
