from __future__ import annotations

import inspect
import logging
from typing import Any

from mobspy.modules.mobspy_expressions import (
    ExpressionDefiner as me_ExpressionDefiner,
)
from mobspy.modules.mobspy_expressions import (
    u,
)

logger = logging.getLogger(__name__)


def expression_compilation_initiation() -> list[me_ExpressionDefiner]:
    """
    First find all ExpressionDefiner objects in the stack.
    Then, sets _ms_active to True to change the behavior of
    quantities and units in expressions
    """
    u._ms_active = True

    expressions_in_stack: list[me_ExpressionDefiner] = []

    for i in range(len(inspect.stack())):
        local_names = inspect.stack()[i][0].f_locals
        global_names = inspect.stack()[i][0].f_globals
        for key, item in global_names.items():  # noqa: B007
            if isinstance(item, me_ExpressionDefiner):
                expressions_in_stack.append(item)
        for key, item in local_names.items():  # noqa: B007
            if isinstance(item, me_ExpressionDefiner):
                expressions_in_stack.append(item)

    for expression in expressions_in_stack:
        expression._ms_active = True

    return expressions_in_stack


def expression_compilation_finish(expressions: list[me_ExpressionDefiner]) -> None:
    """
    Sets _ms_active to false for all expressions found in the
    stack by expression_compilation_initiation()

    :param expressions: all expressions found in stack at the
        moment of reaction compilation
    """
    u._ms_active = False

    for expression in expressions:
        expression._ms_active = False


class Unit_Context_Setter:
    def __enter__(self) -> None:
        self.expressions_in_stack = expression_compilation_initiation()

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: Any,
    ) -> None:
        expression_compilation_finish(self.expressions_in_stack)
