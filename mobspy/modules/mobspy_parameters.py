"""Define and manage symbolic parameters for use in reaction rates and sweeps."""

from __future__ import annotations

import linecache
import sys
from typing import Any

from mobspy.exceptions import ParameterError
from mobspy.mobspy_logging import get_logger

_logger = get_logger(__name__)
from pint import Quantity  # noqa: E402

from mobspy.modules.expression_nodes import (  # noqa: E402
    ParamRefNode,
)
from mobspy.modules.mobspy_expressions import (  # noqa: E402
    ExpressionDefiner as me_ExpressionDefiner,
)
from mobspy.modules.mobspy_expressions import (  # noqa: E402
    QuantityConverter as me_QuantityConverter,
)


class Internal_Parameter_Constructor(me_ExpressionDefiner, me_QuantityConverter):
    """Constructor called by ModelParameters to create a model parameter.

    Not a simulation parameter. The user is not supposed to create
    parameters using this object.
    """

    # convert_received_unit
    parameter_stack: dict[str, Internal_Parameter_Constructor] = {}

    def __init__(self, name: str, value: Any) -> None:
        self._generate_necessary_attributes()

        temp_set = set()
        temp_set.add(self)
        self.name = name
        self.original_value = value
        self.parameter_stack[name] = self

        self._ms_active = True

        self._operation = ParamRefNode(name)
        self._parameter_set.add(self)

        self.process_value(value)

    def unit_process(self, value: Quantity) -> tuple[Any, Any]:  # type: ignore[type-arg]
        """Convert a Pint quantity to MobsPy standard units.

        Returns:
            Tuple of (converted magnitude, original unit).
        """
        converted = self.convert_received_unit(value)
        self.value = converted.magnitude

        self.original_magnitude = value.magnitude
        self.conversion_factor = self.value / self.original_magnitude
        self.original_unit = value.units

        # Store the converted (MobsPy standard) unit for unit tracking.
        # This ensures concentration parameters (e.g. millimolar -> counts/dm³)
        # are compatible with species unit_conc_op (1/dm³) during addition.
        self._unit_count_op = converted
        self._unit_conc_op = converted
        self._unit_operation = converted
        self._has_units = True

        return self.value, self.original_unit

    def process_value(self, value: Any) -> None:
        """Store the parameter value, converting units if present."""
        if isinstance(value, Quantity):
            self.unit_process(value)
        elif isinstance(value, (list, tuple)):
            new_list: list[Any] = []
            first_unit: Any = None
            for i, val in enumerate(value):
                if isinstance(val, Quantity) and i == 0:
                    new_value, first_unit = self.unit_process(val)
                elif isinstance(val, Quantity) and i > 0 and first_unit is None:
                    raise ParameterError("MobsPy parameters must all be the same unit")
                elif isinstance(val, Quantity) and i > 0 and first_unit is not None:
                    new_value, unit = self.unit_process(val)
                    if unit != first_unit:
                        raise ParameterError(
                            "MobsPy parameters must all be the same unit"
                        )
                else:
                    new_value = val

                new_list.append(new_value)

            self.value = new_list
        else:
            self.value = value

            self._unit_count_op = 1
            self._unit_conc_op = 1
            self.conversion_factor = 1
            self._has_units = False

    def convert_to_original_unit(self) -> None:
        """Converts from MobsPy standard unit to the original unit."""
        if self.has_units():
            self.set_value(self.value / self.conversion_factor * self.original_unit)  # pyright: ignore[reportOperatorIssue]

    def rename(self, new_name: str) -> None:
        """Renames a parameter, checking name availability via the parameter stack."""
        if new_name in self.parameter_stack:
            _logger.warning(
                " MobsPy uses a parameter dictionary with parameter"
                " names as keys and the respective object as value"
                " to keep track of created parameters. As there is"
                " a parameter with this name already in the stack,"
                " the old will be deleted and replaced by this one."
            )

        del self.parameter_stack[self.name]
        self.parameter_stack[new_name] = self
        self.name = new_name

    def set_value(self, new_value: Any) -> Internal_Parameter_Constructor:
        """
        Sets value of parameter
        """
        self.value = new_value
        return self

    def has_units(self) -> bool:
        """
        Check if is a unit based parameter or not
        """
        return self._has_units

    def update_value(self, new_value: Any) -> None:
        """Replace the parameter value and reprocess units."""
        temp_set = set()
        temp_set.add(self)
        self.original_value = new_value

        self._ms_active = True

        self._parameter_set.add(self)

        self.process_value(new_value)

    def get_name(self) -> str:
        """Return the parameter name."""
        return self.name

    def __str__(self) -> str:
        return str(self._operation)


def ModelParameters(
    *args: Any,
) -> Internal_Parameter_Constructor | list[Internal_Parameter_Constructor]:
    """
    Creates ModelParameters. Like meta-species, it uses the
    variable names as parameter names
    """
    frame = sys._getframe(1)
    code_line = linecache.getline(frame.f_code.co_filename, frame.f_lineno).rstrip("\n")
    separated_line = code_line.split("=")[-2].replace(" ", "")
    parameter_variable_names = separated_line.split(",")

    if len(args) != len(parameter_variable_names):
        raise ParameterError(
            "You must provide an initial value for every parameter variable declared"
        )

    if len(parameter_variable_names) > 1:
        parameters_to_return: (
            Internal_Parameter_Constructor | list[Internal_Parameter_Constructor]
        ) = [
            Internal_Parameter_Constructor(p, v)
            for p, v in zip(parameter_variable_names, args)
        ]
    else:
        parameters_to_return = Internal_Parameter_Constructor(
            parameter_variable_names[0], args[0]
        )

    return parameters_to_return
