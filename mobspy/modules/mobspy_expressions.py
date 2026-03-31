"""Expression evaluation engine for unit-aware rate functions and species references."""

from __future__ import annotations

import contextlib
import re
from collections.abc import Callable
from copy import deepcopy
from typing import TYPE_CHECKING, Any

from numpy import (
    add as np_add,
)
from numpy import (
    divide as np_divide,
)
from numpy import (
    floating as np_float_,
)
from numpy import (
    integer as np_int_,
)
from numpy import (
    multiply as np_multiply,
)
from numpy import (
    subtract as np_subtract,
)
from pint import DimensionalityError, Quantity, UnitRegistry
from scipy.constants import N_A

from mobspy.constants import (
    CONCENTRATION_PREFIX,
    COUNT_PREFIX,
    NULL_SPECIES,
)
from mobspy.exceptions import CompilationError, UnitError
from mobspy.modules.expression_nodes import (  # noqa: F401
    BinaryOpNode,
    ExprNode,
    FunctionCallNode,
    LiteralNode,
    ParamRefNode,
    SpeciesRefNode,
    _render_resolved,
    _to_expr_node,
)
from mobspy.modules.species_operators import (  # noqa: F401
    Bool_Override,
    Specific_Species_Operator,
    _ms_active_ctx,
)
from mobspy.types import RenderContext

if TYPE_CHECKING:
    from numpy import ufunc as np_ufunc

    from mobspy.modules.meta_class import Species
    from mobspy.modules.model_unit_context import ModelUnitContext


_QUANTITY_OPS: dict[str, Callable[..., Any]] = {
    "__add__": Quantity.__add__,
    "__radd__": Quantity.__radd__,  # type: ignore[dict-item, misc]
    "__sub__": Quantity.__sub__,
    "__rsub__": Quantity.__rsub__,
    "__mul__": Quantity.__mul__,
    "__rmul__": Quantity.__rmul__,  # type: ignore[dict-item, misc]
    "__truediv__": Quantity.__truediv__,
    "__rtruediv__": Quantity.__rtruediv__,
}

_RAW_OPS: dict[str, Callable[..., Any]] = {
    "__add__": lambda a, b: a + b,
    "__radd__": lambda a, b: b + a,
    "__sub__": lambda a, b: a - b,
    "__rsub__": lambda a, b: b - a,
    "__mul__": lambda a, b: a * b,
    "__rmul__": lambda a, b: b * a,
    "__truediv__": lambda a, b: a / b,
    "__rtruediv__": lambda a, b: b / a,
    "__pow__": lambda a, b: a**b,
    "__rpow__": lambda a, b: b**a,
}


class ExpressionDefiner:
    """
    Defines a MobsPy expression. Encompasses what is necessary for
    something to act as an expression, including storing the
    operations an object has gone through.
    """

    _operation: Any
    _unit_count_op: Any
    _unit_conc_op: Any
    _force_expression_mode: bool
    _parameter_set: set[Any]
    _expression_variables: set[Any]
    _has_units: bool
    _count_in_model: bool
    _concentration_in_model: bool
    _count_in_expression: bool
    _concentration_in_expression: bool
    _dimension: int | None
    species_list_operation_order: list[Any]

    def non_expression_add(self, other: Any) -> Any:
        """Perform the operation outside expression mode.

        Same pattern applies to all other ``non_expression_*`` methods:
        they delegate to the underlying type's arithmetic when
        expression-building mode is not active.
        """
        raise NotImplementedError

    def non_expression_radd(self, other: Any) -> Any:
        raise NotImplementedError

    def non_expression_sub(self, other: Any) -> Any:
        raise NotImplementedError

    def non_expression_rsub(self, other: Any) -> Any:
        raise NotImplementedError

    def non_expression_mul(self, other: Any) -> Any:
        raise NotImplementedError

    def non_expression_rmul(self, other: Any) -> Any:
        raise NotImplementedError

    def non_expression_truediv(self, other: Any) -> Any:
        raise NotImplementedError

    def non_expression_rtruediv(self, other: Any) -> Any:
        raise NotImplementedError

    def non_expression_pow(self, other: Any) -> Any:
        raise NotImplementedError

    def non_expression_rpow(self, other: Any) -> Any:
        raise NotImplementedError

    def non_expression_neg(self) -> Any:
        raise NotImplementedError

    @property
    def _ms_active(self) -> bool:
        """True when expression-building mode is active.

        Combines the global context (set during compilation) with an
        instance-level override used by parameters and MobsPyExpression
        instances that are always in expression mode.
        """
        return self._force_expression_mode or _ms_active_ctx.get()

    @_ms_active.setter
    def _ms_active(self, value: bool) -> None:
        self._force_expression_mode = value

    @classmethod
    def execute_op(cls, first: Any, second: Any, operation: str) -> Any:
        """
        Executes a unit operation using the Quantity class

        Args:
            first: First number.
            second: Second number.
            operation: Operation to be executed.
        """
        # Always unwrap OverrideQuantity to plain Quantity to avoid:
        # 1. "different registries" errors in Pint operations
        # 2. ExpressionDefiner overrides when _ms_active is True
        raw_first = first.q_object if isinstance(first, OverrideQuantity) else first
        raw_second = second.q_object if isinstance(second, OverrideQuantity) else second

        if isinstance(raw_first, Quantity) and isinstance(raw_second, Quantity):
            q_op = _QUANTITY_OPS.get(operation)
            if q_op is not None:
                return q_op(raw_first, raw_second)

        raw_op = _RAW_OPS.get(operation)
        if raw_op is not None:
            return raw_op(raw_first, raw_second)

        raise TypeError("Non Valid Operation resulted in no q_object creation")

    @staticmethod
    def _normalize_substance(value: Any) -> Any:
        """Convert [substance] units to counts via QuantityConverter.

        If the value is a Quantity with [substance] dimension, convert
        mol -> N_A counts. Otherwise return as-is.
        """
        if isinstance(value, (Quantity, OverrideQuantity)):
            dim = (
                value.dimensionality
                if not isinstance(value, OverrideQuantity)
                else value.q_object.dimensionality
            )
            if "[substance]" in str(dim):
                try:
                    result = QuantityConverter.convert_received_unit(value)
                    return (
                        result.q_object
                        if isinstance(result, OverrideQuantity)
                        else result
                    )
                except (DimensionalityError, TypeError, ValueError):
                    pass
        return value

    def execute_quantity_op(self, other: Any, operation: str) -> tuple[Any, Any]:
        """Executes the unit based operation. Used to verify the
        unit of a MobsPy expression.

        When an operation fails due to [substance] incompatibility (e.g.,
        mol/L + 1/L), retries after converting mol -> N_A counts.

        Args:
            other: Other number (or expression) to execute the operation on.
            operation: String symbol of the operation to be executed.
        """
        self_count = self._unit_count_op
        self_conc = self._unit_conc_op

        if isinstance(other, ExpressionDefiner):
            other_count = other._unit_count_op
            other_conc = other._unit_conc_op
        else:
            other_count = other
            other_conc = other

        _unit_errors = (DimensionalityError, TypeError, ValueError)
        try:
            count_op = self.execute_op(self_count, other_count, operation)
        except _unit_errors:
            # Retry after normalizing [substance] -> counts
            try:
                count_op = self.execute_op(
                    self._normalize_substance(self_count),
                    self._normalize_substance(other_count),
                    operation,
                )
            except _unit_errors as e:
                count_op = e

        try:
            conc_op = self.execute_op(self_conc, other_conc, operation)
        except _unit_errors:
            # Retry after normalizing [substance] -> counts
            try:
                conc_op = self.execute_op(
                    self._normalize_substance(self_conc),
                    self._normalize_substance(other_conc),
                    operation,
                )
            except _unit_errors as e:
                conc_op = e

        return count_op, conc_op

    _NUMPY_ARRAY_ERR = (
        "MobsPy does not yet support array-wise numpy operations, only element-wise"
    )

    def __array_ufunc__(self, ufunc: np_ufunc, _: str, *inputs: Any) -> Any:
        """
        Handles numpy compatibility with MobsPy expressions.
        ufunc is the operation, and inputs are both the numpy
        element (arrays not yet supported) and the quantity object.

        Args:
            ufunc: Numpy operation.
            _: Method, not used, please don't erase or it becomes an input - will break
                function.
            inputs: Numpy element and the quantity object.
        """
        # Implement all numpy operations
        if ufunc == np_add:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return self.__radd__(inputs[0])
            raise CompilationError(self._NUMPY_ARRAY_ERR)
        if ufunc == np_subtract:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return self.__rsub__(inputs[0])
            raise CompilationError(self._NUMPY_ARRAY_ERR)
        if ufunc == np_multiply:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return self.__rmul__(inputs[0])
            raise CompilationError(self._NUMPY_ARRAY_ERR)
        if ufunc == np_divide:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return self.__rtruediv__(inputs[0])
            raise CompilationError(self._NUMPY_ARRAY_ERR)
        raise CompilationError("Numpy operation not yet supported by MobsPy")

    # T avoids problems with __getattr__ from units/quantities
    def __add__(self, other: Any) -> Any:
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, "__add__")
            return self.create_from_new_operation(other, "+", count_op, conc_op, True)
        return self.non_expression_add(other)

    def __radd__(self, other: Any) -> Any:
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, "__radd__")
            return self.create_from_new_operation(other, "+", count_op, conc_op, False)
        return self.non_expression_radd(other)

    def __sub__(self, other: Any) -> Any:
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, "__sub__")
            return self.create_from_new_operation(other, "-", count_op, conc_op, True)
        return self.non_expression_sub(other)

    def __rsub__(self, other: Any) -> Any:
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, "__rsub__")
            return self.create_from_new_operation(other, "-", count_op, conc_op, False)
        return self.non_expression_rsub(other)

    def __mul__(self, other: Any) -> Any:
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, "__mul__")
            return self.create_from_new_operation(other, "*", count_op, conc_op, True)
        return self.non_expression_mul(other)

    def __rmul__(self, other: Any) -> Any:
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, "__rmul__")
            return self.create_from_new_operation(other, "*", count_op, conc_op, False)
        return self.non_expression_rmul(other)

    def __truediv__(self, other: Any) -> Any:
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, "__truediv__")
            return self.create_from_new_operation(other, "/", count_op, conc_op, True)
        return self.non_expression_truediv(other)

    def __rtruediv__(self, other: Any) -> Any:
        if self._ms_active:
            count_op, conc_op = self.execute_quantity_op(other, "__rtruediv__")
            return self.create_from_new_operation(other, "/", count_op, conc_op, False)
        return self.non_expression_rtruediv(other)

    def __pow__(self, other: Any) -> Any:
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, "__pow__")
            return self.create_from_new_operation(other, "^", count_op, conc_op)
        return self.non_expression_pow(other)

    def __rpow__(self, other: Any) -> Any:
        if self._ms_active:
            other = check_if_non_expression_operated(other)
            count_op, conc_op = self.execute_quantity_op(other, "__rpow__")
            return self.create_from_new_operation(other, "^", count_op, conc_op, False)
        return self.non_expression_rpow(other)

    def __neg__(self) -> Any:
        if self._ms_active:
            other = check_if_non_expression_operated(-1)
            count_op, conc_op = self.execute_quantity_op(-1, "__mul__")
            return self.create_from_new_operation(other, "*", count_op, conc_op, False)
        return self.non_expression_neg()

    def combine_binary_attributes(self, other: Any, attribute: str) -> bool:
        """
        Or gates to binary True or False attributes from self and other

        Args:
            other: Other expression to execute operation.
            attribute: Attributed to be combined.
        """
        to_return = False
        try:
            if self.__dict__[attribute]:
                to_return = True
        except KeyError:
            pass
        except AttributeError:
            pass

        try:
            if other.__dict__[attribute]:
                to_return = True
        except KeyError:
            pass
        except AttributeError:
            pass

        return to_return

    def _generate_necessary_attributes(self) -> None:
        """
        Replacement for __init__ for classes inheriting from
        multiple Pint objects. Gives the object all necessary
        attributes to execute create_from_new_operation.
        """
        # Operation variables
        self._operation = None
        self._unit_count_op = 1
        self._unit_conc_op = 1

        # Per-instance override for always-active expression mode
        # (used by parameters and MobsPyExpression).
        # The global context var handles the compilation context.
        self._force_expression_mode = False

        # Parameter Variables
        self._parameter_set: set[Any] = set()
        self._expression_variables: set[Any] = set()

        # Presence of units
        self._has_units: bool = False

        # Concentration and counts
        self._count_in_model = False
        self._concentration_in_model = False
        self._count_in_expression = False
        self._concentration_in_expression = False

        # Dimension
        self._dimension = None

        self.species_list_operation_order: list[Any] = []

    def create_from_new_operation(  # noqa: PLR0913
        self,
        other: Any,
        symbol: str,
        count_op: Any,
        conc_op: Any,
        direct_sense: bool = True,
        operation: str | ExprNode | None = None,
    ) -> MobsPyExpression:
        """
        Executes the storing based operation. Also executes a unit
        operation to verify the unit of the given expression.

        Args:
            other: Other number (or expression) to execute the operation on.
            symbol: String symbol of the operation.
            count_op: Unit of the expression if the arguments are considered
                dimensionless.
            conc_op: Unit of the expression if the arguments are considered 1/v.
            direct_sense: Sense of the operation.
            operation: Current operation in the stack.
        """
        _has_units = _check_either_has_units(self, other)
        _validate_unit_ops(self, other, count_op, conc_op)

        if isinstance(self, (Quantity, OverrideQuantity)):
            self = QuantityConverter.convert_received_unit(self)  # type: ignore[assignment]  # noqa: PLW0642
        if isinstance(other, (Quantity, OverrideQuantity)):
            other = QuantityConverter.convert_received_unit(other)

        count_op = _safe_convert_unit_op(count_op)
        conc_op = _safe_convert_unit_op(conc_op)

        operation = _build_operation_node(
            self, other, symbol, count_op, direct_sense, operation
        )

        new_parameter_set = self._parameter_set
        new_expression_variables = self._expression_variables

        with contextlib.suppress(AttributeError):
            new_parameter_set = new_parameter_set.union(other._parameter_set)
        with contextlib.suppress(AttributeError):
            new_expression_variables = new_expression_variables.union(
                other._expression_variables
            )

        new_species_list_operation_order = list(
            getattr(self, "species_list_operation_order", [])
        )
        try:
            for spe in other.species_list_operation_order:
                if spe not in new_species_list_operation_order:
                    new_species_list_operation_order.append(spe)
        except AttributeError:
            pass

        _count_in_model = self.combine_binary_attributes(other, "_count_in_model")
        _concentration_in_model = self.combine_binary_attributes(
            other, "_concentration_in_model"
        )
        _count_in_expression = self.combine_binary_attributes(
            other, "_count_in_expression"
        )
        _concentration_in_expression = self.combine_binary_attributes(
            other, "_concentration_in_expression"
        )

        if _count_in_model and _concentration_in_model:
            raise UnitError(
                "A meta-species in a model cannot be both a count and a concentration"
            )

        if _count_in_expression and _concentration_in_expression:
            raise UnitError(
                "A meta-species in an expression cannot be "
                "both a count and a concentration"
            )

        for variable in self._expression_variables:
            variable._operation = SpeciesRefNode(variable.species_string)

        dimension = _resolve_dimension(self, other)

        return MobsPyExpression(
            species_string=NULL_SPECIES,
            species_object=None,
            operation=operation,
            unit_count_op=count_op,
            unit_conc_op=conc_op,
            dimension=dimension,
            expression_variables=new_expression_variables,
            parameter_set=new_parameter_set,
            count_in_model=_count_in_model,
            concentration_in_model=_concentration_in_model,
            count_in_expression=_count_in_expression,
            concentration_in_expression=_concentration_in_expression,
            has_units=_has_units,
            species_list_operation_order=new_species_list_operation_order,
        )


def _check_either_has_units(self_obj: Any, other: Any) -> bool:
    """Return True if either operand has units."""
    _has_units = False
    try:
        if self_obj._has_units:
            _has_units = True
        if other._has_units:
            _has_units = True
    except AttributeError:
        pass
    return _has_units


def _validate_unit_ops(self_obj: Any, other: Any, count_op: Any, conc_op: Any) -> None:
    """Raise UnitError if both count and concentration unit operations failed."""
    try:
        c1 = self_obj._has_units
    except AttributeError:
        c1 = False
    try:
        c2 = other._has_units
    except AttributeError:
        c2 = False
    if (
        (c1 or c2)
        and isinstance(count_op, Exception)
        and isinstance(conc_op, Exception)
    ):
        raise UnitError(
            "Incompatible units in expression. "
            "Both the count and concentration interpretations failed:\n"
            f"  count path: {count_op}\n"
            f"  concentration path: {conc_op}\n"
            "Check that all terms in additions/subtractions "
            "have matching dimensions. "
            "Example of an invalid expression: "
            "(1/u.s) + (1*u.l/u.s)"
        )


def _safe_convert_unit_op(op: Any) -> Any:
    """Convert a Quantity unit op via QuantityConverter.

    Captures conversion failures as exceptions.
    """
    try:
        if isinstance(op, Quantity):
            op = QuantityConverter.convert_received_unit(op)
    except (DimensionalityError, TypeError, ValueError) as e:
        op = e
    return op


def _build_operation_node(  # noqa: PLR0913
    self_obj: Any,
    other: Any,
    symbol: str,
    count_op: Any,
    direct_sense: bool,
    operation: Any,
) -> Any:
    """Build the AST operation node or extract numeric magnitude."""
    self_is_symbolic = isinstance(self_obj._operation, (str, ExprNode))
    other_is_symbolic = isinstance(other, ExpressionDefiner) and isinstance(
        other._operation, (str, ExprNode)
    )

    if self_is_symbolic or other_is_symbolic:
        node_self = (
            _to_expr_node(self_obj.magnitude)
            if isinstance(self_obj, (Quantity, OverrideQuantity))
            else _to_expr_node(self_obj._operation)
            if isinstance(self_obj._operation, ExprNode)
            else _to_expr_node(str(self_obj))
        )
        node_other = (
            _to_expr_node(other.magnitude)
            if isinstance(other, (Quantity, OverrideQuantity))
            else _to_expr_node(other._operation)
            if isinstance(other, ExpressionDefiner)
            and isinstance(other._operation, ExprNode)
            else _to_expr_node(str(other))
        )

        if direct_sense:
            left, right = node_self, node_other
        else:
            left, right = node_other, node_self

        if operation is None:
            operation = BinaryOpNode(left, symbol, right)
    else:
        if isinstance(count_op, Exception):
            raise count_op
        operation = count_op.magnitude
    return operation


def _resolve_dimension(self_obj: Any, other: Any) -> int | None:
    """Merge dimension from two expression operands."""
    dimension_1 = (
        self_obj._dimension if isinstance(self_obj, MobsPyExpression) else None
    )
    dimension_2 = other._dimension if isinstance(other, MobsPyExpression) else None

    if dimension_1 is not None and dimension_2 is None:
        return dimension_1
    if dimension_1 is None and dimension_2 is not None:
        return dimension_2
    if dimension_1 is not None and dimension_2 is not None:
        if dimension_1 != dimension_2:
            raise TypeError(
                "Dimensions are inconsistent between "
                "different MobsPy expression objects"
            )
        return dimension_1
    return None


class OverrideUnitRegistry:
    """
    It was necessary to override Pint's unit registry so it
    behaves one way under context and normally without context
    """

    def __init__(self) -> None:
        self.unit_registry_object = UnitRegistry()

    def __call__(self, *args: Any, **kwargs: Any) -> OverrideQuantity:
        q_object = self.unit_registry_object(*args, **kwargs)
        return OverrideQuantity(q_object)  # pyright: ignore[reportReturnType]

    def __getattr__(self, item: str) -> OverrideQuantity:
        if item == "h":
            item = "hour"

        q_object = 1 * self.unit_registry_object.__getattr__(item)
        return OverrideQuantity(q_object)  # pyright: ignore[reportReturnType]


# u is defined here
u = OverrideUnitRegistry()


class QuantityConverter:
    """
    Class that converts an unit from any to MobsPy L-s-counts conversion
    """

    @classmethod
    def convert_received_unit(
        cls,
        quantity: Quantity | OverrideQuantity,
        model_context: ModelUnitContext | None = None,
    ) -> Quantity | OverrideQuantity:
        """
        Converts a received quantity to L-s-counts, standard MobsPy units

        Args:
            quantity: Received quantity to convert.
        """
        is_override = isinstance(quantity, OverrideQuantity)
        copied_quantity = deepcopy(quantity.q_object if is_override else quantity)

        length_repl, time_repl, has_model_substance = cls._base_units_from_context(
            model_context
        )

        to_convert_into, copied_quantity = cls._apply_dimension_replacements(
            copied_quantity,
            length_repl,
            time_repl,
            has_model_substance,
            model_context,
            quantity,
        )

        copied_quantity.ito(to_convert_into)

        if not is_override:
            return copied_quantity  # type: ignore[no-any-return]  # pyright: ignore[reportReturnType]
        return OverrideQuantity(copied_quantity)  # pyright: ignore[reportReturnType]

    @classmethod
    def _base_units_from_context(
        cls,
        model_context: ModelUnitContext | None,
    ) -> tuple[str, str, bool]:
        """Extract base unit replacement strings from model context."""
        ur = u.unit_registry_object
        length_repl = "dm"
        time_repl = "s"
        has_model_substance = False

        if model_context is not None:
            vol_q = ur.Quantity(1, str(model_context.volume_unit))
            vol_dim = dict(vol_q.dimensionality)
            length_exp = int(vol_dim.get("[length]", 3))
            if length_exp != 0:
                base_length_q = vol_q ** (1.0 / length_exp)
                length_repl = str(base_length_q.units)
            time_repl = str(model_context.time_unit)
            has_model_substance = model_context.substance_is_molar

        return length_repl, time_repl, has_model_substance

    @classmethod
    def _apply_dimension_replacements(  # noqa: PLR0913
        cls,
        copied_quantity: Any,
        length_repl: str,
        time_repl: str,
        has_model_substance: bool,
        model_context: ModelUnitContext | None,
        original_quantity: Any,
    ) -> tuple[str, Any]:
        """Replace dimension placeholders in the dimensionality string."""
        ur = u.unit_registry_object
        to_convert_into = str(copied_quantity.dimensionality)

        _simple_replacements: dict[str, str] = {
            "[length]": length_repl,
            "[time]": time_repl,
            "[temperature]": "K",
            "[mass]": "kg",
        }
        for dim_key, repl in _simple_replacements.items():
            if dim_key in to_convert_into:
                to_convert_into = to_convert_into.replace(dim_key, repl)

        if "[substance]" in to_convert_into:
            mol_power = int(dict(copied_quantity.dimensionality).get("[substance]", 1))
            if has_model_substance:
                sub_repl = str(model_context.substance_unit)  # type: ignore[union-attr]
                to_convert_into = to_convert_into.replace("[substance]", sub_repl)
            else:
                copied_quantity = copied_quantity * (N_A / (1 * ur.mol)) ** mol_power
                if "[substance]" in str(copied_quantity.dimensionality):
                    raise TypeError(
                        f"Could not convert molar quantity {original_quantity} "
                        f"(substance power={mol_power})"
                    )
                to_convert_into = to_convert_into.replace("[substance]", "1")

        return to_convert_into, copied_quantity


class OverrideQuantity(ExpressionDefiner, Quantity):
    """Pint Quantity subclass that participates in MobsPy expression building.

    Wraps a Quantity so that arithmetic operators build symbolic rate
    expressions when expression-building mode is active, while still
    supporting unit conversions when it is not.
    """

    def __array_ufunc__(
        self,
        ufunc: np_ufunc,
        _: str,
        *inputs: Any,
    ) -> OverrideQuantity | None:
        """
        Handles numpy compatibility with MobsPy expressions.
        ufunc is the operation, and inputs are both the numpy
        element (arrays not yet supported) and the quantity
        object.

        Args:
            ufunc: Numpy operation.
            _: Method, not used, please don't erase or it becomes an input - will break
                function.
            inputs: Numpy element and the quantity object.
        """
        # Implement all numpy operations
        if ufunc == np_add:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return OverrideQuantity(float(inputs[0]) + self.q_object)  # pyright: ignore[reportReturnType]
            raise CompilationError(
                "MobsPy does not yet support array-wise "
                "numpy operations, only element-wise"
            )
        if ufunc == np_subtract:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return OverrideQuantity(float(inputs[0]) - self.q_object)  # pyright: ignore[reportReturnType]
            raise CompilationError(
                "MobsPy does not yet support array-wise "
                "numpy operations, only element-wise"
            )
        if ufunc == np_multiply:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return OverrideQuantity(float(inputs[0]) * self.q_object)  # pyright: ignore[reportReturnType]
            raise CompilationError(
                "MobsPy does not yet support array-wise "
                "numpy operations, only element-wise"
            )
        if ufunc == np_divide:
            if isinstance(inputs[0], (np_int_, np_float_)):
                return OverrideQuantity(float(inputs[0]) / self.q_object)  # pyright: ignore[reportReturnType]
            raise CompilationError(
                "MobsPy does not yet support array-wise "
                "numpy operations, only element-wise"
            )
        raise CompilationError("Numpy operation not yet supported by MobsPy")
        return None

    def non_expression_add(self, other: Any) -> OverrideQuantity | Any:
        """Perform unit-aware addition outside expression mode.

        Same pattern applies to all other ``non_expression_*`` methods
        on OverrideQuantity: they delegate to Pint's Quantity arithmetic
        when expression-building mode is not active.
        """
        if isinstance(other, OverrideQuantity):
            # Don't delegate to other.__radd__ to avoid infinite loops
            # Just perform the operation directly
            q_object = Quantity.__add__(self.q_object, other.q_object)
        elif isinstance(other, ExpressionDefiner):
            return other.__radd__(self)
        else:
            q_object = Quantity.__add__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_radd(self, other: Any) -> OverrideQuantity | Any:
        if isinstance(other, OverrideQuantity):
            # Don't delegate to other.__add__ to avoid infinite loops
            # Just perform the operation directly
            q_object = Quantity.__radd__(self.q_object, other.q_object)  # type: ignore[misc]
        elif isinstance(other, ExpressionDefiner):
            return other.__add__(self)
        else:
            q_object = Quantity.__radd__(self.q_object, other)  # type: ignore[misc]
        return OverrideQuantity(q_object)

    def non_expression_sub(self, other: Any) -> OverrideQuantity | Any:
        if isinstance(other, OverrideQuantity):
            # Don't delegate to other.__rsub__ to avoid infinite loops
            # Just perform the operation directly
            q_object = Quantity.__sub__(self.q_object, other.q_object)
        elif isinstance(other, ExpressionDefiner):
            return other.__rsub__(self)
        else:
            q_object = Quantity.__sub__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_rsub(self, other: Any) -> OverrideQuantity | Any:
        if isinstance(other, OverrideQuantity):
            # Don't delegate to other.__sub__ to avoid infinite loops
            # Just perform the operation directly
            q_object = Quantity.__rsub__(self.q_object, other.q_object)
        elif isinstance(other, ExpressionDefiner):
            return other.__sub__(self)
        else:
            q_object = Quantity.__rsub__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_mul(self, other: Any) -> OverrideQuantity | Any:
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__mul__(self.q_object, other.q_object)
        elif isinstance(other, ExpressionDefiner):
            # If multiplying with an expression, delegate to its __rmul__
            return other.__rmul__(self)
        else:
            q_object = Quantity.__mul__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_rmul(self, other: Any) -> OverrideQuantity | Any:
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__rmul__(self.q_object, other.q_object)  # type: ignore[misc]
        elif isinstance(other, ExpressionDefiner):
            # If multiplying with an expression, delegate to its __mul__
            return other.__mul__(self)
        else:
            q_object = Quantity.__rmul__(self.q_object, other)  # type: ignore[misc]
        return OverrideQuantity(q_object)

    def non_expression_truediv(self, other: Any) -> OverrideQuantity | Any:
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__truediv__(self.q_object, other.q_object)
        elif isinstance(other, ExpressionDefiner):
            # If dividing by an expression, delegate to its __rtruediv__
            return other.__rtruediv__(self)
        else:
            q_object = Quantity.__truediv__(self.q_object, other)
        return OverrideQuantity(q_object)

    def non_expression_rtruediv(self, other: Any) -> OverrideQuantity:
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__rtruediv__(self.q_object, other.q_object)
        else:
            q_object = Quantity.__rtruediv__(self.q_object, other)
        return OverrideQuantity(q_object)  # pyright: ignore[reportReturnType]

    def non_expression_pow(self, other: Any) -> OverrideQuantity:
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__pow__(self.q_object, other.q_object)
        else:
            q_object = Quantity.__pow__(self.q_object, other)
        return OverrideQuantity(q_object)  # pyright: ignore[reportReturnType]

    def non_expression_rpow(self, other: Any) -> OverrideQuantity:
        if isinstance(other, OverrideQuantity):
            q_object = Quantity.__rpow__(self.q_object, other.q_object)
        else:
            q_object = Quantity.__rpow__(self.q_object, other)
        return OverrideQuantity(q_object)  # pyright: ignore[reportReturnType]

    def __init__(self, quantity_object: Quantity) -> None:
        self._generate_necessary_attributes()
        self.q_object = quantity_object

        self._unit_count_op = quantity_object
        self._unit_conc_op = quantity_object

        for key, item in quantity_object.__dict__.items():
            self.__dict__[key] = item

        self._operation = self.q_object.magnitude
        self._expression_variables: set[Any] = set()
        self._parameter_set: set[Any] = set()
        self._has_units: bool = True

    def __str__(self) -> str:
        # Always return the full quantity with units for user-facing output
        return str(self.q_object)

    def __repr__(self) -> str:
        return str(self.q_object)

    def convert(self, unit: str) -> OverrideQuantity:
        """Convert to a different unit, returning a new OverrideQuantity."""
        new_q_object = self.q_object.to(unit)
        return OverrideQuantity(new_q_object)  # pyright: ignore[reportReturnType]

    def convert_into(self, unit: str) -> None:
        """Convert this quantity to a different unit in-place."""
        self.q_object.ito(unit)


class MobsPyExpression(Specific_Species_Operator, ExpressionDefiner):
    """
    MobsPy expression objects are passed to rate functions
    to store the expressions they've been through.
    """

    def __str__(self) -> str:
        return str(self._operation)

    def __init__(  # noqa: PLR0913
        self,
        species_string: str,
        species_object: Species | None,
        operation: Any = None,
        unit_count_op: Any = 1,
        unit_conc_op: Any = (1 / u.unit_registry_object.liter),
        dimension: int | None = None,
        expression_variables: set[Any] | None = None,
        parameter_set: set[Any] | None = None,
        count_in_model: bool = True,
        concentration_in_model: bool = False,
        count_in_expression: bool = True,
        concentration_in_expression: bool = False,
        has_units: bool = False,
        species_list_operation_order: list[Any] | None = None,
        model_context: ModelUnitContext | None = None,
    ) -> None:
        super().__init__(species_string, species_object)
        self._generate_necessary_attributes()
        self._model_context = model_context

        if species_list_operation_order is None:
            self.species_list_operation_order = []
        else:
            self.species_list_operation_order = species_list_operation_order

        self._ms_active = True

        self._dimension = dimension

        if expression_variables is None:
            self._expression_variables: set[Any] = set()
            self._expression_variables.add(self)
        else:
            self._expression_variables = expression_variables

        if operation is None:
            self._operation: ExprNode | int | float = SpeciesRefNode(
                self.species_string
            )
        else:
            self._operation = operation

        if parameter_set is None:
            self._parameter_set: set[Any] = set()
        else:
            self._parameter_set = parameter_set

        self._count_in_model = count_in_model
        self._concentration_in_model = concentration_in_model
        self._count_in_expression = count_in_expression
        self._concentration_in_expression = concentration_in_expression

        self._unit_count_op = unit_count_op
        self._unit_conc_op = unit_conc_op

        if self._concentration_in_model:
            self._unit_conc_op = None

        self._has_units: bool = has_units

    def __getattr__(self, item: str) -> Specific_Species_Operator:
        return super().__getattr__(item)

    # Write string operation return and unit
    def generate_string_operation(
        self,
        skip_check: bool = False,
        reaction_order: int | None = None,
        dimension: int | None = None,
        model_context: ModelUnitContext | None = None,
    ) -> tuple[str, bool]:
        """
        Converts the expression from what has been stored
        to a string format for the sbml file.

        Args:
            skip_check: Skip units check or not - always set to False, True only for
                debugging.
            reaction_order: Order of the reaction - to check if the unit is correct.
        """
        if dimension is None:
            dimension = 3 if self._dimension is None else self._dimension

        operation = str(self._operation)

        if skip_check:
            return operation, True

        if not self._has_units:
            self._count_in_expression = True

        _ctx = model_context or getattr(self, "_model_context", None)
        _time_u, _vol_u = _resolve_validation_units(_ctx)

        early = self._check_constant_expression_units(
            operation, reaction_order, dimension, _time_u, _vol_u
        )
        if early is not None:
            return early

        self._resolve_expression_mode(dimension, _time_u, _vol_u)

        if isinstance(self._operation, ExprNode):
            convert_operation = self._convert_ast_path()
        else:
            convert_operation = self._convert_legacy_string_path()

        return convert_operation, self._count_in_expression

    def _check_constant_expression_units(
        self,
        operation: str,
        reaction_order: int | None,
        dimension: int,
        _time_u: Any,
        _vol_u: Any,
    ) -> tuple[str, bool] | None:
        """Check unit match for constant (no-variable) expressions.

        Returns early result or None.
        """
        if (
            self._has_units
            and self._expression_variables == set()
            and reaction_order is not None
        ):
            if self._unit_count_op.units == (1 / _time_u).units:  # pyright: ignore[reportAttributeAccessIssue]
                return operation, True
            if self._unit_conc_op.units == (  # pyright: ignore[reportAttributeAccessIssue]
                _vol_u ** (dimension * (reaction_order - 1)) / _time_u
            ):
                return operation, False
        return None

    def _resolve_expression_mode(
        self,
        dimension: int,
        _time_u: Any,
        _vol_u: Any,
    ) -> None:
        """Determine whether expression uses counts or concentrations."""
        c1 = self._has_units
        c2 = isinstance(self._unit_count_op, Exception)
        c3 = isinstance(self._unit_conc_op, Exception)

        if c1 and (c2 and c3):
            raise TypeError(
                "Unit resolution failed for reaction rate. "
                "Both count and concentration interpretations produced errors:\n"
                f"  count: {self._unit_count_op}\n"
                f"  concentration: {self._unit_conc_op}\n"
                "Verify that all arithmetic operations have compatible dimensions."
            )

        if c1 and not c3:
            conc_units = self._unit_conc_op.units
            is_vol_rate = conc_units == (1 / (_time_u * _vol_u**dimension))
            is_time_rate = conc_units == (1 / _time_u)
            if is_vol_rate or is_time_rate:
                self._concentration_in_expression = True

        if c1 and not c2 and self._unit_count_op.units == (1 / _time_u):
            self._count_in_expression = True

        if (
            self._has_units
            and not self._count_in_expression
            and not self._concentration_in_expression
        ):
            raise TypeError(
                "Could not determine whether the rate expression "
                "uses counts or concentrations.\n"
                f"  count interpretation: {self._unit_count_op}\n"
                f"  concentration interpretation: {self._unit_conc_op}\n"
                "The rate must resolve to 1/[time] for counts "
                "or 1/([time]*[volume]) for concentrations."
            )

    def _convert_ast_path(self) -> str:
        """Resolve expression via AST node tree."""
        if (
            not self._count_in_model
            and not self._concentration_in_model
            and not self._count_in_expression
            and not self._concentration_in_expression
        ):
            raise ValueError(
                "The expression did not resolve for "
                "lack of concentration/count "
                "specifications"
            )

        convert_operation = _render_resolved(
            self._operation,
            self._expression_variables,
            RenderContext(
                count_in_model=self._count_in_model,
                concentration_in_model=self._concentration_in_model,
                count_in_expression=self._count_in_expression,
                concentration_in_expression=self._concentration_in_expression,
            ),
        )

        expr_vars = self._expression_variables
        n_vars = len(expr_vars) if expr_vars else 0
        if (
            self._count_in_model
            and self._concentration_in_expression
            and not self._count_in_expression
        ):
            for _ in range(max(n_vars, 1)):
                convert_operation = "(" + convert_operation + ")" + "*volume"
        elif (
            self._concentration_in_model
            and self._count_in_expression
            and not self._concentration_in_expression
        ):
            for _ in range(max(n_vars, 1)):
                convert_operation = "(" + convert_operation + ")" + "/volume"
        return convert_operation

    def _convert_legacy_string_path(self) -> str:
        """Resolve expression via legacy string replacement."""
        convert_operation = str(self._operation)

        for variable in self._expression_variables:
            convert_operation = _apply_legacy_variable_conversion(
                convert_operation,
                variable,
                self._count_in_model,
                self._concentration_in_model,
                self._count_in_expression,
                self._concentration_in_expression,
            )

        return convert_operation


def _resolve_validation_units(
    ctx: ModelUnitContext | None,
) -> tuple[Any, Any]:
    """Return (time_unit, volume_base_unit) from a model context or defaults."""
    ur = u.unit_registry_object
    _time_u = ur.second
    _vol_u = ur.decimeter
    if ctx is not None:
        _time_u = ur.Quantity(1, str(ctx.time_unit)).units  # type: ignore[assignment]
        vol_q = ur.Quantity(1, str(ctx.volume_unit))
        vol_dim = dict(vol_q.dimensionality)
        length_exp = int(vol_dim.get("[length]", 3))
        if length_exp != 0:
            _vol_u = (vol_q ** (1.0 / length_exp)).units
        else:
            _vol_u = ur.decimeter
    return _time_u, _vol_u


def _apply_legacy_variable_conversion(  # noqa: PLR0913
    convert_operation: str,
    variable: Any,
    count_in_model: bool,
    concentration_in_model: bool,
    count_in_expression: bool,
    concentration_in_expression: bool,
) -> str:
    """Apply count/concentration name substitutions for one variable.

    Used in the legacy string-based expression resolution mode.
    """
    count_name = COUNT_PREFIX + variable.species_string
    concentration_name = CONCENTRATION_PREFIX + variable.species_string
    default_name = variable.species_string
    replace_name = variable.species_string

    if count_in_model:
        convert_operation = replace_spe_in_expr(
            convert_operation, count_name, replace_name
        )
        convert_operation = replace_spe_in_expr(
            convert_operation,
            concentration_name,
            "(" + replace_name + "/volume)",
        )
        if count_in_expression:
            pass
        elif concentration_in_expression:
            convert_operation = replace_spe_in_expr(
                convert_operation,
                default_name,
                "(" + replace_name + "/volume)",
            )
            convert_operation = "(" + convert_operation + ")" + "*volume"
        else:
            raise ValueError(
                "The expression did not resolve for "
                "lack of concentration/count "
                "specifications"
            )
    elif concentration_in_model:
        convert_operation = replace_spe_in_expr(
            convert_operation, concentration_name, replace_name
        )
        convert_operation = convert_operation.replace(
            count_name, "(" + replace_name + "*volume)"
        )
        convert_operation = replace_spe_in_expr(
            convert_operation,
            count_name,
            "(" + replace_name + "*volume)",
        )
        if count_in_expression:
            convert_operation = convert_operation.replace(
                default_name, "(" + replace_name + "*volume)"
            )
            convert_operation = replace_spe_in_expr(
                convert_operation,
                default_name,
                "(" + replace_name + "*volume)",
            )
            convert_operation = "(" + convert_operation + ")" + "/volume"
        elif concentration_in_expression:
            pass
        else:
            raise ValueError(
                "The expression did not resolve for "
                "lack of concentration/count "
                "specifications"
            )

    return convert_operation


def check_if_non_expression_operated(other: Any) -> Any:
    """Wrap non-numeric, non-expression operands as MobsPyExpression."""
    if (
        not isinstance(other, ExpressionDefiner)
        and not isinstance(other, int)
        and not isinstance(other, float)
        and not isinstance(other, Quantity)
        and not isinstance(other, (np_int_, np_float_))
    ):
        asg_node = SpeciesRefNode(str(other), mode="assignment")
        other = MobsPyExpression(
            str(asg_node),
            None,
            operation=asg_node,
            dimension=None,
            count_in_model=True,
            concentration_in_model=False,
            count_in_expression=False,
            concentration_in_expression=False,
            species_list_operation_order=[other]
            if hasattr(other, "get_spe_object")
            else [],
        )
    return other


def replace_spe_in_expr(string: str, to_replace: str, replacement: str) -> str:
    """Replace a species name in an expression string at word boundaries.

    Examples:
        >>> replace_spe_in_expr("A + B * A", "A", "X")
        'X + B * X'
    """
    pattern = re.compile(re.escape(to_replace) + r"(?![a-zA-Z0-9_])")
    return pattern.sub(replacement, string)


def _set_species_mode_in_tree(operation: Any, species_string: str, mode: str) -> Any:
    """Walk an expression tree/string and set the mode on matching SpeciesRefNodes."""
    if isinstance(operation, SpeciesRefNode):
        if operation.species_string == species_string:
            operation.mode = mode
    elif isinstance(operation, BinaryOpNode):
        _set_species_mode_in_tree(operation.left, species_string, mode)
        _set_species_mode_in_tree(operation.right, species_string, mode)
    elif isinstance(operation, FunctionCallNode):
        _set_species_mode_in_tree(operation.arg, species_string, mode)
    elif isinstance(operation, str):
        # Fallback for legacy string operations
        operation = operation.replace(species_string, f"${mode}${species_string}")
    return operation


class _Count_Base:  # noqa: N801
    def __getitem__(self, item: MobsPyExpression) -> MobsPyExpression | None:
        try:
            for v in item._expression_variables:
                if isinstance(item._operation, ExprNode):
                    _set_species_mode_in_tree(
                        item._operation,
                        v.species_string,
                        "count",
                    )
                else:
                    item._operation = str(item._operation).replace(  # type: ignore[assignment]
                        v.species_string, COUNT_PREFIX + v.species_string
                    )
            return item
        except AttributeError as e:
            raise CompilationError(
                "Count[] operator can only be used in MobsPy expressions"
            ) from e


Count = _Count_Base()


class _Conc_Base:  # noqa: N801
    def __getitem__(self, item: MobsPyExpression) -> MobsPyExpression | None:
        try:
            for v in item._expression_variables:
                if isinstance(item._operation, ExprNode):
                    _set_species_mode_in_tree(
                        item._operation, v.species_string, "concentration"
                    )
                else:
                    item._operation = str(item._operation).replace(  # type: ignore[assignment]
                        v.species_string, CONCENTRATION_PREFIX + v.species_string
                    )
            return item
        except AttributeError as e:
            raise CompilationError(
                "Concentration[] operator can only be used in MobsPy expressions"
            ) from e
        return None


Concentration = _Conc_Base()
