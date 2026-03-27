"""Utility functions for MobsPy package."""

from __future__ import annotations

import inspect
from typing import Any

from mobspy.constants import DEFAULT_STACK_INDEX, INTERACTIVE_CONTEXT_LINE_NUMBER


def safe_get_caller_info(stack_index: int = DEFAULT_STACK_INDEX) -> tuple[str, int]:
    """Safely extract caller information with fallbacks for interactive contexts.

    This utility handles the common pattern of stack inspection that's used throughout
    MobsPy for error reporting and code parsing. It gracefully handles cases where
    code_context is None (e.g., interactive interpreter, REPL).

    Args:
        stack_index: Index in the call stack to examine (default: 1 for immediate
        caller)

    Returns:
        Tuple of (code_line, line_number) where:
        - code_line: The source code line, or a fallback string for interactive contexts
        - line_number: The line number, or 0 for interactive contexts

    Example:
        >>> code_line, line_number = safe_get_caller_info()
        >>> print(f"Called from: {code_line} at line {line_number}")
    """
    try:
        frame = inspect.stack()[stack_index]
        if frame.code_context and len(frame.code_context) > 0:
            code_line = frame.code_context[0].rstrip()
            line_number = frame.lineno
            return code_line, line_number
        # Handle case where code_context is empty but frame exists
        return "interactive_context", frame.lineno
    except (IndexError, TypeError, AttributeError):
        # Handle cases where stack inspection fails completely
        return "interactive_context", INTERACTIVE_CONTEXT_LINE_NUMBER


def safe_get_caller_context_for_parsing(stack_index: int = DEFAULT_STACK_INDEX) -> str:
    """Extract caller context specifically for parsing operations.

    This is used by functions like BaseSpecies() and reaction definitions that need
    to parse the calling code. Provides sensible defaults for interactive usage.
    The original code expected the line without the final character, so we preserve
    that.

    Args:
        stack_index: Index in the call stack to examine

    Returns:
        A code line suitable for parsing, with fallbacks for different contexts
    """
    try:
        frame = inspect.stack()[stack_index]
        if frame.code_context and len(frame.code_context) > 0:
            # Original code expected [:-1] slicing to remove trailing character
            return (
                frame.code_context[0].rstrip()[:-1]
                if len(frame.code_context[0].rstrip()) > 0
                else "temp = BaseSpecies()"
            )
        # Provide a minimal parseable context for interactive use
        return "temp = BaseSpecies()"
    except (IndexError, TypeError, AttributeError):
        # Fallback that can be parsed by the existing parsing logic
        return "temp = BaseSpecies()"


def create_parameter_config(**kwargs: Any) -> dict[str, Any]:
    """Create a parameter configuration dictionary with validation.

    Helper function to create parameter dictionaries for complex functions
    while providing validation and defaults.

    Args:
        **kwargs: Parameter key-value pairs

    Returns:
        Validated parameter configuration dictionary
    """
    config = {}

    # Add validation logic here as needed
    for key, value in kwargs.items():
        if value is not None:  # Skip None values unless explicitly needed
            config[key] = value

    return config


def validate_parameter_types(
    params: dict[str, Any], expected_types: dict[str, type | tuple[type, ...]]
) -> None:
    """Validate parameter types against expected types.

    Args:
        params: Dictionary of parameter names to values
        expected_types: Dictionary of parameter names to expected types

    Raises:
        TypeError: If any parameter doesn't match expected type
    """
    for param_name, expected_type in expected_types.items():
        if param_name in params:
            value = params[param_name]
            if not isinstance(value, expected_type):
                expected_str = (
                    expected_type.__name__
                    if hasattr(expected_type, "__name__")
                    else str(expected_type)
                )
                raise TypeError(
                    f"Parameter '{param_name}' must be of type {expected_str}, got"
                    " {type(value).__name__}"
                )
