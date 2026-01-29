import logging
import sys
import traceback

# Configure the logger for MobsPy
logger = logging.getLogger("mobspy")
logger.setLevel(logging.INFO)  # Logger accepts all levels

# Create console handler with formatting
_console_handler = logging.StreamHandler(sys.stderr)
_console_handler.setLevel(logging.DEBUG)  # Default: show debug and above


class ColoredFormatter(logging.Formatter):
    """Custom formatter that adds color to warning and error messages"""

    COLORS: dict[str, str] = {
        "WARNING": "\033[93m",  # Yellow
        "ERROR": "\033[91m",  # Red
        "CRITICAL": "\033[91m",  # Red
        "RESET": "\033[0m",
    }

    def format(self, record: logging.LogRecord) -> str:
        log_message = super().format(record)
        if record.levelname in self.COLORS:
            return f"{self.COLORS[record.levelname]}{log_message}{self.COLORS['RESET']}"
        return log_message


_formatter = ColoredFormatter("%(levelname)s: %(message)s")
_console_handler.setFormatter(_formatter)
logger.addHandler(_console_handler)

# Prevent propagation to root logger
logger.propagate = False

# For debugging - enable full exception logging
_PRINT_FULL_EXCEPTION_LOG: bool = False


def set_log_level(level: int | str) -> None:
    """
    Set the logging level.

    Args:
        level: Either a logging level constant (logging.DEBUG, logging.INFO, etc.)
               or an integer (DEBUG=10, INFO=20, WARNING=30, ERROR=40, CRITICAL=50)
               or a string ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL')
    """
    if isinstance(level, str):
        level = getattr(logging, level.upper())
    _console_handler.setLevel(level)


def error(message: str, full_exception_log: bool = False) -> None:
    """
    Log an error message and exit the program.

    Args:
        message: The error message to log
        full_exception_log: If True, include full traceback
    """
    # Include full traceback if requested
    if _PRINT_FULL_EXCEPTION_LOG or full_exception_log:
        traceback_details: str = "".join(traceback.format_stack())
        message = f"Full Exception Log:\n{traceback_details}\n{message}"

    logger.error(message, stacklevel=2)
    sys.exit(1)


def debug(message: str = "") -> None:
    """
    Log a debug message.

    Args:
        message: The debug message to log
    """
    # Replace _dot_ with . for better readability
    formatted_message: str = str(message).replace("_dot_", ".")
    logger.debug(formatted_message, stacklevel=2)


def info(message: str) -> None:
    """
    Log an info message.

    Args:
        message: The info message to log
    """
    # Replace _dot_ with . for better readability
    formatted_message: str = str(message).replace("_dot_", ".")
    logger.info(formatted_message, stacklevel=2)


def warning(message: str) -> None:
    """
    Log a warning message.

    Args:
        message: The warning message to log
    """
    # Replace _dot_ with . for better readability
    formatted_message: str = str(message).replace("_dot_", ".")
    logger.warning(formatted_message, stacklevel=2)


# Backward compatibility: map old global_simlog_level to logging levels
# Old system: 0=ERROR, 1=WARNING, 2+=DEBUG
# Note: This is a module attribute, setting it won't automatically update the logger
# Use set_log_level() instead for dynamic changes
# global_simlog_level: int = 1  # Default level
