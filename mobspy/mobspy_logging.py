"""
Centralized logging system for MobsPy.

This module provides a consistent logging interface throughout the MobsPy codebase,
allowing for better control over log levels, formatting, and output destinations.
"""

import logging
import sys
import traceback
from typing import Optional, Any


class ColoredFormatter(logging.Formatter):
    """Custom formatter that adds color to warning and error messages."""

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


class MobsPyLogger:
    """
    Centralized logger for MobsPy with configurable log levels and formatting.

    This class provides a consistent logging interface that can be used throughout
    the MobsPy codebase, making it easier to control logging behavior and output.
    """

    def __init__(self, name: str = "mobspy"):
        """
        Initialize the MobsPy logger.

        Args:
            name: Name of the logger (default: "mobspy")
        """
        self.logger = logging.getLogger(name)
        self.logger.setLevel(
            logging.DEBUG
        )  # Set to lowest level, filtering happens at handlers

        # Prevent duplicate handlers if logger already exists
        if not self.logger.handlers:
            self._configure_default_logging()

    def _configure_default_logging(self) -> None:
        """Set up default logging configuration with console handler."""
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(logging.DEBUG)

        formatter = ColoredFormatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        console_handler.setFormatter(formatter)

        self.logger.addHandler(console_handler)
        self.logger.propagate = False

    def _format_message(self, message: str) -> str:
        """
        Format message by replacing _dot_ with . for readability.

        Args:
            message: The message to format

        Returns:
            Formatted message
        """
        return str(message).replace("_dot_", ".")

    def set_log_level(self, level: int | str) -> None:
        """
        Set the logging level for all handlers.

        Args:
            level: Logging level constant (logging.DEBUG, logging.INFO, etc.),
                   integer (10, 20, 30, 40, 50), or string ('DEBUG', 'INFO', etc.)
        """
        if isinstance(level, str):
            level = getattr(logging, level.upper())
        for handler in self.logger.handlers:
            handler.setLevel(level)

    def add_file_handler(self, filename: str, level: int = logging.DEBUG) -> None:
        """
        Add a file handler to the logger.

        Args:
            filename: Path to the log file
            level: Logging level for the file handler
        """
        file_handler = logging.FileHandler(filename)
        file_handler.setLevel(level)

        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        file_handler.setFormatter(formatter)

        self.logger.addHandler(file_handler)

    def debug(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Log a debug message."""
        formatted_message = self._format_message(message)
        self.logger.debug(formatted_message, *args, stacklevel=2, **kwargs)

    def info(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Log an info message."""
        formatted_message = self._format_message(message)
        self.logger.info(formatted_message, *args, stacklevel=2, **kwargs)

    def warning(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Log a warning message."""
        formatted_message = self._format_message(message)
        self.logger.warning(formatted_message, *args, stacklevel=2, **kwargs)

    def error(
        self, message: str, *args: Any, full_exception_log: bool = False, **kwargs: Any
    ) -> None:
        """
        Log an error message and exit the program.

        Args:
            message: The error message to log
            full_exception_log: If True, include full traceback
        """
        formatted_message = self._format_message(message)

        # Include full traceback if requested
        if full_exception_log:
            traceback_details: str = "".join(traceback.format_stack())
            formatted_message = (
                f"Full Exception Log:\n{traceback_details}\n{formatted_message}"
            )

        self.logger.error(formatted_message, *args, stacklevel=2, **kwargs)
        sys.exit(1)

    def critical(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Log a critical message."""
        formatted_message = self._format_message(message)
        self.logger.critical(formatted_message, *args, stacklevel=2, **kwargs)

    def exception(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Log an exception with traceback."""
        formatted_message = self._format_message(message)
        self.logger.exception(formatted_message, *args, stacklevel=2, **kwargs)


# Global logger instance
logger = MobsPyLogger()


def get_logger(name: Optional[str] = None) -> MobsPyLogger:
    """
    Get a logger instance.

    Args:
        name: Optional name for the logger. If None, returns the default logger.

    Returns:
        MobsPyLogger instance
    """
    if name is None:
        return logger
    return MobsPyLogger(name)
