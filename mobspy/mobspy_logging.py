"""
Centralized logging system for MobsPy.

This module provides a consistent logging interface throughout the MobsPy codebase,
allowing for better control over log levels, formatting, and output destinations.
"""

import logging
import sys
from typing import Optional, Dict, Any


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
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)

        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        console_handler.setFormatter(formatter)

        self.logger.addHandler(console_handler)

    def set_log_level(self, level: int) -> None:
        """
        Set the logging level for all handlers.

        Args:
            level: Logging level (e.g., logging.DEBUG, logging.INFO, logging.WARNING)
        """
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
        self.logger.debug(message, *args, **kwargs)

    def info(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Log an info message."""
        self.logger.info(message, *args, **kwargs)

    def warning(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Log a warning message."""
        self.logger.warning(message, *args, **kwargs)

    def error(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Log an error message."""
        self.logger.error(message, *args, **kwargs)

    def critical(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Log a critical message."""
        self.logger.critical(message, *args, **kwargs)

    def exception(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Log an exception with traceback."""
        self.logger.exception(message, *args, **kwargs)


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
