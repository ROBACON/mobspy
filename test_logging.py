"""
Tests for the logging module (mobspy/mobspy_logging.py)
"""

import logging
import pytest
from mobspy.mobspy_logging import get_logger

# Get a test logger instance
log_scripts = get_logger("test_logger")


def test_debug_function_exists():
    """Test debug() function exists and can be called"""
    log_scripts.set_log_level(logging.DEBUG)
    log_scripts.debug("Test debug message")
    # If no exception, test passes


def test_info_function_exists():
    """Test info() function exists and can be called"""
    log_scripts.set_log_level(logging.INFO)
    log_scripts.info("Test info message")
    # If no exception, test passes


def test_warning_function_exists():
    """Test warning() function exists and can be called"""
    log_scripts.set_log_level(logging.WARNING)
    log_scripts.warning("Test warning message")
    # If no exception, test passes


def test_error_function_exits():
    """Test error() function exits the program with code 1"""
    with pytest.raises(SystemExit) as exc_info:
        log_scripts.error("Test error message")
    assert exc_info.value.code == 1


def test_set_log_level_with_string():
    """Test set_log_level() accepts string argument"""
    log_scripts.set_log_level("DEBUG")
    log_scripts.set_log_level("INFO")
    log_scripts.set_log_level("WARNING")
    log_scripts.set_log_level("ERROR")
    log_scripts.set_log_level("CRITICAL")
    # If no exception, test passes


def test_set_log_level_with_constant():
    """Test set_log_level() accepts logging constants"""
    log_scripts.set_log_level(logging.DEBUG)
    log_scripts.set_log_level(logging.INFO)
    log_scripts.set_log_level(logging.WARNING)
    log_scripts.set_log_level(logging.ERROR)
    log_scripts.set_log_level(logging.CRITICAL)
    # If no exception, test passes


def test_set_log_level_with_integer():
    """Test set_log_level() accepts integer argument"""
    log_scripts.set_log_level(10)  # DEBUG
    log_scripts.set_log_level(20)  # INFO
    log_scripts.set_log_level(30)  # WARNING
    log_scripts.set_log_level(40)  # ERROR
    log_scripts.set_log_level(50)  # CRITICAL
    # If no exception, test passes


def test_error_with_full_exception_log():
    """Test error() with full_exception_log=True"""
    with pytest.raises(SystemExit):
        log_scripts.error("Error with trace", full_exception_log=True)
    # Test passes if SystemExit is raised



def test_log_level_affects_output():
    """Test that changing log level affects which messages are shown"""
    # Set to ERROR level
    log_scripts.set_log_level(logging.ERROR)
    # These should be silently ignored (no exception)
    log_scripts.debug("Should be filtered")
    log_scripts.info("Should be filtered")
    log_scripts.warning("Should be filtered")

    # Reset to DEBUG
    log_scripts.set_log_level(logging.DEBUG)
    # These should all work
    log_scripts.debug("Should work")
    log_scripts.info("Should work")
    log_scripts.warning("Should work")
    # If no exception, test passes


def test_colored_formatter_exists():
    """Test that ColoredFormatter class exists"""
    import mobspy.mobspy_logging
    assert hasattr(mobspy.mobspy_logging, 'ColoredFormatter')


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
