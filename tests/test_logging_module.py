"""Tests for the logging module."""

from __future__ import annotations

import logging

import pytest

from mobspy.exceptions import MobsPyError
from mobspy.mobspy_logging import get_logger

log_scripts = get_logger("test_logger")


def test_debug_function_exists():
    log_scripts.set_log_level(logging.DEBUG)
    log_scripts.debug("Test debug message")


def test_info_function_exists():
    log_scripts.set_log_level(logging.INFO)
    log_scripts.info("Test info message")


def test_warning_function_exists():
    log_scripts.set_log_level(logging.WARNING)
    log_scripts.warning("Test warning message")


def test_error_raises_mobspy_error():
    with pytest.raises(MobsPyError, match="Test error message"):
        log_scripts.error("Test error message")


def test_set_log_level_with_string():
    log_scripts.set_log_level("DEBUG")
    log_scripts.set_log_level("INFO")
    log_scripts.set_log_level("WARNING")
    log_scripts.set_log_level("ERROR")
    log_scripts.set_log_level("CRITICAL")


def test_set_log_level_with_constant():
    log_scripts.set_log_level(logging.DEBUG)
    log_scripts.set_log_level(logging.INFO)
    log_scripts.set_log_level(logging.WARNING)
    log_scripts.set_log_level(logging.ERROR)
    log_scripts.set_log_level(logging.CRITICAL)


def test_set_log_level_with_integer():
    log_scripts.set_log_level(10)
    log_scripts.set_log_level(20)
    log_scripts.set_log_level(30)
    log_scripts.set_log_level(40)
    log_scripts.set_log_level(50)


def test_error_with_full_exception_log():
    with pytest.raises(MobsPyError):
        log_scripts.error("Error with trace", full_exception_log=True)


def test_log_level_affects_output():
    log_scripts.set_log_level(logging.ERROR)
    log_scripts.debug("Should be filtered")
    log_scripts.info("Should be filtered")
    log_scripts.warning("Should be filtered")
    log_scripts.set_log_level(logging.DEBUG)
    log_scripts.debug("Should work")
    log_scripts.info("Should work")
    log_scripts.warning("Should work")


def test_colored_formatter_exists():
    import mobspy.mobspy_logging

    assert hasattr(mobspy.mobspy_logging, "ColoredFormatter")
