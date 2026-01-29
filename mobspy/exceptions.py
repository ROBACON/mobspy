"""
Custom exceptions for MobsPy to provide better error handling and debugging.
"""


class MobsPyError(Exception):
    """Base exception for all MobsPy-related errors."""

    pass


class CompilationError(MobsPyError):
    """Errors that occur during model compilation."""

    pass


class SimulationError(MobsPyError):
    """Errors that occur during simulation execution."""

    pass


class ParameterError(MobsPyError):
    """Errors related to invalid or missing parameters."""

    pass


class ReactionError(MobsPyError):
    """Errors related to reaction definition and processing."""

    pass


class EventError(MobsPyError):
    """Errors related to event handling and processing."""

    pass


class ValidationError(MobsPyError):
    """Errors related to input validation and model constraints."""

    pass


class ImportError(MobsPyError):
    """Errors related to module imports and dependencies."""

    pass


class UnitError(MobsPyError):
    """Errors related to unit handling and conversion."""

    pass


class SBMLError(MobsPyError):
    """Errors related to SBML generation and processing."""

    pass


class AntimonyError(MobsPyError):
    """Errors related to Antimony model generation."""

    pass
