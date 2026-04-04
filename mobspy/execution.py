"""Standalone simulation execution and SBML generation functions.

These functions provide a composition-friendly alternative to the
Simulation class methods.  They can be used independently of the
Simulation object for programmatic workflows.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from mobspy.sbml_simulator.builder import build as sbml_build
from mobspy.sbml_simulator.run import simulate as sbml_simulate
from mobspy.types import SBMLModelData

if TYPE_CHECKING:
    from mobspy.types import CompilerResult


def generate_sbml_from_compiled(
    compiled: CompilerResult,
    model_context: Any = None,
) -> str:
    """Generate an SBML string from a CompilerResult.

    Standalone function that does not require a Simulation instance.

    Args:
        compiled: Result from ``compile_model()``.
        model_context: Optional ModelUnitContext for unit-aware output.

    Returns:
        SBML XML string.
    """
    sbml_data = SBMLModelData(
        species_for_sbml=compiled.species_for_sbml,
        parameters_for_sbml=compiled.parameters_for_sbml,
        reactions_for_sbml=compiled.reactions_for_sbml,
        events_for_sbml=compiled.events_for_sbml,
        assignments_for_sbml=compiled.assignments_for_sbml,
    )
    return sbml_build(sbml_data, model_context=model_context)


def run_sbml(
    sbml_models: list[list[Any]],
    parameters: list[Any],
    jobs: int = -1,
) -> list[Any]:
    """Run SBML simulations via BasiCO/COPASI.

    Standalone wrapper around the simulation engine.

    Args:
        sbml_models: Nested list of SBMLModelData (parameter sweeps).
        parameters: List of simulation parameter dicts.
        jobs: Number of parallel jobs (-1 for auto).

    Returns:
        List of raw simulation results.
    """
    import joblib  # noqa: PLC0415

    def _sim_one(x: Any) -> Any:
        return sbml_simulate(jobs, parameters, x)

    results: list[Any] = list(
        joblib.Parallel(n_jobs=jobs, prefer="threads")(
            joblib.delayed(_sim_one)(sbml) for sbml in sbml_models
        )
    )
    return results
