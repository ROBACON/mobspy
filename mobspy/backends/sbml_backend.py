"""SBML/COPASI simulation backend.

Default backend for MobsPy simulations. Generates SBML XML
and runs simulations via BasiCO/COPASI.
"""

from __future__ import annotations

from typing import Any

from mobspy.execution import generate_sbml_from_compiled, run_sbml
from mobspy.types import ConcreteModel


class SBMLBackend:
    """SBML/COPASI simulation backend.

    Wraps the existing standalone functions ``generate_sbml_from_compiled``
    and ``run_sbml`` behind the ``SimulationBackend`` protocol.
    """

    def generate_model(
        self,
        model: ConcreteModel,
        model_context: Any = None,  # noqa: ARG002
    ) -> str:
        """Generate SBML XML from compiled model data.

        Args:
            model: A ConcreteModel instance.
            model_context: Unused. The ConcreteModel carries its own
                unit_context. Kept for protocol compatibility.
        """
        return generate_sbml_from_compiled(
            model.to_compiler_result(), model.unit_context
        )

    def run(
        self,
        model_strings: list[Any],
        parameters: list[Any],
        jobs: int = -1,
    ) -> list[Any]:
        """Execute simulations via BasiCO/COPASI."""
        return run_sbml(model_strings, parameters, jobs)
