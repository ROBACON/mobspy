"""Translate backend-independent execution plans to SBML/COPASI."""

from __future__ import annotations

from mobspy.execution import generate_sbml_from_compiled
from mobspy.sbml.runner import simulate
from mobspy.types import BackendResults, ConcreteModel, ExecutionPlan


class SBMLBackend:
    """Execute compiled model chains using isolated COPASI runs."""

    def generate_model(self, model: ConcreteModel) -> str:
        return generate_sbml_from_compiled(model)

    def run(self, plan: ExecutionPlan, jobs: int = -1) -> BackendResults:
        return [
            simulate(
                jobs,
                list(plan.settings),
                [
                    model.to_compiled_model(
                        species_not_mapped=dict(model.species),
                        mappings=dict(model.mappings),
                    )
                    for model in chain
                ],
            )
            for chain in plan.models
        ]
