"""Compatibility views derived from one authoritative compiled model."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from mobspy.constants import DOT_SEPARATOR

if TYPE_CHECKING:
    from mobspy.types import (
        AssignmentsForSbml,
        CompiledModel,
        ConcreteModel,
        EventsForSbml,
        MappingsForSbml,
        ParametersForSbml,
        ParameterUsedInfo,
        ReactionsForSbml,
        SpeciesForSbml,
    )
    from mobspy.units.model_context import ModelUnitContext


class CompiledState:
    """Read-through legacy attributes; none owns a second model copy."""

    _concrete_model: ConcreteModel | None

    @property
    def _species_for_sbml(self) -> SpeciesForSbml | None:
        return (
            self._concrete_model.species if self._concrete_model is not None else None
        )

    @property
    def _reactions_for_sbml(self) -> ReactionsForSbml | None:
        return (
            self._concrete_model.reactions if self._concrete_model is not None else None
        )

    @property
    def _parameters_for_sbml(self) -> ParametersForSbml | None:
        return (
            self._concrete_model.parameters
            if self._concrete_model is not None
            else None
        )

    @property
    def _mappings_for_sbml(self) -> MappingsForSbml | None:
        return (
            self._concrete_model.mappings if self._concrete_model is not None else None
        )

    @property
    def _events_for_sbml(self) -> EventsForSbml | None:
        return self._concrete_model.events if self._concrete_model is not None else None

    @property
    def _assignments_for_sbml(self) -> AssignmentsForSbml:
        return (
            self._concrete_model.assignments if self._concrete_model is not None else {}
        )

    @property
    def _model_context(self) -> ModelUnitContext | None:
        return (
            self._concrete_model.unit_context
            if self._concrete_model is not None
            else None
        )

    @property
    def model_string(self) -> str:
        return (
            self._concrete_model.model_string
            if self._concrete_model is not None
            else ""
        )

    @property
    def mappings(self) -> MappingsForSbml:
        return self._concrete_model.mappings if self._concrete_model is not None else {}

    @property
    def model_parameters(self) -> dict[str, ParameterUsedInfo]:
        return (
            self._concrete_model.parameters_used
            if self._concrete_model is not None
            else {}
        )

    @property
    def model_parameter_objects_dict(self) -> dict[str, Any] | None:
        return (
            self._concrete_model.parameter_objects
            if self._concrete_model is not None
            else None
        )

    @property
    def _has_mole(self) -> bool:
        return (
            self._concrete_model.has_mole if self._concrete_model is not None else False
        )

    @property
    def _assigned_species_list(self) -> list[str]:
        return (
            list(self._concrete_model.assigned_species) if self._concrete_model else []
        )

    @property
    def all_species_not_mapped(self) -> dict[str, int | float]:
        return {
            name.replace(DOT_SEPARATOR, "."): value
            for name, value in (self._species_for_sbml or {}).items()
        }

    @property
    def _list_of_models(self) -> list[CompiledModel]:
        if self._concrete_model is None:
            return []
        return [
            self._concrete_model.to_compiled_model(
                species_not_mapped=self.all_species_not_mapped,
                mappings=self.mappings,
            )
        ]
