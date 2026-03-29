"""Model-wide unit context for MobsPy simulations.

Resolves user-provided units into a consistent model unit system and provides
conversion methods for rates, counts, volumes, and time.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import libsbml as sbml
from pint import Quantity, Unit
from scipy.constants import N_A

from mobspy.exceptions import UnitError
from mobspy.modules.mobspy_expressions import OverrideQuantity
from mobspy.modules.mobspy_expressions import u as _mobspy_u

_ur = _mobspy_u.unit_registry_object


# ---------------------------------------------------------------------------
# Pint unit -> libsbml UnitDefinition mapping
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _SbmlUnitComponent:
    """One <unit> element inside an SBML <unitDefinition>."""

    kind: int
    exponent: int = 1
    scale: int = 0
    multiplier: float = 1.0


# Lookup table: Pint unit string -> list of SBML unit components.
# Composite units (e.g. molar = mol/L) are expressed as multiple components.
_PINT_TO_SBML: dict[str, list[_SbmlUnitComponent]] = {
    # Time
    "second": [_SbmlUnitComponent(sbml.UNIT_KIND_SECOND)],
    "minute": [_SbmlUnitComponent(sbml.UNIT_KIND_SECOND, multiplier=60)],
    "hour": [_SbmlUnitComponent(sbml.UNIT_KIND_SECOND, multiplier=3600)],
    "day": [_SbmlUnitComponent(sbml.UNIT_KIND_SECOND, multiplier=86400)],
    # Substance
    "mole": [_SbmlUnitComponent(sbml.UNIT_KIND_MOLE)],
    "millimole": [_SbmlUnitComponent(sbml.UNIT_KIND_MOLE, scale=-3)],
    "micromole": [_SbmlUnitComponent(sbml.UNIT_KIND_MOLE, scale=-6)],
    "nanomole": [_SbmlUnitComponent(sbml.UNIT_KIND_MOLE, scale=-9)],
    # Volume / length^3
    "liter": [_SbmlUnitComponent(sbml.UNIT_KIND_LITRE)],
    "litre": [_SbmlUnitComponent(sbml.UNIT_KIND_LITRE)],
    "milliliter": [_SbmlUnitComponent(sbml.UNIT_KIND_LITRE, scale=-3)],
    "millilitre": [_SbmlUnitComponent(sbml.UNIT_KIND_LITRE, scale=-3)],
    "microliter": [_SbmlUnitComponent(sbml.UNIT_KIND_LITRE, scale=-6)],
    "microlitre": [_SbmlUnitComponent(sbml.UNIT_KIND_LITRE, scale=-6)],
    # Length (for non-3D or when user specifies length directly)
    "meter": [_SbmlUnitComponent(sbml.UNIT_KIND_METRE)],
    "metre": [_SbmlUnitComponent(sbml.UNIT_KIND_METRE)],
    "decimeter": [_SbmlUnitComponent(sbml.UNIT_KIND_METRE, scale=-1)],
    "decimetre": [_SbmlUnitComponent(sbml.UNIT_KIND_METRE, scale=-1)],
    "centimeter": [_SbmlUnitComponent(sbml.UNIT_KIND_METRE, scale=-2)],
    "centimetre": [_SbmlUnitComponent(sbml.UNIT_KIND_METRE, scale=-2)],
    "millimeter": [_SbmlUnitComponent(sbml.UNIT_KIND_METRE, scale=-3)],
    "millimetre": [_SbmlUnitComponent(sbml.UNIT_KIND_METRE, scale=-3)],
    "micrometer": [_SbmlUnitComponent(sbml.UNIT_KIND_METRE, scale=-6)],
    "micrometre": [_SbmlUnitComponent(sbml.UNIT_KIND_METRE, scale=-6)],
}


# ---------------------------------------------------------------------------
# ModelUnitContext
# ---------------------------------------------------------------------------


@dataclass
class ModelUnitContext:
    """Resolved model-wide unit system for a single compilation.

    When no user units are provided, defaults to second / item / dm^dimension,
    which reproduces the legacy MobsPy behavior identically.
    """

    time_unit: Unit = field(default_factory=lambda: _ur.second)  # type: ignore[assignment]
    substance_unit: Unit | None = field(default=None)  # None = counts (item)
    volume_unit: Unit = field(default_factory=lambda: _ur.decimeter**3)  # type: ignore[assignment]
    dimension: int = 3

    # ---- internal cache set by from_simulation ----
    _time_set: bool = field(default=False, repr=False)
    _substance_set: bool = field(default=False, repr=False)
    _volume_set: bool = field(default=False, repr=False)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _to_plain_quantity(q: Any) -> Quantity:  # type: ignore[type-arg]
        """Convert OverrideQuantity to a plain Pint Quantity for safe .to() calls."""
        if isinstance(q, OverrideQuantity):
            return q.q_object  # type: ignore[no-any-return]
        return q  # type: ignore[no-any-return]

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def substance_is_molar(self) -> bool:
        """True when the model tracks substance in moles (not counts)."""
        return self.substance_unit is not None

    @property
    def volume_unit_str(self) -> str:
        return str(self.volume_unit)

    @property
    def time_unit_str(self) -> str:
        return str(self.time_unit)

    @property
    def substance_unit_str(self) -> str:
        if self.substance_unit is None:
            return "item"
        return str(self.substance_unit)

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @classmethod
    def from_simulation(
        cls,
        volume: Any = 1,
        duration: Any = 60,
        species_counts: list[dict[str, Any]] | None = None,
        dimension: int | None = None,
    ) -> ModelUnitContext:
        """Scan user-provided inputs and resolve the model unit system.

        Implements first-unit-wins with volume/duration priority.

        Non-default volume/time units are only activated when species counts
        include substance units (e.g. moles, molar). This prevents changing
        the numerical behavior of models where volume is provided with units
        but rates/counts are dimensionless.
        """
        ctx = cls()

        # Track volume/time units from input, but don't commit yet
        _vol_unit_from_input = None
        _time_unit_from_input = None

        # 1. Volume (extract dimension, tentatively store volume unit)
        if isinstance(volume, Quantity):
            vol_q = cls._to_plain_quantity(volume)
            vol_ur = _ur.Quantity(1, str(vol_q.units))
            dim_dict = dict(vol_ur.dimensionality)
            if "[length]" in dim_dict:
                length_exp = int(dim_dict["[length]"])
                if dimension is None:
                    dimension = length_exp
                ctx.dimension = dimension
                _vol_unit_from_input = vol_ur.units

        if dimension is not None:
            ctx.dimension = dimension

        # 2. Duration (tentatively store time unit)
        if isinstance(duration, Quantity):
            dur_q = cls._to_plain_quantity(duration)
            dur_ur = _ur.Quantity(1, str(dur_q.units))
            dim_dict = dict(dur_ur.dimensionality)
            if "[time]" in dim_dict and len(dim_dict) == 1:
                _time_unit_from_input = dur_ur.units

        # 3. Species counts (substance detection + volume from concentrations)
        has_substance_units = False
        if species_counts is not None:
            for count in species_counts:
                q = count.get("quantity")
                if not isinstance(q, Quantity):
                    continue
                dim_dict = dict(q.dimensionality)

                if "[substance]" in dim_dict:
                    has_substance_units = True
                    if not ctx._substance_set:
                        ctx._substance_set = True
                        ctx.substance_unit = _extract_substance_unit(q)

                    if "[length]" in dim_dict and not ctx._volume_set:
                        ctx._volume_set = True
                        ctx.volume_unit = _extract_volume_unit(q, ctx.dimension)

        # Update default volume unit to match resolved dimension
        if not ctx._volume_set:
            ctx.volume_unit = _ur.decimeter**ctx.dimension  # type: ignore[assignment]

        # Only activate non-default volume/time units when substance units
        # are present. This ensures backward compatibility for models with
        # unit volumes but dimensionless rates/counts.
        if has_substance_units:
            if _vol_unit_from_input is not None and not ctx._volume_set:
                ctx.volume_unit = _vol_unit_from_input  # type: ignore[assignment]
                ctx._volume_set = True
            if _time_unit_from_input is not None and not ctx._time_set:
                ctx.time_unit = _time_unit_from_input  # type: ignore[assignment]
                ctx._time_set = True

        return ctx

    # ------------------------------------------------------------------
    # Conversion methods
    # ------------------------------------------------------------------

    def convert_rate(
        self,
        quantity: int | float | Quantity,
        reaction_order: int,
        dimension: int | None = None,
    ) -> tuple[float | int | Any, int, bool]:
        """Convert a rate constant to model units.

        Returns (magnitude, dimension, is_count) matching the signature
        of unit_handler.convert_rate().
        """
        if dimension is None:
            dimension = self.dimension

        if not isinstance(quantity, Quantity):
            return quantity, dimension, False

        quantity = self._to_plain_quantity(quantity)
        volume_power = reaction_order - 1
        dim_dict = dict(quantity.dimensionality)
        has_time = "[time]" in dim_dict
        has_substance = "[substance]" in dim_dict
        has_length = "[length]" in dim_dict

        try:
            if has_time and not has_substance and not has_length:
                # Pure 1/[time] rate (count-based, order 0 or 1)
                target = 1 / self.time_unit
                converted = quantity.to(target)
                return converted.magnitude, dimension, True

            elif has_substance and not has_length:
                # [substance]/[time] rate (e.g. mol/s)
                if self.substance_is_molar:
                    assert self.substance_unit is not None
                    target = self.substance_unit / self.time_unit
                    converted = quantity.to(target)
                    return converted.magnitude, dimension, True
                else:
                    # Convert to moles/time_unit then multiply by N_A for counts
                    target = _ur.mole / self.time_unit
                    converted = quantity.to(target)
                    return converted.magnitude * N_A, dimension, True

            elif has_substance:
                # Concentration rate with moles: [length]^n/([substance]^m*[time])
                if self.substance_is_molar:
                    assert self.substance_unit is not None
                    target = self.volume_unit**volume_power / (
                        self.substance_unit**volume_power * self.time_unit
                    )
                    converted = quantity.to(target)
                    return converted.magnitude, dimension, False
                else:
                    # Convert substance to moles then to counts via N_A
                    target = self.volume_unit**volume_power / (
                        _ur.mole**volume_power * self.time_unit
                    )
                    converted = quantity.to(target)
                    return (
                        converted.magnitude / (N_A**volume_power),
                        dimension,
                        False,
                    )

            else:
                # [length]^n/[time] concentration rate (no substance)
                target = self.volume_unit**volume_power / self.time_unit
                converted = quantity.to(target)
                return converted.magnitude, dimension, False

        except Exception as e:
            raise UnitError(
                str(e) + "\n"
                f"Problem converting rate {quantity}\n"
                f"Is the rate in the form [volume]**{volume_power}/[time]?"
            ) from e

    def convert_counts(
        self,
        quantity: int | float | Quantity | Any,
        volume: int | float,
    ) -> Any:
        """Convert species count/concentration to model units."""
        if not isinstance(quantity, Quantity):
            return quantity

        quantity = self._to_plain_quantity(quantity)
        dim_dict = dict(quantity.dimensionality)
        has_length = "[length]" in dim_dict
        has_substance = "[substance]" in dim_dict
        is_dimensionless = len(dim_dict) == 0

        if not has_length and not has_substance and not is_dimensionless:
            raise UnitError(
                f"The assigned quantity {quantity} is neither a count or concentration"
            )

        if is_dimensionless:
            return quantity.magnitude

        try:
            if has_substance:
                if has_length:
                    # Concentration (e.g. molar, millimolar)
                    if self.substance_is_molar:
                        assert self.substance_unit is not None
                        target = self.substance_unit / self.volume_unit
                        converted = quantity.to(target)
                        # Multiply by volume to get amount
                        return converted.magnitude * volume
                    else:
                        # Convert to moles/volume_unit then to counts
                        target = _ur.mole / self.volume_unit
                        converted = quantity.to(target)
                        return converted.magnitude * volume * N_A
                # Pure substance amount (e.g. 5 * u.mole)
                elif self.substance_is_molar:
                    assert self.substance_unit is not None
                    converted = quantity.to(self.substance_unit)
                    return converted.magnitude
                else:
                    converted = quantity.to(_ur.mole)
                    return converted.magnitude * N_A
            else:
                if has_length:
                    # Count concentration (e.g. 1/L)
                    target = 1 / self.volume_unit
                    converted = quantity.to(target)
                    return converted.magnitude * volume
                return quantity.magnitude

        except Exception as e:
            raise UnitError(
                str(e) + "\n"
                f"Problem converting count {quantity}\n"
                "Is it really a count or concentration?"
            ) from e

    def convert_volume(self, volume: int | float | Quantity) -> int | float:  # type: ignore[type-arg]
        """Convert a volume quantity to model volume units."""
        if isinstance(volume, Quantity):
            pq = self._to_plain_quantity(volume)
            return pq.to(self.volume_unit).magnitude  # type: ignore[no-any-return]
        return volume  # pyright: ignore[reportReturnType]

    def convert_time(self, time: int | float | Quantity) -> int | float:  # type: ignore[type-arg]
        """Convert a time quantity to model time units."""
        if isinstance(time, Quantity):
            pq = self._to_plain_quantity(time)
            dim = dict(pq.dimensionality)
            if dim.get("[time]") and len(dim) == 1:
                return pq.to(self.time_unit).magnitude  # type: ignore[no-any-return]
        return time  # type: ignore[return-value]

    # ------------------------------------------------------------------
    # SBML unit definition generation
    # ------------------------------------------------------------------

    def get_sbml_time_units_id(self) -> str:
        """SBML id string for the model's time unit."""
        s = str(self.time_unit)
        # Normalize common names
        return _sanitize_sbml_id(s)

    def get_sbml_substance_units_id(self) -> str:
        """SBML id string for the model's substance unit."""
        if self.substance_unit is None:
            return "item"
        return _sanitize_sbml_id(str(self.substance_unit))

    def get_sbml_volume_units_id(self) -> str:
        """SBML id string for the model's volume unit."""
        return _sanitize_sbml_id(str(self.volume_unit))

    def create_sbml_unit_definitions(self, model: Any) -> None:
        """Create UnitDefinition elements on an SBML model object."""
        # Time unit
        time_id = self.get_sbml_time_units_id()
        if time_id != "second":
            _create_unit_def(model, time_id, self.time_unit)

        # Substance unit
        if self.substance_is_molar:
            assert self.substance_unit is not None
            sub_id = self.get_sbml_substance_units_id()
            _create_unit_def(model, sub_id, self.substance_unit)

        # Volume unit
        vol_id = self.get_sbml_volume_units_id()
        if vol_id != "dimensionless":
            _create_unit_def(model, vol_id, self.volume_unit)

        # Rate unit (per model time) - always useful
        rate_id = f"per_{time_id}"
        _create_rate_unit_def(model, rate_id, self.time_unit)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _extract_substance_unit(quantity: Quantity) -> Unit:  # type: ignore[type-arg]
    """Extract the substance component from a Quantity's units.

    For example, from millimolar (mmol/L), extract millimole.
    From mol/s, extract mole.
    """
    # Ensure we work with the shared registry
    quantity = ModelUnitContext._to_plain_quantity(quantity)
    # Re-create in _ur to avoid cross-registry issues
    q = _ur.Quantity(1, str(quantity.units))
    dim = dict(q.dimensionality)
    if "[substance]" not in dim:
        msg = f"No substance dimension in {quantity}"
        raise UnitError(msg)

    # Build complementary quantity to neutralize non-substance dimensions
    complement = _ur.Quantity(1.0, "dimensionless")
    for dim_name, power in dim.items():
        if dim_name == "[substance]":
            continue
        if dim_name == "[length]":
            complement = complement * _ur.liter ** (-int(power) // 3)
        elif dim_name == "[time]":
            complement = complement * _ur.second ** (-int(power))
        elif dim_name == "[mass]":
            complement = complement * _ur.kilogram ** (-int(power))

    isolated = q * complement

    # Try to match against standard substance units
    candidates: list[Any] = [
        _ur.nanomole,
        _ur.micromole,
        _ur.millimole,
        _ur.mole,
    ]
    for target in candidates:
        try:
            converted = isolated.to(target)
            if 0.99 < abs(converted.magnitude) < 1.01:
                return target  # type: ignore[no-any-return]
        except Exception:
            continue

    return _ur.mole  # type: ignore[no-any-return]


def _extract_volume_unit(quantity: Quantity, dimension: int) -> Unit:  # type: ignore[type-arg]
    """Extract the volume component from a Quantity.

    For example, from millimolar (mmol/L), extract liter.
    From 5 * u.mL, extract milliliter.
    """
    # Ensure we work with the shared registry
    quantity = ModelUnitContext._to_plain_quantity(quantity)
    q = _ur.Quantity(1, str(quantity.units))
    dim = dict(q.dimensionality)
    length_power = int(dim.get("[length]", 0))
    if length_power == 0:
        msg = f"No length dimension in {quantity}"
        raise UnitError(msg)

    abs_power = abs(length_power)

    # For pure volume quantities (e.g. 5 * u.mL), try direct matching
    if len(dim) == 1:
        q_norm = q if length_power > 0 else (1 / q)
        candidates: list[Any] = [_ur.microliter, _ur.milliliter, _ur.liter]
        for target in candidates:
            try:
                converted = q_norm.to(target)
                if 0.99 < abs(converted.magnitude) < 1.01:
                    return target  # type: ignore[no-any-return]
            except Exception:
                continue

    # For compound units (e.g. millimolar = mmol/L), convert to base units
    q.to_base_units()

    if abs_power == dimension:
        # It's a volume: figure out which standard volume from the base magnitude
        # Extract just the length contribution by checking base conversion
        # against known volume units
        # Use a reference: 1 liter = 0.001 m^3
        candidates_vol: list[tuple[Any, float]] = [
            (_ur.microliter, 1e-9),  # m^3
            (_ur.milliliter, 1e-6),  # m^3
            (_ur.liter, 1e-3),  # m^3
        ]
        # The base quantity has the volume contribution embedded.
        # For molar: base = 1e3 mol/m^3, volume part is m^-3 -> 1/m^3 -> 1000 L
        # We can't easily extract it from compound units.
        # Simpler: for concentrations (substance + length), just use liter as default.
        if "[substance]" in dim:
            return _ur.liter  # type: ignore[no-any-return]

        # For pure volume quantities, match by converting
        for target, _base_m3 in candidates_vol:
            try:
                test = q.to(target)
                if 0.99 < abs(test.magnitude) < 1.01:
                    return target  # type: ignore[no-any-return]
            except Exception:
                continue

    return _ur.decimeter**dimension  # type: ignore[no-any-return, return-value]


def _sanitize_sbml_id(unit_str: str) -> str:
    """Convert a Pint unit string to a valid SBML unit definition ID."""
    # Replace spaces and special chars
    result = unit_str.replace(" ", "_").replace("**", "").replace("^", "")
    result = result.replace("/", "_per_").replace("*", "_")
    # Remove leading/trailing underscores
    result = result.strip("_")
    # SBML IDs must start with a letter
    if result and not result[0].isalpha():
        result = "unit_" + result
    return result


def _create_unit_def(model: Any, unit_id: str, pint_unit: Unit) -> None:
    """Create an SBML UnitDefinition from a Pint unit, using the lookup table."""
    unit_str = str(pint_unit)

    # Try direct lookup first
    components = _PINT_TO_SBML.get(unit_str)

    if components is None:
        # Try to decompose: convert to base units and build from there
        components = _decompose_pint_unit(pint_unit)

    if components is None:
        # Fallback: dimensionless
        components = [_SbmlUnitComponent(sbml.UNIT_KIND_DIMENSIONLESS)]

    ud = model.createUnitDefinition()
    ud.setId(unit_id)
    for comp in components:
        u = ud.createUnit()
        u.setKind(comp.kind)
        u.setExponent(comp.exponent)
        u.setScale(comp.scale)
        u.setMultiplier(comp.multiplier)


def _create_rate_unit_def(model: Any, rate_id: str, time_unit: Unit) -> None:
    """Create a per-time-unit definition for rates."""
    time_str = str(time_unit)
    time_components = _PINT_TO_SBML.get(time_str)

    if time_components is None:
        time_components = _decompose_pint_unit(time_unit)
    if time_components is None:
        time_components = [_SbmlUnitComponent(sbml.UNIT_KIND_SECOND)]

    ud = model.createUnitDefinition()
    ud.setId(rate_id)
    for comp in time_components:
        u = ud.createUnit()
        u.setKind(comp.kind)
        u.setExponent(-comp.exponent)  # Invert for "per"
        u.setScale(comp.scale)
        u.setMultiplier(comp.multiplier)


def _decompose_pint_unit(pint_unit: Unit) -> list[_SbmlUnitComponent] | None:
    """Decompose a Pint unit into SBML unit components via base units."""
    try:
        # Get the quantity in base SI units
        q = (1 * pint_unit).to_base_units()  # pyright: ignore[reportAttributeAccessIssue]
        magnitude = q.magnitude
        dim = dict(q.dimensionality)
    except Exception:
        return None

    components: list[_SbmlUnitComponent] = []

    # Map Pint base dimensions to SBML unit kinds
    dim_to_kind = {
        "[time]": sbml.UNIT_KIND_SECOND,
        "[length]": sbml.UNIT_KIND_METRE,
        "[substance]": sbml.UNIT_KIND_MOLE,
        "[mass]": sbml.UNIT_KIND_KILOGRAM,
        "[temperature]": sbml.UNIT_KIND_KELVIN,
    }

    first = True
    for dim_name, power in dim.items():
        kind = dim_to_kind.get(dim_name)
        if kind is None:
            continue
        comp = _SbmlUnitComponent(
            kind=kind,
            exponent=int(power),
            scale=0,
            # Put the magnitude multiplier on the first component only
            multiplier=magnitude if first else 1.0,
        )
        components.append(comp)
        if first:
            first = False

    if not components:
        return None

    return components
