"""
This module is responsible for converting a model_str into a SBML format
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import libsbml as sbml

from mobspy.exceptions import SBMLError

if TYPE_CHECKING:
    from mobspy.types import (
        AssignmentsForSbml,
        EventsForSbml,
        ParametersForSbml,
        ReactionsForSbml,
        SpeciesForSbml,
    )
    from mobspy.units.model_context import ModelUnitContext


def check(value: Any, message: str = "error") -> None:
    """If 'value' is None, prints an error message constructed using
    'message' and then exits with status code 1.  If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.
    """
    if value is None:
        raise SBMLError(f"LibSBML returned a null value trying to {message}.")

    if isinstance(value, int):
        if value == sbml.LIBSBML_OPERATION_SUCCESS:
            return
        err_msg = (
            "Error encountered trying to "
            + message
            + "."
            + "LibSBML returned error code "
            + str(value)
            + ': "'
            + sbml.OperationReturnValue_toString(value).strip()
            + '"'
        )
        raise SBMLError(err_msg)
    return


def create_model(  # noqa: PLR0913  # complex function signature
    species: SpeciesForSbml | None = None,
    parameters: ParametersForSbml | None = None,
    reactions: ReactionsForSbml | None = None,
    events: EventsForSbml | None = None,
    assignments: AssignmentsForSbml | None = None,
    model_context: ModelUnitContext | None = None,
) -> Any:
    """
    Returns an SBML Level 3 model.

    Example:

        species = { \
            'E': 1,
            'EM': 0,
            'EM2': 0,
            'F': 100, }

        parameters = { \
            'k': (1e-06, 'per_min'), }

        reactions = { \
            'Production_E': { \
                're': [(1, 'E'), (1, 'F')],
                'pr': [(2, 'E')],
                'kin': 'k * E * F', } }

        events = { \
            'e': { \
                'trigger': 'true',
                'delay': '10',
                'assignments': [('M', '1')], } }

    """

    if species is None:
        species = {}
    if parameters is None:
        parameters = {}
    if reactions is None:
        reactions = {}
    if events is None:
        events = {}
    if assignments is None:
        assignments = {}

    try:
        document = sbml.SBMLDocument(3, 1)
    except ValueError as exc:
        raise SBMLError("Could not create SBMLDocumentation object") from exc

    model = document.createModel()
    check(model, "create model")

    time_id, substance_id, spatial_dim = _setup_units(model, model_context)

    check(model.setTimeUnits(time_id), "set model-wide time units")
    check(model.setExtentUnits(substance_id), "set model units of extent")
    check(model.setSubstanceUnits(substance_id), "set model substance units")

    # Extract volume from parameters for compartment (remove from SBML params)
    vol_size = 1.0
    vol_units = "dimensionless"
    if "volume" in parameters:
        vol_size = float(parameters["volume"][0])
        vol_units = parameters["volume"][1]
        parameters = {k: v for k, v in parameters.items() if k != "volume"}

    _create_compartment(model, spatial_dim, vol_size, vol_units)
    _create_species(model, species, substance_id)
    _create_parameters(model, parameters)
    _create_reactions(model, reactions)
    _create_events(model, events)
    _create_assignments(model, assignments)

    return document


def _setup_units(
    model: Any,
    model_context: ModelUnitContext | None,
) -> tuple[str, str, int]:
    """Configure unit definitions and return (time_id, substance_id, spatial_dim)."""
    if model_context is not None:
        model_context.create_sbml_unit_definitions(model)
        return (
            model_context.get_sbml_time_units_id(),
            model_context.get_sbml_substance_units_id(),
            model_context.dimension,
        )

    per_second = model.createUnitDefinition()
    check(per_second, "create unit definition")
    check(per_second.setId("per_second"), "set unit definition id")
    unit = per_second.createUnit()
    check(unit, "create unit")
    check(unit.setKind(sbml.UNIT_KIND_SECOND), "set unit kind")
    check(unit.setExponent(-1), "set unit exponent")
    check(unit.setScale(0), "set unit scale")
    check(unit.setMultiplier(1), "set unit multiplier")
    return "second", "item", 3


def _create_compartment(
    model: Any,
    spatial_dim: int,
    size: float = 1.0,
    units: str = "dimensionless",
) -> None:
    """Create the compartment with the model's real volume."""
    c1 = model.createCompartment()
    check(c1, "create compartment")
    check(c1.setId("c1"), "set compartment id")
    check(c1.setConstant(True), 'set compartment "constant"')
    check(c1.setSize(size), 'set compartment "size"')
    check(c1.setSpatialDimensions(spatial_dim), "set compartment dimensions")
    check(c1.setUnits(units), "set compartment size units")


def _create_species(
    model: Any,
    species: SpeciesForSbml,
    substance_id: str,
) -> None:
    """Add species to the SBML model."""
    for s_str, s_val in species.items():
        s = model.createSpecies()
        check(s, "create species")
        check(s.setId(s_str), "set species id")
        check(s.setCompartment("c1"), "set species compartment")
        check(s.setConstant(False), 'set "constant" attribute')
        check(s.setInitialAmount(float(s_val)), "set initial amount")
        check(s.setSubstanceUnits(substance_id), "set substance units")
        check(s.setBoundaryCondition(False), 'set "boundaryCondition"')
        check(s.setHasOnlySubstanceUnits(False), 'set "hasOnlySubstanceUnits"')


def _create_parameters(model: Any, parameters: ParametersForSbml) -> None:
    """Add parameters to the SBML model."""
    for k_str in parameters:
        k = model.createParameter()
        check(k, "create parameter k")
        check(k.setId(k_str), "set parameter id")
        check(k.setConstant(True), 'set parameter "constant"')
        check(k.setValue(parameters[k_str][0]), "set parameter value")
        check(k.setUnits(parameters[k_str][1]), "set parameter units")


def _create_reactions(model: Any, reactions: ReactionsForSbml) -> None:
    """Add reactions to the SBML model."""
    for r_str in reactions:
        r = model.createReaction()
        check(r, "create reaction")
        check(r.setId(r_str), "set reaction id")
        check(r.setReversible(False), "set reaction reversibility flag")
        check(r.setFast(False), 'set reaction "fast" attribute')

        for re_val, re_str in reactions[r_str].reactants:
            species_ref = r.createReactant()
            check(species_ref, "create reactant")
            check(species_ref.setSpecies(re_str), "assign reactant species")
            check(species_ref.setStoichiometry(re_val), "set stoichiometry")
            check(species_ref.setConstant(True), 'set "constant" on species')

        for pr_val, pr_str in reactions[r_str].products:
            species_ref = r.createProduct()
            check(species_ref, "create product")
            check(species_ref.setSpecies(pr_str), "assign product species")
            check(species_ref.setStoichiometry(pr_val), "set stoichiometry")
            check(species_ref.setConstant(True), 'set "constant" on species')

        math_ast = sbml.parseL3Formula(reactions[r_str].kinetics)
        kinetic_law = r.createKineticLaw()
        check(math_ast, "create AST for rate expression")
        check(kinetic_law, "create kinetic law")
        check(kinetic_law.setMath(math_ast), "set math on kinetic law")


def _create_events(model: Any, events: EventsForSbml) -> None:
    """Add events to the SBML model."""
    for e_str in events:
        e = model.createEvent()
        check(e, "create event")
        check(e.setId(e_str), "set id")
        check(e.setUseValuesFromTriggerTime(False), "set use values from trigger time")

        t = e.createTrigger()
        check(t, "create trigger")
        check(
            t.setMath(sbml.parseL3Formula(events[e_str].trigger)),
            "set trigger condition",
        )
        check(t.setPersistent(False), "default not persistent")
        check(t.setInitialValue(False), "default not initially true")
        check(
            e.getTrigger().getMath(),
            "Problem when creating the trigger condition. The trigger will not work.",
        )
        d = e.createDelay()
        check(d, "create delay")
        check(d.setMath(sbml.parseL3Formula(str(events[e_str].delay))), "set math")

        for ass in events[e_str].assignments:
            ea = e.createEventAssignment()
            check(ea, "create event assignment")
            check(ea.setVariable(str(ass[0])), "set variable")
            check(ea.setMath(sbml.parseL3Formula(str(ass[1]))), "set math")


def _create_assignments(model: Any, assignments: AssignmentsForSbml) -> None:
    """Add assignment rules to the SBML model."""
    for asg in assignments.values():
        assignment_rule = model.createAssignmentRule()
        check(assignment_rule, "create assignment rule")
        check(assignment_rule.setVariable(asg.species), "set assignment variable")
        math_ast = sbml.parseL3Formula(asg.expression)
        check(math_ast, f"create AST for assignment expression: {asg.expression}")
        check(assignment_rule.setMath(math_ast), "set math on assignment rule")
