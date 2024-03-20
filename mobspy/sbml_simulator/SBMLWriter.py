"""
    This module is responsible for converting a model_str into a SBML format
"""
import libsbml as sbml


def check(value, message="error"):
    """If 'value' is None, prints an error message constructed using
    'message' and then exits with status code 1.  If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.
    """
    if value is None:
        raise RuntimeError(f"LibSBML returned a null value trying to {message}.")

    elif type(value) is int:
        if value == sbml.LIBSBML_OPERATION_SUCCESS:
            return
        else:
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
            raise RuntimeError(err_msg)
    else:
        return


def create_model(species={}, parameters={}, reactions={}, events={}, assignments={}):
    """
    Returns an SBML Level 3 model.
    Example:
    species = { 'E': 1,
                'EM': 0,
                'EM2': 0,
                'F': 100,
                },
    parameters = {'k': (1e-06,'per_min'),
                 }
    reactions = { 'Production_E':
                        { 're': [(1,'E'),(1,'F')],
                          'pr': [(2,'E')],
                          'kin' : 'k * E * F'
                        },
                },
    events = {'e':
              { 'trigger': 'true',
                'delay': '10',
                'assignments': [('M','1'),],
              },
    }
    """

    # Create an empty SBMLDocument object.  It's a good idea to check for
    # possible errors.  Even when the parameter values are hardwired like
    # this, it is still possible for a failure to occur (e.g., if the
    # operating system runs out of memory).

    try:
        document = sbml.SBMLDocument(3, 1)
    except ValueError:
        raise RuntimeError("Could not create SBMLDocumention object")

    # Create the basic Model object inside the SBMLDocument object.

    model = document.createModel()
    check(model, "create model")
    check(model.setTimeUnits("second"), "set model-wide time units")
    check(model.setExtentUnits("item"), "set model units of extent")
    check(
        model.setSubstanceUnits("item"), "set model substance units"
    )  # mole, item, gram, kilogram, dimensionless

    # Create a unit definition we will need later.

    per_second = model.createUnitDefinition()
    check(per_second, "create unit definition")
    check(per_second.setId("per_min"), "set unit definition id")
    unit = per_second.createUnit()
    check(unit, "create unit")
    check(unit.setKind(sbml.UNIT_KIND_SECOND), "set unit kind")
    check(unit.setExponent(-1), "set unit exponent")
    check(unit.setScale(0), "set unit scale")
    check(
        unit.setMultiplier(1), "set unit multiplier"
    )

    # Create a compartment inside this model

    c1 = model.createCompartment()
    check(c1, "create compartment")
    check(c1.setId("c1"), "set compartment id")
    check(c1.setConstant(True), 'set compartment "constant"')
    check(c1.setSize(1), 'set compartment "size"')
    check(c1.setSpatialDimensions(3), "set compartment dimensions")
    check(
        c1.setUnits("dimensionless"), "set compartment size units"
    )

    # Create species inside this model, set the required attributes
    # for each species in SBML Level 3 (which are the 'id', 'compartment',
    # 'constant', 'hasOnlySubstanceUnits', and 'boundaryCondition'
    # attributes), and initialize the amount of the species along with the
    # units of the amount.

    for s_str, s_val in species.items():
        s = model.createSpecies()
        check(s, "create species")
        check(s.setId(s_str), "set species id")
        check(s.setCompartment("c1"), "set species compartment")
        check(s.setConstant(False), 'set "constant" attribute')
        check(s.setInitialAmount(float(s_val)), "set initial amount")
        check(s.setSubstanceUnits("item"), "set substance units")
        check(s.setBoundaryCondition(False), 'set "boundaryCondition"')
        check(s.setHasOnlySubstanceUnits(False), 'set "hasOnlySubstanceUnits"')

    # Create a parameter object inside this model, set the required
    # attributes 'id' and 'constant' for a parameter in SBML Level 3, and
    # initialize the parameter with a value along with its units.

    for k_str in parameters:
        k = model.createParameter()
        check(k, "create parameter k")
        check(k.setId(k_str), "set parameter id")
        check(k.setConstant(True), 'set parameter "constant"')
        check(k.setValue(parameters[k_str][0]), "set parameter value")
        check(k.setUnits(parameters[k_str][1]), "set parameter units")

    # Create a reaction inside this model, set the reactants and products,
    # and set the reaction rate expression (the SBML "kinetic law").  We
    # set the minimum required attributes for all of these objects.  The
    # units of the reaction rate are determined from the 'timeUnits' and
    # 'extentUnits' attributes on the Model object.

    for r_str in reactions:
        r = model.createReaction()
        check(r, "create reaction")
        check(r.setId(r_str), "set reaction id")
        check(r.setReversible(False), "set reaction reversibility flag")
        check(r.setFast(False), 'set reaction "fast" attribute')

        reactants = reactions[r_str]["re"]
        for re_val, re_str in reactants:
            species_ref = r.createReactant()
            check(species_ref, "create reactant")
            check(species_ref.setSpecies(re_str), "assign reactant species")
            check(species_ref.setStoichiometry(re_val), "set set stoichiometry")
            check(species_ref.setConstant(True), 'set "constant" on species')

        products = reactions[r_str]["pr"]
        for pr_val, pr_str in products:
            species_ref = r.createProduct()
            check(species_ref, "create product")
            check(species_ref.setSpecies(pr_str), "assign product species")
            check(species_ref.setStoichiometry(pr_val), "set set stoichiometry")
            check(species_ref.setConstant(True), 'set "constant" on species')

        math_ast = sbml.parseL3Formula(reactions[r_str]["kin"])
        kinetic_law = r.createKineticLaw()
        check(math_ast, f"create AST for rate expression")
        check(kinetic_law, "create kinetic law")
        check(kinetic_law.setMath(math_ast), "set math on kinetic law")

    # create events
    for e_str in events:
        e = model.createEvent()
        check(e, "create event")
        check(e.setId(e_str), "set id")
        check(e.setUseValuesFromTriggerTime(False), "?")

        t = model.createTrigger()
        check(t, "create trigger")
        check(
            t.setMath(sbml.parseL3Formula(events[e_str]["trigger"])),
            "set trigger condition",
        )
        check(t.setPersistent(False), "default not persistent")
        check(t.setInitialValue(False), "default not initially true")
        check(e.getTrigger().getMath(), 'Problem when creating the trigger condition. The trigger will not work.')
        # print( '> ' + sbml.formulaToString(e.getTrigger().getMath()) )

        d = model.createDelay()
        check(d, "create delay")
        check(d.setMath(sbml.parseFormula(events[e_str]["delay"])), "set math")
        check(e.setDelay(d), "set delay")

        for ass in events[e_str]["assignments"]:
            ea = model.createEventAssignment()
            check(ea, "check event assignment")
            check(ea.setVariable(ass[0]), "set variable")
            check(ea.setMath(sbml.parseL3Formula(ass[1])), "set math")

    for _, asg in assignments.items():

        species_id = asg['species']
        expression = asg['expression']

        # Create an AssignmentRule for species
        assignment_rule = model.createAssignmentRule()
        check(assignment_rule, "create assignment rule")
        check(assignment_rule.setVariable(species_id), "set assignment variable")
        math_ast = sbml.parseL3Formula(expression)
        check(math_ast, f"create AST for assignment expression: {expression}")
        check(assignment_rule.setMath(math_ast), "set math on assignment rule")

    return document  # sbml.writeSBMLToString(document)
