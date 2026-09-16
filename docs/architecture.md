# Architecture

MobsPy separates model declarations, compilation, execution, and presentation.

1. **DSL and Model:** Species own their reactions and initial counts. The session
   registry is a weak discovery index. `Model` captures only the selected species'
   declarations and their parameter namespace.
2. **Compiler:** `compile_model` expands declarations into a `ConcreteModel` with
   species, reactions, events, parameters, and a resolved unit context.
3. **Simulation:** User configuration retains its units. Compilation resolves a
   private configuration copy and stores one authoritative compiled model.
   Compatibility attributes derive from that model.
4. **Execution plan:** `build_execution_plan` creates independent parameter sweep
   variants and numeric `RunSettings`. A composition owns its plan and results;
   it does not use its first child as an execution-state container.
5. **Backend:** `generate_model(model)` exports a compiled model. `run(plan, jobs)`
   returns `BackendResults`: sweep → repetition → time-series dictionary. Each
   dictionary contains `Time` and species trajectories in native concentration
   units. Across a chain, time and concentration use the first stage's units.
6. **Results:** Result processing converts backend output into requested units,
   then wraps it in `SimulationResults`. Plotting consumes those results.

## Ownership rules

- A simulation's compilation and parameter updates cannot change another
  simulation's configuration or parameter objects.
- Backend execution must not mutate the plan's models or settings. Repetitions
  start from the same initial conditions.
- DSL declarations and event expressions retain their original units.
- `delete()` releases the receiver's state; shared models and species remain
  usable.
- Species characteristics should be completed before taking a Model snapshot.
  Do not mutate species or share a Simulation instance while compiling it in
  another thread.

The default backend coordinates the complete BasiCO/COPASI data-model lifecycle
within a process because the solver maintains a global model registry. Separate
processes can execute independent simulations concurrently.

## Export boundaries

Sequential runs convert time, volume, and substance units between stages. Exporting
an entire chain as one SBML or Antimony model with `compose=True` requires identical
unit systems, volumes, and assignment rules. Unsupported combinations raise an
error; `generate_sbml()` and `generate_antimony()` export each stage separately.

Result dictionaries and JSON exports include serializable model metadata and
parameter-sweep values. Saving a composition writes its own results and respects
`output_dir` and `output_file`.
