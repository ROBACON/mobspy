# Changelog

All notable changes to MobsPy will be documented in this file.

This project follows [Semantic Versioning](https://semver.org/). The public API
consists of everything exported from `mobspy.__init__` (listed in `__all__`).

## [Unreleased]

### Changed
- Package versions are derived from Git tags; releases no longer require editing a version field.

## [3.0.0rc1] — release candidate

### Added
- Reusable `Model` snapshots with species-owned declarations and weak session indexes.
- Typed `ExecutionPlan`, `RunSettings`, and backend result contracts.
- CI tests the built wheel on Python 3.11–3.14 and builds the documentation.
- CI executes all tutorial notebooks and example scripts, including plots and the full XOR sweep.
- Separate CI coverage for the BioCRNpyler, PySB, and BioNetGen reference models.
- Tag and artifact version checks, reusable CI, and Trusted Publishing release automation.

### Fixed
- Repeated compilation preserves user units and numerical results.
- Count and parameter updates modify the authoritative compiled model only.
- Species updates accept dot notation and internal species names.
- Composed simulations own their execution plans; repeated runs and exports preserve child simulations.
- Composed stages convert time, volume, and substance units at their boundaries.
- Unit-bearing event times survive compilation and recompilation.
- COPASI model lifecycle and repetition state are isolated during execution.
- Unit-bearing parameter updates and sweeps resolve independently in each stage.
- Antimony exports include every stage; composed SBML preserves transition times and amounts.
- JSON saving includes serializable model metadata and respects output filenames.
- Inheritance, assignment, and plotting tutorials now demonstrate and check their stated behavior.
- The phage example includes its plot configuration; the BioNetGen reference uses the required `.bngl` extension.
- Plot configuration rejects invalid positional arguments and copies nested species styles.

### Migration
- This is a major release. See `docs/migration_guide.md` for renamed module paths.
- User configuration retains its original values; normalized settings live in execution plans.
- `Simulation.delete()` releases its own state. Shared species and Model objects are unchanged.
- Compositions own their plot configuration. Compiled initial counts are numeric.
- Custom backends implement `generate_model(ConcreteModel)` and `run(ExecutionPlan, jobs)`.


### Changed
- Codebase rewrite: strict linting, type safety, thread safety, and code cleanup.
- `@` is now the preferred rate syntax; `[]` emits `DeprecationWarning`.
- `Rev[]` emits `DeprecationWarning`; use tuple rates `A >> B @ (k_fwd, k_rev)`.
- `Simulation` accepts a pluggable `backend` parameter (default: `SBMLBackend`).
- `PlotConfigProxy` replaced by `PlotConfig` (dict subclass with attribute access).
- `Simulation_Utils` and `Experimental_Data_Holder` inlined into `Simulation`.
- `OverrideUnitRegistry` and `u` extracted to `mobspy.units.registry`.
- Expression engine operator boilerplate collapsed via dispatch tables.
- All tutorials, examples, and documentation updated to use `@` syntax.

## [2.8.0]

### Added
- Initial ODE syntax (experimental, not yet compatible with all features).

## [2.5.0]

### Added
- `generate_antimony()` function to translate MobsPy models to Antimony format.

## [2.4.4]

### Added
- Full assignment rule support (both notations).

## [2.3.0]

### Changed
- Standard output is now always concentration.

## [2.2.0]

### Added
- MobsPy expressions for complex rate functions.

## [2.1.0]

### Added
- Model parameters support.
- `MySim.fres` attribute for accessing only the first time series.

### Changed
- `MySim.results` now returns a list of all resulting time series regardless of repetition count.

## [2.0.1]

### Added
- Events.
- Simulation concatenation via `+` operator.

### Changed
- Simplified results access: `MySim.results["MetaSpeciesName"]` replaces
  `MySim.results["data"]["runs"]["MetaSpeciesName"]`.

## [1.1.0]

### Added
- Automated testing on git pushes.
