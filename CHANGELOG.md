# Changelog

All notable changes to MobsPy will be documented in this file.

This project follows [Semantic Versioning](https://semver.org/). The public API
consists of everything exported from `mobspy.__init__` (listed in `__all__`).

## [Unreleased]

### Changed
- Codebase rewrite: strict linting, type safety, thread safety, and code cleanup.

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
