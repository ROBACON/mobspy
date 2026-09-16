"""Write simulation results using the same rules for single and chained runs."""

from __future__ import annotations

import json
from pathlib import Path

from mobspy.exceptions import SimulationError
from mobspy.results.time_series import SimulationResults


def save_results(results: object, file: str | None, default_file: str | None) -> None:
    """Save result JSON, creating the requested output directory if needed."""
    if not isinstance(results, SimulationResults):
        raise SimulationError("No simulation results available to save")
    filename = file if file is not None else default_file
    if not filename:
        raise SimulationError("No default output file specified in parameters")
    if not filename.endswith(".json"):
        filename += ".json"
    path = Path(filename)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as stream:
            json.dump(results.to_dict(), stream, indent=4)
    except OSError as error:
        raise SimulationError(f"Error saving data to file: {error}") from error
    except (TypeError, ValueError) as error:
        raise SimulationError(f"Error serializing simulation data: {error}") from error
