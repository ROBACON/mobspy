"""Run doctests from mobspy modules in isolated subprocesses.

MobsPy uses global state (species registries, reaction sets) that gets
polluted when doctests run in the same process as regular tests. Each
doctest module is tested via subprocess to avoid cross-contamination.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

DOCTEST_MODULES = [
    "mobspy/constants.py",
    "mobspy/exceptions.py",
    "mobspy/simulation_config.py",
    "mobspy/types.py",
    "mobspy/modules/species.py",
    "mobspy/modules/species_constructors.py",
    "mobspy/simulation.py",
]


@pytest.mark.parametrize("module", DOCTEST_MODULES)
def test_module_doctests(module: str) -> None:
    """Run doctests for a single module in a clean subprocess."""
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "--doctest-modules",
            module,
            "--no-header",
            "-q",
            "--no-cov",
            "--override-ini=addopts=",
        ],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert result.returncode == 0, (
        f"Doctests failed for {module}:\n{result.stdout}\n{result.stderr}"
    )
