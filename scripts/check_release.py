"""Validate release tags and built distribution versions before publication."""

from __future__ import annotations

import argparse
import email
import os
import tarfile
import zipfile
from pathlib import Path

from packaging.version import Version
from setuptools_scm import get_version


def validate_version(tag: str, version: str) -> Version:
    """Require an explicit, canonical v-prefixed tag matching package metadata."""
    parsed = Version(version)
    if parsed.local is not None or tag != f"v{parsed}" or version != str(parsed):
        raise ValueError(f"Tag {tag!r} must match the package version v{parsed}")
    return parsed


def validate_artifacts(directory: Path, version: str) -> None:
    """Check the actual metadata of the wheel and source distribution."""
    wheels = list(directory.glob("*.whl"))
    sources = list(directory.glob("*.tar.gz"))
    if len(wheels) != 1 or len(sources) != 1:
        raise ValueError("Expected one portable wheel and one source distribution")
    with zipfile.ZipFile(wheels[0]) as wheel:
        metadata_name = next(
            n for n in wheel.namelist() if n.endswith(".dist-info/METADATA")
        )
        wheel_metadata = wheel.read(metadata_name).decode()
    with tarfile.open(sources[0]) as source:
        member = next(
            m
            for m in source.getmembers()
            if m.name.count("/") == 1 and m.name.endswith("/PKG-INFO")
        )
        metadata_file = source.extractfile(member)
        if metadata_file is None:
            raise ValueError("Source distribution has no package metadata")
        source_metadata = metadata_file.read().decode()
    for metadata in (wheel_metadata, source_metadata):
        package = email.message_from_string(metadata)
        if package["Name"] != "mobspy" or package["Version"] != version:
            raise ValueError("Distribution metadata does not match the release")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag")
    parser.add_argument("--artifacts", type=Path)
    args = parser.parse_args()
    version = Version(get_version(root=Path(__file__).resolve().parents[1]))
    if args.tag is not None:
        validate_version(args.tag, str(version))
    if args.artifacts is not None:
        validate_artifacts(args.artifacts, str(version))
    if output := os.environ.get("GITHUB_OUTPUT"):
        with Path(output).open("a", encoding="utf-8") as file:
            prerelease = version.is_prerelease or version.is_devrelease
            file.write(f"prerelease={str(prerelease).lower()}\n")


if __name__ == "__main__":
    main()
