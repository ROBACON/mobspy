"""Release validation rejects mismatched tags and distribution metadata."""

import io
import runpy
import tarfile
import zipfile
from pathlib import Path

import pytest

release = runpy.run_path(str(Path(__file__).parents[1] / "scripts/check_release.py"))


@pytest.mark.parametrize("version", ["3.0.0rc1", "3.0.0", "3.1.2"])
def test_release_accepts_matching_canonical_tags(version):
    assert str(release["validate_version"]("v" + version, version)) == version


@pytest.mark.parametrize(
    ("tag", "version"),
    [
        ("v3.0.0", "3.0.0rc1"),
        ("3.0.0", "3.0.0"),
        ("v3.0.0+local", "3.0.0+local"),
        ("v3.0.0rc1", "3.0.0-rc1"),
    ],
)
def test_release_rejects_mismatched_or_noncanonical_tags(tag, version):
    with pytest.raises(ValueError):
        release["validate_version"](tag, version)


def distribution_pair(directory, wheel_version="3.0.0rc1", source_version="3.0.0rc1"):
    with zipfile.ZipFile(directory / "mobspy.whl", "w") as wheel:
        wheel.writestr(
            "mobspy.dist-info/METADATA",
            f"Name: mobspy\nVersion: {wheel_version}\n",
        )
    with tarfile.open(directory / "mobspy.tar.gz", "w:gz") as source:
        metadata = f"Name: mobspy\nVersion: {source_version}\n".encode()
        member = tarfile.TarInfo("mobspy/PKG-INFO")
        member.size = len(metadata)
        source.addfile(member, io.BytesIO(metadata))


def test_release_checks_metadata_inside_both_archives(tmp_path):
    distribution_pair(tmp_path)
    release["validate_artifacts"](tmp_path, "3.0.0rc1")
    distribution_pair(tmp_path, source_version="2.8.0")
    with pytest.raises(ValueError, match="metadata"):
        release["validate_artifacts"](tmp_path, "3.0.0rc1")
    distribution_pair(tmp_path, wheel_version="2.8.0")
    with pytest.raises(ValueError, match="metadata"):
        release["validate_artifacts"](tmp_path, "3.0.0rc1")


def test_release_rejects_missing_or_extra_archives(tmp_path):
    with pytest.raises(ValueError, match="one portable wheel"):
        release["validate_artifacts"](tmp_path, "3.0.0rc1")
    distribution_pair(tmp_path)
    (tmp_path / "old.whl").touch()
    with pytest.raises(ValueError, match="one portable wheel"):
        release["validate_artifacts"](tmp_path, "3.0.0rc1")
