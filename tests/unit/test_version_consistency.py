"""Package version consistency: zt005.__version__ vs pyproject.toml.

Guards against the recurring drift that left src/zt005/__init__.py at
"0.2.0" while pyproject was at 12.x for dozens of releases.
"""
import tomllib
from pathlib import Path

import zt005

PYPROJECT = Path(__file__).resolve().parents[2] / "pyproject.toml"


def test_package_version_matches_pyproject():
    with PYPROJECT.open("rb") as f:
        pyproject_version = tomllib.load(f)["project"]["version"]
    assert zt005.__version__ == pyproject_version, (
        f"src/zt005/__init__.py says {zt005.__version__!r}, "
        f"pyproject says {pyproject_version!r}"
    )
