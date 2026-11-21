"""Versions of the packages a result depends on.

Molecule identity is the RDKit canonical isomeric SMILES, so another RDKit release can change identities and with them
molecule ids and cached planner records. Reports and run manifests therefore record the versions, and a run prints a
note when RDKit differs from the release the results of the paper were made with."""
from __future__ import annotations
import platform
from importlib import metadata
from typing import Dict, Optional

VALIDATED_RDKIT = "2023.09"          # release series behind the reported results

PACKAGES = ("rdkit", "numpy", "pandas", "aizynthfinder", "rxn-insight", "rxnmapper")


def _version(pkg: str) -> Optional[str]:
    try:
        return metadata.version(pkg)
    except metadata.PackageNotFoundError:
        return None


def environment_info() -> Dict[str, Optional[str]]:
    from . import __version__
    info = {"synspeq": __version__, "python": platform.python_version()}
    for pkg in PACKAGES:
        info[pkg] = _version(pkg)
    if info["rdkit"] is None:
        import rdkit
        info["rdkit"] = rdkit.__version__
    return info


def rdkit_note() -> Optional[str]:
    """A one-line note when RDKit is not the validated release series, else None."""
    import rdkit
    v = rdkit.__version__
    parts = v.replace("_", ".").split(".")
    series = ".".join(parts[:2])
    if series == VALIDATED_RDKIT:
        return None
    return (f"RDKit {v} differs from the release series {VALIDATED_RDKIT} used for the reported results; canonical SMILES, and with "
            "them molecule ids and cached planner records, can differ. Use environment.yml or constraints.txt for comparable runs.")
