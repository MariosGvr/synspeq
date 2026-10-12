"""Planner stocks: AiZynthFinder loads a stock as a set of InChIKeys. A catalogue the user imported as SMILES or SDF must be
converted before planning; passing the text file itself would load as a set of SMILES strings and match nothing."""
from __future__ import annotations
import json
import os
from typing import Optional
from ..hashing import content_id
from .store import validate
from .union import build_union


def is_compiled_stock(path: str) -> bool:
    """HDF5 with an `inchi_key` table, or a text file of InChIKeys."""
    low = path.lower()
    if not (low.endswith((".h5", ".hdf5")) or low.endswith(".txt")):
        return False
    return validate(path, "aizynthfinder_hdf5").get("status") in ("ok", "unchecked") and (
        low.endswith((".h5", ".hdf5")) or validate(path, "aizynthfinder_hdf5").get("status") == "ok")


def compiled_planner_stock(rid: str, path: str, cache_dir: str, hash_cache: Optional[str] = None, flags: Optional[dict] = None) -> str:
    """Path of an InChIKey stock for resource `rid`. Already compiled stocks are used as they are; other catalogues are
    converted once and the result is cached under cache_dir, keyed by the content hash of the source."""
    if os.path.isdir(path):
        for n in ("zinc_stock.hdf5", "stock.hdf5"):
            if os.path.exists(os.path.join(path, n)):
                return os.path.join(path, n)
    elif is_compiled_stock(path):
        return path
    cid = content_id(path, hash_cache)
    out_dir = os.path.join(cache_dir, "stocks")
    name = f"{rid}_{cid[:16]}"
    man = os.path.join(out_dir, f"{name}_manifest.json")
    if os.path.exists(man):
        try:
            m = json.load(open(man))
            if os.path.exists(m["file"]) and m.get("source_sha256") == cid:
                return m["file"]
        except (OSError, ValueError, KeyError):
            pass
    info = build_union({rid: path}, out_dir, flags={rid: flags or {}}, name=name)
    info["source_sha256"] = cid
    if info["n_unique_inchikeys"] == 0:
        raise ValueError(f"{rid}: {path} yielded no usable molecules")
    with open(man, "w") as fh:
        json.dump(info, fh, indent=1)
    return info["file"]
