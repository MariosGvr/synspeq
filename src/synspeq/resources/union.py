"""Union planner stock from several catalogues.

Every member is reduced to InChIKeys (AiZynthFinder's stock key): HDF5 stocks are read as they are, SMILES/SDF/CSV
catalogues are parsed with RDKit. Records RDKit cannot parse are counted, never dropped silently. The result is one
HDF5 (or, when PyTables is missing, a plain text file with one key per line, which AiZynthFinder also reads) plus a
manifest with member counts, overlaps and licence flags. The union is for the user's own runs: it is marked
redistributable only if every member is.
"""
from __future__ import annotations
import hashlib
import json
import os
from typing import Dict, Iterator, Optional, Set, Tuple
from rdkit import Chem


def _keys_from_hdf5(path: str) -> Iterator[str]:
    import pandas as pd
    df = pd.read_hdf(path, key="table")
    yield from (str(k) for k in df["inchi_key"].values)


def member_keys(path: str, smiles_col: str = "smiles") -> Tuple[Set[str], dict]:
    from ..purchasing.catalogue import _iter_smiles, catalogue_files, is_catalogue_file
    keys: Set[str] = set()
    n_in = n_bad = 0
    if os.path.isdir(path):
        files = sorted(os.path.join(path, f) for f in os.listdir(path)
                       if os.path.isfile(os.path.join(path, f)) and not f.startswith(".") and (is_catalogue_file(f) or f.lower().endswith((".h5", ".hdf5"))))
    else:
        files = [path]
    for f in files:
        if f.lower().endswith((".h5", ".hdf5")):
            for k in _keys_from_hdf5(f):
                n_in += 1
                keys.add(k)
            continue
        for s in _iter_smiles(f, smiles_col):
            n_in += 1
            m = Chem.MolFromSmiles(s)
            if m is None:
                n_bad += 1
                continue
            keys.add(Chem.MolToInchiKey(m))
    return keys, dict(n_records=n_in, n_unparseable=n_bad, n_unique_keys=len(keys))


def build_union(members: Dict[str, str], out_dir: str, flags: Optional[Dict[str, dict]] = None, name: str = "union_stock") -> dict:
    os.makedirs(out_dir, exist_ok=True)
    allk: Set[str] = set()
    per: Dict[str, dict] = {}
    seen_before: Set[str] = set()
    for rid, path in members.items():
        k, info = member_keys(path)
        info["new_keys_added"] = len(k - allk)
        info["overlap_with_previous"] = len(k & allk)
        per[rid] = info
        allk |= k
    try:
        import tables  # noqa: F401
        import pandas as pd
        out = os.path.join(out_dir, f"{name}.hdf5")
        tmp = f"{out}.{os.getpid()}.tmp"
        pd.DataFrame({"inchi_key": sorted(allk)}).to_hdf(tmp, key="table", mode="w")
        os.replace(tmp, out)
        fmt = "hdf5"
    except ImportError:
        out = os.path.join(out_dir, f"{name}.txt")
        tmp = f"{out}.{os.getpid()}.tmp"
        with open(tmp, "w") as fh:
            fh.write("\n".join(sorted(allk)) + "\n")
        os.replace(tmp, out)
        fmt = "text"
    h = hashlib.sha256()
    with open(out, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    fl = {rid: (flags or {}).get(rid, {}) for rid in members}
    man = dict(file=out, format=fmt, n_unique_inchikeys=len(allk), sha256=h.hexdigest(), members=per, member_flags=fl,
               noncommercial=any(f.get("noncommercial") for f in fl.values()),
               sharealike=any(f.get("sharealike") for f in fl.values()),
               redistributable=all(f.get("redistributable", False) for f in fl.values()) if fl else False,
               note="derived from the members above; follow each member's terms before sharing this file")
    with open(os.path.join(out_dir, f"{name}_manifest.json"), "w") as fh:
        json.dump(man, fh, indent=1)
    return man
