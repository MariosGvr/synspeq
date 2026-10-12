"""Local resource store: a registry of where the user's own copies are, with size, hash and validation.

Nothing is copied, uploaded or shared. A registered resource is a reference (absolute path) plus the facts
needed for the run manifest: size, sha256 (files; directories get a file count and a listing hash),
record count, validation result, snapshot date (file modification time unless given), import time.
"""
from __future__ import annotations
import datetime as dt
import gzip
import hashlib
import json
import os
from typing import Dict, Optional

ENV = "SYNSPEQ_DATA_DIR"


def default_dir() -> str:
    return os.environ.get(ENV) or os.path.join(os.getcwd(), ".synspeq")   # project-local, like .git


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _open_text(path: str):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path, "rt", errors="replace")


def _count_sdf(path: str, sample: int = 2000):
    from rdkit import Chem
    n = ok = bad = 0
    with _open_text(path) as fh:
        buf = []
        for line in fh:
            if line.startswith("$$$$"):
                n += 1
                if n <= sample:
                    blk = "".join(buf)
                    cut = blk.find("M  END")
                    m = Chem.MolFromMolBlock(blk[: cut + 6] + "\n", sanitize=True) if cut >= 0 else None
                    ok, bad = (ok + 1, bad) if m is not None else (ok, bad + 1)
                buf = []
            else:
                buf.append(line)
    return n, {"sampled": min(n, sample), "parse_failures_in_sample": bad}


def _count_smiles(path: str, sample: int = 2000):
    from rdkit import Chem
    n = bad = 0
    with _open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            n += 1
            if n <= sample and Chem.MolFromSmiles(line.split()[0].split("\t")[0]) is None:
                bad += 1
    return n, {"sampled": min(n, sample), "parse_failures_in_sample": bad}


def validate(path: str, fmt: Optional[str]) -> dict:
    """Light, format-aware validation. Returns {status: ok|invalid|unchecked, ...}."""
    if not os.path.exists(path):
        return {"status": "invalid", "reason": "path_not_found"}
    if os.path.isdir(path):
        files = sorted(f for f in os.listdir(path) if not f.startswith("."))
        if not files:
            return {"status": "invalid", "reason": "empty_directory"}
        return {"status": "ok", "n_files": len(files)}
    if os.path.getsize(path) == 0:
        return {"status": "invalid", "reason": "empty_file"}
    low = path.lower()
    try:
        if fmt == "aizynthfinder_hdf5":
            with open(path, "rb") as fh:
                head = fh.read(8)
            if low.endswith((".h5", ".hdf5")):
                return {"status": "ok" if head == b"\x89HDF\r\n\x1a\n" else "invalid", **({} if head == b"\x89HDF\r\n\x1a\n" else {"reason": "not_an_hdf5_file"}), "size_bytes": os.path.getsize(path)}
            with open(path) as fh:                                  # text stock: one InChIKey per line
                lines = [l.strip() for l in fh if l.strip()][:2000]
            good = sum(1 for l in lines if len(l) == 27 and l[14] == "-" and l[25] == "-")
            return {"status": "ok" if lines and good / len(lines) > 0.8 else "invalid", "n_sampled": len(lines)}
        if low.endswith((".sdf", ".sdf.gz")):
            n, info = _count_sdf(path)
        elif low.endswith((".smi", ".smi.gz", ".txt", ".txt.gz", ".csv", ".csv.gz")):
            n, info = _count_smiles(path)
        elif fmt == "parquet" or low.endswith((".parquet", ".gzip")):
            with open(path, "rb") as fh:
                head = fh.read(4)
                fh.seek(-4, 2)
                tail = fh.read(4)
            if head != b"PAR1" or tail != b"PAR1":
                return {"status": "invalid", "reason": "not_a_parquet_file"}
            return {"status": "ok", "size_bytes": os.path.getsize(path)}
        elif fmt == "sdf":
            n, info = _count_sdf(path)
        elif fmt in ("smiles", "sdf_or_smiles"):
            n, info = _count_smiles(path)
        else:
            return {"status": "unchecked", "reason": f"no validator for format {fmt!r}"}
    except Exception as e:                   # unreadable / wrong format
        return {"status": "invalid", "reason": f"{type(e).__name__}: {e}"}
    if n == 0:
        return {"status": "invalid", "reason": "no_records"}
    if info["sampled"] and info["parse_failures_in_sample"] / info["sampled"] > 0.2:
        return {"status": "invalid", "reason": "more_than_20pct_parse_failures_in_sample", "n_records": n, **info}
    return {"status": "ok", "n_records": n, **info}


class Store:
    def __init__(self, root: Optional[str] = None):
        self.root = root or default_dir()
        self.path = os.path.join(self.root, "registry.json")

    def _load(self) -> Dict[str, dict]:
        if not os.path.exists(self.path):
            return {}
        with open(self.path) as fh:
            return json.load(fh)

    def _save(self, reg: Dict[str, dict]) -> None:
        os.makedirs(self.root, exist_ok=True)
        tmp = self.path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(reg, fh, indent=1, sort_keys=True)
        os.replace(tmp, self.path)           # atomic

    def get(self, rid: str) -> Optional[dict]:
        return self._load().get(rid)

    def all(self) -> Dict[str, dict]:
        return self._load()

    def register(self, rid: str, path: str, fmt: Optional[str] = None, hash_it: bool = True, coverage: Optional[dict] = None,
                 snapshot_date: Optional[str] = None, terms_acknowledged: bool = False) -> dict:
        path = os.path.abspath(os.path.expanduser(path))
        v = validate(path, fmt)
        entry = dict(path=path, format=fmt, validation=v, coverage=coverage, terms_acknowledged=bool(terms_acknowledged),
                     imported_at=dt.datetime.now().isoformat(timespec="seconds"))
        if os.path.exists(path):
            entry["snapshot_date"] = snapshot_date or dt.datetime.fromtimestamp(os.path.getmtime(path)).date().isoformat()
            if os.path.isdir(path):
                names = sorted(os.listdir(path))
                entry["size_bytes"] = sum(os.path.getsize(os.path.join(path, f)) for f in names if os.path.isfile(os.path.join(path, f)))
                entry["n_files"] = len(names)
                entry["listing_sha256"] = hashlib.sha256("\n".join(f"{f}:{os.path.getsize(os.path.join(path, f))}" for f in names if os.path.isfile(os.path.join(path, f))).encode()).hexdigest()
            else:
                entry["size_bytes"] = os.path.getsize(path)
                if hash_it:
                    entry["sha256"] = _sha256(path)
        reg = self._load()
        reg[rid] = entry
        self._save(reg)
        return entry
