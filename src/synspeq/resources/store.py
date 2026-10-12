"""Local resource store: a registry of where the user's own copies are, with size, hash and validation.

Nothing is copied, uploaded or shared. A registered resource is a reference (absolute path) plus the facts
needed for the run manifest: size, modification time, sha256 (files; directories get a file count and a listing
hash), validation result, snapshot date, import time. A resource whose file changed since import is not trusted.

Validation reads the same way the loaders do (same SDF record splitting, same SMILES readers) and requires at least one
usable record; a magic number alone does not make a file valid.
"""
from __future__ import annotations
import datetime as dt
import json
import os
import uuid
from typing import Dict, Optional
from ..hashing import file_sha256, listing_id

ENV = "SYNSPEQ_DATA_DIR"


def default_dir() -> str:
    return os.environ.get(ENV) or os.path.join(os.getcwd(), ".synspeq")   # project-local, like .git


# ------------------------------------------------------------------ validation
def _sample_smiles(path: str, sample: int = 2000):
    """-> (n_records, n_sampled, n_parse_failures) using the loaders' own readers."""
    from rdkit import Chem
    from ..purchasing.catalogue import _iter_smiles
    low = path.lower()
    n = nbad = ns = 0
    for s in _iter_smiles(path):
        n += 1
        if ns < sample:
            ns += 1
            if Chem.MolFromSmiles(s) is None:
                nbad += 1
    if low.endswith((".sdf", ".sdf.gz")):
        pass                                   # SDF records that fail to parse are skipped by the reader; count them below
    return n, ns, nbad


def _count_sdf_records(path: str) -> int:
    from ..io.sdf import iter_sdf_blocks
    return sum(1 for _ in iter_sdf_blocks(path))


def _validate_file(path: str, fmt: Optional[str], check: dict) -> dict:
    low = path.lower()
    if os.path.getsize(path) == 0:
        return {"status": "invalid", "reason": "empty_file"}
    if fmt is None and low.endswith(".txt"):               # a text file of InChIKeys is a compiled stock, not a SMILES list
        with open(path, errors="replace") as fh:
            head = [l.strip() for l in fh if l.strip()][:20]
        if head and all(len(l) == 27 and l[14] == "-" and l[25] == "-" for l in head):
            fmt = "aizynthfinder_hdf5"
    if fmt == "aizynthfinder_hdf5" or (low.endswith((".h5", ".hdf5")) and fmt in (None, "aizynthfinder_hdf5")):
        if low.endswith((".h5", ".hdf5")):
            try:
                import pandas as pd
                with pd.HDFStore(path, mode="r") as hs:
                    if "/table" not in hs.keys():
                        return {"status": "invalid", "reason": "no 'table' dataset in the HDF5 file"}
                    node = hs.get_node("table/axis1")          # row index of a fixed-format frame
                    nrows = int(node.shape[0]) if node is not None else 0
                head = pd.read_hdf(path, key="table", start=0, stop=5)
                if "inchi_key" not in head.columns:
                    return {"status": "invalid", "reason": f"no 'inchi_key' column (columns: {list(head.columns)})"}
                if not nrows:
                    return {"status": "invalid", "reason": "no_records"}
                return {"status": "ok", "n_records": int(nrows), "size_bytes": os.path.getsize(path)}
            except ImportError:
                return {"status": "unchecked", "reason": "PyTables not installed; HDF5 content not verified"}
            except Exception as e:
                return {"status": "invalid", "reason": f"unreadable HDF5 stock: {type(e).__name__}: {str(e)[:120]}"}
        with open(path) as fh:                     # text stock: one InChIKey per line
            lines = [l.strip() for l in fh if l.strip()][:2000]
        good = sum(1 for l in lines if len(l) == 27 and l[14] == "-" and l[25] == "-")
        if not lines or good / len(lines) <= 0.8:
            return {"status": "invalid", "reason": "not_inchikeys", "n_sampled": len(lines)}
        return {"status": "ok", "n_sampled": len(lines)}
    if fmt == "parquet" or low.endswith((".parquet", ".gzip")):
        try:
            import pyarrow.parquet as pq
        except ImportError:
            return {"status": "unchecked", "reason": "pyarrow not installed; Parquet content not verified"}
        try:
            md = pq.ParquetFile(path)
            names = list(md.schema_arrow.names)
            nrows = md.metadata.num_rows
        except Exception as e:
            return {"status": "invalid", "reason": f"unreadable Parquet file: {type(e).__name__}: {str(e)[:120]}"}
        need = check.get("columns")
        if need:
            miss = [c for c in need if c not in names]
            if miss:
                return {"status": "invalid", "reason": f"missing columns {miss}"}
        elif not any(n.lower() == "smiles" for n in names):
            return {"status": "invalid", "reason": f"no smiles column (columns: {names[:8]})"}
        if not nrows:
            return {"status": "invalid", "reason": "no_records"}
        return {"status": "ok", "n_records": int(nrows)}
    try:
        if low.endswith((".sdf", ".sdf.gz")) or fmt == "sdf":
            n_blocks = _count_sdf_records(path)
            n, ns, nbad = _sample_smiles(path) if n_blocks else (0, 0, 0)
            if n_blocks == 0 or n == 0:
                return {"status": "invalid", "reason": "no_readable_records", "n_records": n_blocks}
            # records the reader could not parse were skipped; flag a mostly unreadable file
            if n_blocks and n / n_blocks < 0.8:
                return {"status": "invalid", "reason": f"only {n} of {n_blocks} SDF records are readable", "n_records": n_blocks}
            return {"status": "ok", "n_records": n_blocks, "sampled": ns, "parse_failures_in_sample": nbad}
        if fmt in ("smiles", "sdf_or_smiles", "smiles_tranche_files", "json_or_csv_with_smiles") or low.endswith((".smi", ".smi.gz", ".txt", ".txt.gz", ".csv", ".csv.gz", ".tsv")):
            n, ns, nbad = _sample_smiles(path)
            if n == 0:
                return {"status": "invalid", "reason": "no_records"}
            if ns and nbad / ns > 0.2:
                return {"status": "invalid", "reason": "more_than_20pct_parse_failures_in_sample", "n_records": n, "sampled": ns, "parse_failures_in_sample": nbad}
            return {"status": "ok", "n_records": n, "sampled": ns, "parse_failures_in_sample": nbad}
    except Exception as e:                   # unreadable / wrong format
        return {"status": "invalid", "reason": f"{type(e).__name__}: {str(e)[:160]}"}
    return {"status": "unchecked", "reason": f"no validator for format {fmt!r}"}


def validate(path: str, fmt: Optional[str], check: Optional[dict] = None) -> dict:
    """Light, format-aware validation. Returns {status: ok|invalid|unchecked, ...}."""
    check = check or {}
    if not os.path.exists(path):
        return {"status": "invalid", "reason": "path_not_found"}
    if os.path.isdir(path):
        names = sorted(f for f in os.listdir(path) if not f.startswith(".") and os.path.isfile(os.path.join(path, f)))
        if not names:
            return {"status": "invalid", "reason": "empty_directory"}
        if check.get("dir_requires") or check.get("dir_requires_ext"):
            miss = [f for f in check.get("dir_requires", []) if f not in names]
            if miss:
                return {"status": "invalid", "reason": f"missing files {miss}"}
            ext = check.get("dir_requires_ext")
            if ext and not any(f.lower().endswith(ext) for f in names):
                return {"status": "invalid", "reason": f"no *{ext} files"}
            return {"status": "ok", "n_files": len(names)}
        from ..purchasing.catalogue import is_catalogue_file
        rec = [f for f in names if is_catalogue_file(f) or f.lower().endswith((".h5", ".hdf5"))]
        if not rec:
            return {"status": "invalid", "reason": "no_recognised_catalogue_files", "files_found": names[:5]}
        picks = sorted({rec[0], rec[len(rec) // 2], rec[-1]})
        n_total = 0
        for f in picks:
            v = _validate_file(os.path.join(path, f), fmt, check)
            if v["status"] == "invalid":
                return {"status": "invalid", "reason": f"{f}: {v.get('reason')}", "n_files": len(rec)}
            n_total += v.get("n_records", 0)
        return {"status": "ok", "n_files": len(rec), "n_records_in_sampled_files": n_total, "sampled_files": picks}
    return _validate_file(path, fmt, check)


# ------------------------------------------------------------------ registry
from ..locking import file_lock as _locked


class Store:
    def __init__(self, root: Optional[str] = None):
        self.root = os.path.abspath(root or default_dir())
        self.path = os.path.join(self.root, "registry.json")
        self.hash_cache = os.path.join(self.root, "file_hashes.json")

    def _load(self) -> Dict[str, dict]:
        if not os.path.exists(self.path):
            return {}
        with open(self.path) as fh:
            return json.load(fh)

    def _save(self, reg: Dict[str, dict]) -> None:
        os.makedirs(self.root, exist_ok=True)
        tmp = f"{self.path}.{os.getpid()}.{uuid.uuid4().hex[:8]}.tmp"
        with open(tmp, "w") as fh:
            json.dump(reg, fh, indent=1, sort_keys=True)
        os.replace(tmp, self.path)           # atomic

    def get(self, rid: str) -> Optional[dict]:
        return self._load().get(rid)

    def all(self) -> Dict[str, dict]:
        return self._load()

    def register(self, rid: str, path: str, fmt: Optional[str] = None, hash_it: bool = True, coverage: Optional[dict] = None,
                 snapshot_date: Optional[str] = None, terms_acknowledged: bool = False, check: Optional[dict] = None) -> dict:
        path = os.path.abspath(os.path.expanduser(path))
        v = validate(path, fmt, check)
        entry = dict(path=path, format=fmt, validation=v, coverage=coverage, terms_acknowledged=bool(terms_acknowledged),
                     imported_at=dt.datetime.now().isoformat(timespec="seconds"), check=check or None)
        if os.path.exists(path):
            st = os.stat(path)
            entry["snapshot_date"] = snapshot_date or dt.datetime.fromtimestamp(st.st_mtime).date().isoformat()
            if os.path.isdir(path):
                names = [f for f in sorted(os.listdir(path)) if os.path.isfile(os.path.join(path, f)) and not f.startswith(".")]
                entry["size_bytes"] = sum(os.path.getsize(os.path.join(path, f)) for f in names)
                entry["n_files"] = len(names)
                entry["listing_sha256"] = listing_id(path)
                entry["mtime_ns"] = max([os.stat(os.path.join(path, f)).st_mtime_ns for f in names] or [0])
            else:
                entry["size_bytes"] = st.st_size
                entry["mtime_ns"] = st.st_mtime_ns
                if hash_it:
                    entry["sha256"] = file_sha256(path, self.hash_cache)
        with _locked(os.path.join(self.root, "registry.lock")):
            reg = self._load()                   # read-modify-write under the lock
            reg[rid] = entry
            self._save(reg)
        return entry

    @staticmethod
    def changed_since_import(entry: dict) -> Optional[str]:
        """Reason string when the registered file no longer matches what was imported, else None."""
        p = entry.get("path")
        if not p or not os.path.exists(p):
            return "path_not_found"
        if os.path.isdir(p):
            names = [f for f in sorted(os.listdir(p)) if os.path.isfile(os.path.join(p, f)) and not f.startswith(".")]
            size = sum(os.path.getsize(os.path.join(p, f)) for f in names)
            if entry.get("listing_sha256") and listing_id(p) != entry["listing_sha256"]:
                return "directory contents changed since import (file names or sizes differ)"
            if entry.get("size_bytes") is not None and size != entry["size_bytes"]:
                return "directory size changed since import"
            mt = max([os.stat(os.path.join(p, f)).st_mtime_ns for f in names] or [0])
            if entry.get("mtime_ns") and mt != entry["mtime_ns"]:
                return "directory files were modified since import"
            return None
        st = os.stat(p)
        if entry.get("size_bytes") is not None and st.st_size != entry["size_bytes"]:
            return "file size changed since import"
        if entry.get("mtime_ns") and st.st_mtime_ns != entry["mtime_ns"]:
            if entry.get("sha256"):
                from ..hashing import file_sha256 as _h
                if _h(p) != entry["sha256"]:
                    return "file content changed since import (sha256 differs)"
                return None                      # touched but identical
            return "file was modified since import (no hash recorded)"
        return None
