"""Content identity of resource files. A (size, mtime) pair is only used to avoid re-hashing an unchanged file; the
identity itself is always the sha256 of the bytes."""
from __future__ import annotations
import hashlib
import json
import os
from typing import Dict, Iterable, Optional

_CONTENT_EXT = (".onnx", ".gz", ".yml", ".yaml", ".csv", ".json", ".txt", ".smi", ".sdf", ".hdf5", ".h5", ".parquet", ".gzip", ".pkl")


def _load(cache_file: Optional[str]) -> Dict[str, dict]:
    if not cache_file or not os.path.exists(cache_file):
        return {}
    try:
        with open(cache_file) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def _save(cache_file: Optional[str], d: Dict[str, dict]) -> None:
    if not cache_file:
        return
    os.makedirs(os.path.dirname(os.path.abspath(cache_file)), exist_ok=True)
    tmp = f"{cache_file}.{os.getpid()}.tmp"
    with open(tmp, "w") as fh:
        json.dump(d, fh)
    os.replace(tmp, cache_file)


def file_sha256(path: str, cache_file: Optional[str] = None) -> str:
    path = os.path.abspath(path)
    st = os.stat(path)
    cache = _load(cache_file)
    e = cache.get(path)
    if e and e.get("size") == st.st_size and e.get("mtime_ns") == st.st_mtime_ns:
        return e["sha256"]
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    cache[path] = dict(size=st.st_size, mtime_ns=st.st_mtime_ns, sha256=h.hexdigest())
    _save(cache_file, cache)
    return h.hexdigest()


def content_id(path: str, cache_file: Optional[str] = None, exclude: Iterable[str] = ()) -> str:
    """sha256 of a file, or of the sorted (name, sha256) pairs of the top-level data files of a directory."""
    if os.path.isfile(path):
        return file_sha256(path, cache_file)
    skip = set(exclude)
    parts = []
    for f in sorted(os.listdir(path)):
        fp = os.path.join(path, f)
        if f in skip or f.startswith(".") or not os.path.isfile(fp) or not f.lower().endswith(_CONTENT_EXT):
            continue
        parts.append((f, file_sha256(fp, cache_file)))
    return hashlib.sha256(json.dumps(parts).encode()).hexdigest()


def listing_id(path: str) -> str:
    """Cheap identity of a large directory (names and sizes only; contents are not read)."""
    parts = [(f, os.path.getsize(os.path.join(path, f))) for f in sorted(os.listdir(path))
             if os.path.isfile(os.path.join(path, f)) and not f.startswith(".")]
    return hashlib.sha256(json.dumps(parts).encode()).hexdigest()
