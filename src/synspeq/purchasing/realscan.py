"""Sharded REAL membership scan with completion tracking.

A shard file starts with a header line `#SYNSPEQ-REAL-SCAN {json}` that identifies the query set (sha256), the catalogue
(name, size and sha256 of every file the shard scanned; the sha256 of the sorted listing for the whole catalogue), the shard
number, the number of shards, and the number and sha256 of the hit lines that follow, so that a truncated or edited shard file
is detected. Hit lines are `smiles<TAB>file`.

`merge_shards` marks queried molecules as NON-hits only when every shard 0..n-1 of one scan is present and all headers agree.
With a partial set of shards only the hits are kept (a hit is a hit, wherever it was found); every other query stays unknown.
"""
from __future__ import annotations
import hashlib
import json
import os
from typing import Dict, List, Optional, Set
from .catalogue import RealLookup, catalogue_files, scan_membership

FORMAT = "synspeq-real-scan-v1"
HEADER = "#SYNSPEQ-REAL-SCAN "


def read_queries(path: str) -> List[str]:
    with open(path) as fh:
        return sorted({l.strip() for l in fh if l.strip()})


def queries_sha256(queries: List[str]) -> str:
    return hashlib.sha256("\n".join(queries).encode()).hexdigest()


def catalogue_id(files: List[str]) -> str:
    """Cheap identity of a catalogue: names and sizes. Content is checked per file (`file_sha256`) in the shard headers."""
    return hashlib.sha256(json.dumps([(os.path.basename(f), os.path.getsize(f)) for f in files]).encode()).hexdigest()


def hits_sha256(hit_lines: List[str]) -> str:
    return hashlib.sha256("\n".join(sorted(hit_lines)).encode()).hexdigest()


def scan_shard(queries_file: str, real_dir: str, shard: int, n_shards: int, out: str, smiles_col: str = "smiles",
               recanonicalise: bool = False) -> dict:
    if n_shards < 1 or not (0 <= shard < n_shards):
        raise ValueError(f"shard must satisfy 0 <= shard < n_shards (got {shard}, {n_shards})")
    import rdkit
    files = catalogue_files(real_dir)
    if not files:
        raise ValueError(f"{real_dir}: no catalogue files")
    queries = read_queries(queries_file)
    mine = files[shard::n_shards]
    hits = scan_membership(set(queries), mine, smiles_col=smiles_col, already_canonical=not recanonicalise)
    from ..hashing import file_sha256
    file_map = {os.path.basename(f): [os.path.getsize(f), file_sha256(f)] for f in mine}
    lines = [f"{k}\t{f}" for k, f in sorted(hits.items())]
    meta = dict(format=FORMAT, queries_sha256=queries_sha256(queries), n_queries=len(queries), catalogue_sha256=catalogue_id(files),
                n_catalogue_files=len(files), shard=shard, n_shards=n_shards, n_files_in_shard=len(mine), n_hits=len(hits),
                hit_lines_sha256=hits_sha256(lines), files=file_map,
                recanonicalised=bool(recanonicalise), rdkit_version=rdkit.__version__)
    tmp = f"{out}.{os.getpid()}.tmp"
    with open(tmp, "w") as fh:
        fh.write(HEADER + json.dumps(meta) + "\n")
        for l in lines:
            fh.write(l + "\n")
    os.replace(tmp, out)
    return meta


def _read_shard(path: str):
    with open(path) as fh:
        first = fh.readline()
        if not first.startswith(HEADER):
            raise ValueError(f"{path}: not a SYNSPEQ shard file (no header). Re-run `synspeq real-scan`; hit lists without a "
                             "header cannot show that the scan was complete.")
        meta = json.loads(first[len(HEADER):])
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if len(lines) != meta.get("n_hits") or hits_sha256(lines) != meta.get("hit_lines_sha256"):
        raise ValueError(f"{path}: the hit list does not match its header ({len(lines)} lines, header says {meta.get('n_hits')}); "
                         "the file is truncated or was edited. Re-run `synspeq real-scan` for this shard.")
    return meta, {l.split("\t")[0].strip() for l in lines}


def merge_shards(queries_file: str, hit_files: List[str], out_dir: str) -> dict:
    queries = read_queries(queries_file)
    qsha = queries_sha256(queries)
    metas, hits = [], set()
    for f in hit_files:
        m, h = _read_shard(f)
        if m.get("format") != FORMAT:
            raise ValueError(f"{f}: unsupported shard format {m.get('format')!r}")
        if m["queries_sha256"] != qsha:
            raise ValueError(f"{f}: scanned a different query set than {queries_file}")
        metas.append(m)
        hits |= h
    if not metas:
        raise ValueError("no shard files given")
    for key in ("catalogue_sha256", "n_shards", "recanonicalised", "rdkit_version", "n_catalogue_files"):
        if len({m[key] for m in metas}) != 1:
            raise ValueError(f"shard files disagree on {key}: {sorted({str(m[key]) for m in metas})}")
    n = metas[0]["n_shards"]
    file_map: Dict[str, list] = {}
    for m in metas:
        for name, v in m.get("files", {}).items():
            if name in file_map and file_map[name] != v:
                raise ValueError(f"catalogue file {name} has different content in two shard files")
            file_map[name] = v
    for m in metas:
        if not isinstance(m.get("files"), dict) or len(m["files"]) != m.get("n_files_in_shard"):
            raise ValueError("shard is missing file content provenance")
        if not isinstance(m.get("shard"), int) or not isinstance(n, int) or n < 1 or not 0 <= m["shard"] < n:
            raise ValueError("invalid shard number or shard count")
    present = sorted({m["shard"] for m in metas})
    dup = [s for s in present if sum(1 for m in metas if m["shard"] == s) > 1]
    complete = present == list(range(n))
    if complete and len(file_map) != metas[0]["n_catalogue_files"]:
        raise ValueError("complete scan does not cover every catalogue file")
    lk = RealLookup.load(out_dir) if os.path.isdir(out_dir) else RealLookup()
    prev = {}
    mp = os.path.join(out_dir, "meta.json")
    if os.path.exists(mp):
        with open(mp) as fh:
            prev = json.load(fh)
        for scan in prev.get("scans", []):
            for key in ("recanonicalised", "rdkit_version"):
                if scan.get(key) != metas[0][key]:
                    raise ValueError(f"lookup matching policy changed ({key}); use a new --out folder")
        old_cat = {s.get("catalogue_sha256") for s in prev.get("scans", [])}
        if old_cat and metas[0]["catalogue_sha256"] not in old_cat:
            raise ValueError(f"{out_dir} already holds lookups from a different catalogue; use a new --out folder")
        for name, v in file_map.items():           # same names and sizes but other content is a different catalogue
            if name in prev.get("files", {}) and prev["files"][name] != v:
                raise ValueError(f"{out_dir} was built from another version of catalogue file {name} (content differs); "
                                 "use a new --out folder")
    unknown_hits = hits - set(queries)
    if unknown_hits:
        raise ValueError(f"{len(unknown_hits)} hit(s) are not in the query set; shard files and queries do not belong together")
    if complete:
        lk.add_scan(set(queries), hits)
    else:
        lk.hits |= hits
        lk.queried |= hits
    lk.save(out_dir)
    scan = dict(queries_sha256=qsha, catalogue_sha256=metas[0]["catalogue_sha256"], n_shards=n, shards_present=present, duplicates=dup,
                complete=complete, n_queries=len(queries), n_hits=len(hits), rdkit_version=metas[0]["rdkit_version"],
                recanonicalised=metas[0]["recanonicalised"])
    scans = [s for s in prev.get("scans", []) if s.get("queries_sha256") != qsha] + [scan]
    files_all = dict(prev.get("files", {}), **file_map)
    with open(mp, "w") as fh:
        json.dump(dict(scans=scans, files=files_all, n_queried=len(lk.queried), n_hits=len(lk.hits)), fh, indent=1)
    return scan


def check_lookup_matches_catalogue(lookup_dir: str, catalogue_path: str) -> None:
    """Refuse a lookup that was not built from the catalogue that is configured now.

    File names and sizes are compared first, then the content of every recorded shard file is checked against the SHA-256
    recorded at scan time. Hashing a full REAL copy is slow, so a file whose size and modification time are unchanged since
    it was last verified at the same path against the same recorded hash is not hashed again (kept in `<lookup>/.verified.json`).
    Unscanned shards may stay unknown; an unverified lookup never yields measured values."""
    mp = os.path.join(lookup_dir, "meta.json")
    if not os.path.exists(mp):
        raise ValueError(f"{lookup_dir} has no meta.json: it was not written by `synspeq real-merge`, so it cannot be tied to a catalogue")
    with open(mp) as fh:
        meta = json.load(fh)
    scans = meta.get("scans", [])
    if not scans:
        raise ValueError(f"{lookup_dir}/meta.json records no scan")
    files = catalogue_files(catalogue_path)
    now = catalogue_id(files)
    bad = [s for s in scans if s.get("catalogue_sha256") != now]
    if bad:
        raise ValueError(f"{lookup_dir} was built from a different catalogue than {catalogue_path} (file names or sizes differ); "
                         "re-run `real-scan` and `real-merge` for the configured catalogue")
    from ..hashing import file_sha256
    rec = meta.get("files")
    if not isinstance(rec, dict) or not rec:
        raise ValueError(f"{lookup_dir}: missing catalogue content provenance; rebuild the lookup")
    current = {os.path.basename(f): f for f in files}
    if set(rec) - set(current):
        raise ValueError("lookup refers to files absent from the configured catalogue")
    if any(s.get("complete") for s in scans) and set(rec) != set(current):
        raise ValueError("complete lookup lacks content provenance for some catalogue files")
    cache_path = os.path.join(lookup_dir, ".verified.json")
    try:
        with open(cache_path) as fh:
            cache = json.load(fh)
        if not isinstance(cache, dict):
            cache = {}
    except (OSError, ValueError):
        cache = {}
    changed = False
    for name, identity in rec.items():
        if not isinstance(identity, list) or len(identity) != 2:
            raise ValueError(f"invalid catalogue content provenance for {name}")
        st = os.stat(current[name])
        if st.st_size != identity[0]:
            raise ValueError(f"catalogue file {name} differs in content from the version this lookup was built from")
        where = os.path.realpath(current[name])
        seen = cache.get(name)
        if (isinstance(seen, list) and len(seen) == 4 and seen[0] == st.st_size and seen[1] == st.st_mtime_ns
                and seen[2] == identity[1] and seen[3] == where):
            continue
        if file_sha256(current[name]) != identity[1]:
            raise ValueError(f"catalogue file {name} differs in content from the version this lookup was built from")
        cache[name] = [st.st_size, st.st_mtime_ns, identity[1], where]
        changed = True
    if changed:
        try:
            tmp = cache_path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(cache, fh)
            os.replace(tmp, cache_path)
        except OSError:
            pass                                               # a read-only lookup folder simply re-verifies next time
