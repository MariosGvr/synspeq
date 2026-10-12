"""Sharded REAL membership scan with completion tracking.

A shard file starts with a header line `#SYNSPEQ-REAL-SCAN {json}` that identifies the query set (sha256), the catalogue
(names and sizes of its files), the shard number and the number of shards, followed by `smiles<TAB>file` hit lines.

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
    return hashlib.sha256(json.dumps([(os.path.basename(f), os.path.getsize(f)) for f in files]).encode()).hexdigest()


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
    meta = dict(format=FORMAT, queries_sha256=queries_sha256(queries), n_queries=len(queries), catalogue_sha256=catalogue_id(files),
                n_catalogue_files=len(files), shard=shard, n_shards=n_shards, n_files_in_shard=len(mine), n_hits=len(hits),
                recanonicalised=bool(recanonicalise), rdkit_version=rdkit.__version__)
    tmp = f"{out}.{os.getpid()}.tmp"
    with open(tmp, "w") as fh:
        fh.write(HEADER + json.dumps(meta) + "\n")
        for k, f in sorted(hits.items()):
            fh.write(f"{k}\t{f}\n")
    os.replace(tmp, out)
    return meta


def _read_shard(path: str):
    with open(path) as fh:
        first = fh.readline()
        if not first.startswith(HEADER):
            raise ValueError(f"{path}: not a SYNSPEQ shard file (no header). Re-run `synspeq real-scan`; hit lists without a "
                             "header cannot show that the scan was complete.")
        meta = json.loads(first[len(HEADER):])
        hits = {l.split("\t")[0].strip() for l in fh if l.strip()}
    return meta, hits


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
    present = sorted({m["shard"] for m in metas})
    dup = [s for s in present if sum(1 for m in metas if m["shard"] == s) > 1]
    complete = present == list(range(n))
    lk = RealLookup.load(out_dir) if os.path.isdir(out_dir) else RealLookup()
    prev = {}
    mp = os.path.join(out_dir, "meta.json")
    if os.path.exists(mp):
        with open(mp) as fh:
            prev = json.load(fh)
        old_cat = {s.get("catalogue_sha256") for s in prev.get("scans", [])}
        if old_cat and metas[0]["catalogue_sha256"] not in old_cat:
            raise ValueError(f"{out_dir} already holds lookups from a different catalogue; use a new --out folder")
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
    with open(mp, "w") as fh:
        json.dump(dict(scans=scans, n_queried=len(lk.queried), n_hits=len(lk.hits)), fh, indent=1)
    return scan
