"""Purchasing catalogues for Pillar 3.

Matching is strict: exact RDKit canonical isomeric SMILES (stereochemistry kept; no desalting, no
normalisation). Tiers: 'in_stock' = a stocked building block catalogue (Enamine BB); 'on_demand' =
a make-on-demand catalogue (Enamine REAL). A molecule found in both is 'in_stock'.

REAL is far too large to hold in memory, so it is used through lookups: the intermediates of the routes are
collected (`frontier_queries`), scanned once against the REAL files (`scan_membership`), and the answers are
kept in a RealLookup with THREE states per queried molecule: hit, non-hit, and unknown (never queried).
Unknown is not a non-hit: results that depend on unknown molecules are reported as bounds.
"""
from __future__ import annotations
import gzip
import os
from typing import Dict, Iterable, Iterator, Optional, Set
from rdkit import Chem


def canonical(smiles: str) -> Optional[str]:
    m = Chem.MolFromSmiles(smiles) if smiles else None
    return Chem.MolToSmiles(m) if m is not None else None


def _open(path: str):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path, "rt", errors="replace")


CATALOGUE_EXT = (".smi", ".txt", ".tsv", ".csv", ".sdf", ".parquet")


def is_catalogue_file(name: str) -> bool:
    n = name.lower()
    if n.endswith(".gz"):
        n = n[:-3]
    return n.endswith(CATALOGUE_EXT)


def _iter_smiles(path: str, smiles_col: str = "smiles") -> Iterator[str]:
    low = path.lower()
    if low.endswith((".sdf", ".sdf.gz")):
        from ..io.sdf import iter_sdf_blocks, molblock_of
        for blk in iter_sdf_blocks(path):             # same record splitting as the molecule reader (incl. a last record without $$$$)
            mb = molblock_of(blk)
            if mb is None:
                continue
            m = Chem.MolFromMolBlock(mb, sanitize=True)
            if m is not None:
                yield Chem.MolToSmiles(Chem.RemoveHs(m))
    elif low.endswith(".parquet"):
        try:
            import pyarrow.parquet as pq
        except ImportError as e:                      # optional dependency, only for parquet catalogues
            raise RuntimeError("reading .parquet needs pyarrow: pip install pyarrow") from e
        for s in pq.read_table(path, columns=[smiles_col]).column(smiles_col).to_pylist():
            yield s
    elif low.endswith((".csv", ".csv.gz")):
        import csv
        with _open(path) as fh:
            rd = csv.DictReader(fh)
            col = smiles_col if smiles_col in (rd.fieldnames or []) else next((c for c in (rd.fieldnames or []) if c.lower() == "smiles"), None)
            if col is None:
                raise ValueError(f"{os.path.basename(path)}: no '{smiles_col}' column (columns: {rd.fieldnames})")
            for row in rd:
                if row.get(col):
                    yield row[col]
    else:                                              # .smi / .txt / .tsv: first whitespace-separated token
        with _open(path) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    yield line.split()[0]


def catalogue_files(path: str):
    if os.path.isdir(path):
        return sorted(os.path.join(path, f) for f in os.listdir(path)
                      if not f.startswith(".") and os.path.isfile(os.path.join(path, f)) and is_catalogue_file(f))
    return [path]


def load_building_blocks(path: str, smiles_col: str = "smiles", already_canonical: bool = False, stats: Optional[dict] = None) -> Set[str]:
    """Stocked catalogue -> set of canonical isomeric SMILES. Pass a dict as `stats` to receive n_records, n_unparseable,
    n_files. A directory with no recognised catalogue file raises ValueError (it is not an empty catalogue)."""
    out: Set[str] = set()
    files = [f for f in catalogue_files(path) if is_catalogue_file(os.path.basename(f))]
    if not files:
        raise ValueError(f"{path}: no catalogue files (*.smi, *.txt, *.csv, *.sdf, *.parquet, optionally .gz)")
    n = bad = 0
    for f in files:
        for s in _iter_smiles(f, smiles_col):
            n += 1
            c = s if already_canonical else canonical(s)
            if c:
                out.add(c)
            else:
                bad += 1
    if stats is not None:
        stats.update(n_records=n, n_unparseable=bad, n_files=len(files), n_unique=len(out))
    if not out:
        raise ValueError(f"{path}: no parseable molecules in {len(files)} file(s)")
    return out


def scan_membership(queries: Set[str], files: Iterable[str], smiles_col: str = "smiles", already_canonical: bool = True) -> Dict[str, str]:
    """One pass over catalogue files: returns {query: source file} for every query present. `queries` are canonical SMILES.
    REAL is pre-canonicalised by RDKit, so catalogue strings are compared as given unless already_canonical=False."""
    hits: Dict[str, str] = {}
    for f in files:
        for s in _iter_smiles(f, smiles_col):
            k = s if already_canonical else canonical(s)
            if k in queries and k not in hits:
                hits[k] = os.path.basename(f)
    return hits


class RealLookup:
    """Recorded REAL answers: `hits` (members) and `queried` (everything asked, hit or not)."""

    def __init__(self, hits: Optional[Set[str]] = None, queried: Optional[Set[str]] = None):
        self.hits = set(hits or ())
        self.queried = set(queried or ()) | self.hits

    def state(self, canonical_smiles: str) -> str:
        if canonical_smiles in self.hits:
            return "hit"
        return "non_hit" if canonical_smiles in self.queried else "unknown"

    def add_scan(self, queries: Set[str], hits: Iterable[str]) -> None:
        self.queried |= set(queries)
        self.hits |= set(hits)

    def save(self, path: str) -> None:
        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, "hits.txt"), "w") as fh:
            fh.write("\n".join(sorted(self.hits)) + ("\n" if self.hits else ""))
        with open(os.path.join(path, "queried.txt"), "w") as fh:
            fh.write("\n".join(sorted(self.queried)) + ("\n" if self.queried else ""))

    @classmethod
    def load(cls, path: str) -> "RealLookup":
        def rd(n):
            p = os.path.join(path, n)
            return {l.strip().split("\t")[0] for l in open(p) if l.strip()} if os.path.exists(p) else set()
        return cls(rd("hits.txt"), rd("queried.txt"))


class Catalogue:
    """tier_of(smiles) -> ('in_stock' | 'on_demand' | None, real_state) with real_state in hit | non_hit | unknown | 'n/a' (no REAL configured)."""

    def __init__(self, building_blocks: Optional[Set[str]] = None, real: Optional[RealLookup] = None):
        self.bb = building_blocks
        self.real = real

    @property
    def has_stocked(self) -> bool:
        return self.bb is not None

    @property
    def has_on_demand(self) -> bool:
        return self.real is not None

    def tier_of(self, smiles: str):
        c = canonical(smiles)
        if c is None:
            return None, "n/a"
        if self.bb is not None and c in self.bb:
            return "in_stock", "n/a"
        if self.real is not None:
            st = self.real.state(c)
            return ("on_demand" if st == "hit" else None), st
        return None, "n/a"
