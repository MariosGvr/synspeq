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


def _iter_smiles(path: str, smiles_col: str = "smiles") -> Iterator[str]:
    low = path.lower()
    if low.endswith((".sdf", ".sdf.gz")):
        with _open(path) as fh:
            buf = []
            for line in fh:
                if line.startswith("$$$$"):
                    blk = "".join(buf)
                    buf = []
                    cut = blk.find("M  END")
                    if cut >= 0:
                        m = Chem.MolFromMolBlock(blk[: cut + 6] + "\n", sanitize=True)
                        if m is not None:
                            yield Chem.MolToSmiles(Chem.RemoveHs(m))
                else:
                    buf.append(line)
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
            for row in csv.DictReader(fh):
                if row.get(smiles_col):
                    yield row[smiles_col]
    else:                                              # .smi / .txt / .tsv: first whitespace-separated token
        with _open(path) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    yield line.split()[0]


def catalogue_files(path: str):
    if os.path.isdir(path):
        return sorted(os.path.join(path, f) for f in os.listdir(path) if not f.startswith(".") and os.path.isfile(os.path.join(path, f)))
    return [path]


def load_building_blocks(path: str, smiles_col: str = "smiles", already_canonical: bool = False) -> Set[str]:
    """Stocked catalogue -> set of canonical isomeric SMILES (records RDKit cannot parse are skipped and counted by the caller via len)."""
    out: Set[str] = set()
    for f in catalogue_files(path):
        for s in _iter_smiles(f, smiles_col):
            c = s if already_canonical else canonical(s)
            if c:
                out.add(c)
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
