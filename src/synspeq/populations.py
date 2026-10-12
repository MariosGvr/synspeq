"""Populations: record -> molecule mapping, windows, raw and distinct accounting.

Evaluation unit ("molecule"): (pocket, identity) when a pocket field is given, otherwise identity alone, so
identical molecules in different pockets stay separate and a single ungrouped collection is deduplicated
globally. The source record of a molecule is its first occurrence in input order (first-occurrence rule); duplicates keep their multiplicity so raw and distinct views share one set of outcomes.
"""
from __future__ import annotations
import hashlib
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional
from .records import InputRecord


@dataclass
class Molecule:
    mol_id: str
    identity: str
    pocket: Optional[str]
    source_record: InputRecord          # first occurrence
    n_records: int = 1                  # raw rows sharing this unit (multiplicity)


@dataclass
class Population:
    label: str
    records: List[InputRecord]
    pocket_field: Optional[str] = None
    window: Optional[str] = None        # human-readable description of the applied window
    molecules: List[Molecule] = field(default_factory=list)
    record_to_mol: Dict[str, str] = field(default_factory=dict)

    def pocket_of(self, r: InputRecord) -> Optional[str]:
        if self.pocket_field is None:
            return None
        v = r.props.get(self.pocket_field)
        return None if v is None or str(v).strip() == "" else str(v)

    def build(self) -> "Population":
        units: Dict[tuple, Molecule] = {}
        self.molecules, self.record_to_mol = [], {}
        for r in self.records:
            if r.status != "ok":
                continue
            p = self.pocket_of(r)
            k = (p, r.identity)
            m = units.get(k)
            if m is None:
                mid = "m_" + hashlib.sha1(f"{p}|{r.identity}".encode()).hexdigest()[:14]
                m = Molecule(mid, r.identity, p, r, 1)
                units[k] = m
                self.molecules.append(m)
            else:
                m.n_records += 1
            self.record_to_mol[r.record_id] = m.mol_id
        return self

    def accounting(self) -> dict:
        n_in = len(self.records)
        valid = sum(1 for r in self.records if r.status == "ok")
        reasons: Dict[str, int] = {}
        for r in self.records:
            if r.status != "ok":
                reasons[r.reason or "unknown"] = reasons.get(r.reason or "unknown", 0) + 1
        distinct = len(self.molecules)
        return dict(
            label=self.label, window=self.window, pocket_field=self.pocket_field,
            n_input_records=n_in, n_valid_records=valid, n_invalid_records=n_in - valid, invalid_reasons=reasons,
            n_distinct_molecules=distinct,
            n_duplicate_records=valid - distinct,
            duplicate_rate=(valid - distinct) / valid if valid else None,   # share of valid raw rows repeating an earlier unit
            n_multi_fragment=sum(1 for m in self.molecules if "." in m.identity),
            n_pockets=len({m.pocket for m in self.molecules if m.pocket is not None}) or None,
        )


def make_population(records: Iterable[InputRecord], label: str, pocket_field: Optional[str] = None,
                    first: Optional[int] = None, last: Optional[int] = None) -> Population:
    """first/last: keep the first/last N *input rows* (per pocket when pocket_field is given), in file order, BEFORE
    deduplication; membership of the window is preserved by the returned records."""
    recs = list(records)
    window = None
    if first is not None or last is not None:
        if first is not None and last is not None:
            raise ValueError("choose either first or last")
        n = first if first is not None else last
        by: Dict[Optional[str], List[InputRecord]] = {}
        for r in recs:
            p = None if pocket_field is None else (str(r.props.get(pocket_field)) if r.props.get(pocket_field) not in (None, "") else None)
            by.setdefault(p, []).append(r)
        keep = set()
        for rs in by.values():
            sel = rs[:n] if first is not None else rs[-n:]
            keep.update(x.record_id for x in sel)
        recs = [r for r in recs if r.record_id in keep]
        window = f"{'first' if first is not None else 'last'} {n} rows" + (" per pocket" if pocket_field else "") + ", applied before deduplication"
    return Population(label, recs, pocket_field, window).build()
