"""CSV / SMILES-text adapters to the same record schema."""
from __future__ import annotations
import csv
import gzip
import os
from typing import Iterator, Optional
from ..identity import identity_from_smiles, n_fragments
from ..records import InputRecord


def _open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.endswith(".gz") else open(path, "rt", encoding="utf-8", newline="")


def read_table(path: str, smiles_col: str = "smiles", name_col: Optional[str] = "mol_id", source: str | None = None) -> Iterator[InputRecord]:
    """CSV with a header (columns kept as props) or a .smi/.txt file ('SMILES [name]' per line)."""
    src = source or os.path.basename(path)
    base = path[:-3] if path.endswith(".gz") else path
    if base.lower().endswith((".smi", ".txt")):
        with _open_text(path) as fh:
            for i, line in enumerate(l for l in fh if l.strip()):
                parts = line.split()
                smi, name = parts[0], (parts[1] if len(parts) > 1 else None)
                ident, why = identity_from_smiles(smi)
                yield InputRecord(f"{src}#{i}", src, i, name, {}, ident, "ok" if ident else "invalid", why, n_fragments(ident) if ident else 1)
        return
    with _open_text(path) as fh:
        for i, row in enumerate(csv.DictReader(fh)):
            smi = row.get(smiles_col)
            ident, why = identity_from_smiles(smi)
            name = row.get(name_col) if name_col else None
            props = {k: v for k, v in row.items() if k not in (smiles_col,)}
            yield InputRecord(f"{src}#{i}", src, i, name, props, ident, "ok" if ident else "invalid", why, n_fragments(ident) if ident else 1)
