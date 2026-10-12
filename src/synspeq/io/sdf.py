"""SDF / SDF.GZ ingestion. Records are split on '$$$$' so the position of every record, readable or not,
is known; failures are returned with a reason instead of being dropped or repaired."""
from __future__ import annotations
import gzip
import os
from typing import Iterator
from rdkit import Chem
from ..identity import identity_from_mol, n_fragments
from ..records import InputRecord

_SANITIZE_FLAGS = [(getattr(Chem.SanitizeFlags, n), n) for n in dir(Chem.SanitizeFlags) if n.startswith("SANITIZE_") and n not in ("SANITIZE_NONE", "SANITIZE_ALL")]


def _open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.endswith(".gz") else open(path, "rt", encoding="utf-8", errors="replace")


def _blocks(fh) -> Iterator[str]:
    buf = []
    for line in fh:
        if line.startswith("$$$$"):
            yield "".join(buf)
            buf = []
        else:
            buf.append(line)
    if "".join(buf).strip():
        yield "".join(buf)          # final record without a terminator


def iter_sdf_blocks(path: str) -> Iterator[str]:
    """Raw text of every SDF record (also a last record without '$$$$'). Shared by the molecule reader, the catalogue
    reader and resource validation so that all three see the same records."""
    with _open_text(path) as fh:
        yield from _blocks(fh)


def molblock_of(block: str):
    """The molfile part of a record (up to and including 'M  END') or None."""
    mb, _tail, _title = _split(block)
    return mb


def _split(block: str):
    lines = block.splitlines()
    for i, l in enumerate(lines):
        if l.startswith("M  END"):
            return "\n".join(lines[: i + 1]) + "\n", lines[i + 1:], (lines[0].strip() if lines else "")
    return None, lines, (lines[0].strip() if lines else "")


def _props(tail):
    out, key = {}, None
    for l in tail:
        if l.startswith(">") and "<" in l and ">" in l[l.index("<"):]:
            key = l[l.index("<") + 1: l.rindex(">")]
            out[key] = []
        elif key is not None:
            if l.strip() == "":
                key = None
            else:
                out[key].append(l)
    return {k: "\n".join(v) for k, v in out.items()}


def _sanitize_failure(mol) -> str:
    probe = Chem.Mol(mol)
    res = Chem.SanitizeMol(probe, catchErrors=True)
    for flag, name in _SANITIZE_FLAGS:
        if res & flag:
            return "sanitize_failed:" + name.replace("SANITIZE_", "").lower()
    return "sanitize_failed"


def read_sdf(path: str, source: str | None = None) -> Iterator[InputRecord]:
    src = source or os.path.basename(path)
    for i, block in enumerate(iter_sdf_blocks(path)):
        molblock, tail, title = _split(block)
        rid = f"{src}#{i}"
        props = _props(tail)
        if molblock is None:
            yield InputRecord(rid, src, i, title or None, props, None, "invalid", "no_molblock_end")
            continue
        mol = Chem.MolFromMolBlock(molblock, sanitize=False, removeHs=False)
        if mol is None:
            yield InputRecord(rid, src, i, title or None, props, None, "invalid", "molblock_unreadable")
            continue
        if mol.GetNumAtoms() == 0:
            yield InputRecord(rid, src, i, title or None, props, None, "invalid", "no_atoms")
            continue
        reason = None
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveHs(mol)
            if mol.GetNumConformers() and mol.GetConformer().Is3D():
                Chem.AssignStereochemistryFrom3D(mol)
        except Exception:                       # RDKit raises on sanitisation failure
            reason = _sanitize_failure(mol)
        if reason:
            yield InputRecord(rid, src, i, title or None, props, None, "invalid", reason)
        else:
            ident = identity_from_mol(mol)
            yield InputRecord(rid, src, i, title or None, props, ident, "ok", None, n_fragments(ident))
