"""Input record schema. Record identity is separate from molecular identity (identity.py) and from the
evaluation unit (populations.py); the user-supplied name may repeat and is never used as the key."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class InputRecord:
    record_id: str                 # "<source>#<index>", unique inside one run
    source: str                    # file name or label
    index: int                     # 0-based position in the source
    name: Optional[str]            # SDF title line / name column; may be empty or repeated
    props: dict = field(default_factory=dict)
    identity: Optional[str] = None  # canonical isomeric SMILES; None for invalid records
    status: str = "ok"             # ok | invalid
    reason: Optional[str] = None   # why invalid (parse error, sanitisation failure, empty, ...)
    n_fragments: int = 1
