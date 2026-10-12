"""Reaction-precedent index over the Rxn-INSIGHT USPTO database (Zenodo 10.5281/zenodo.10171745, CC BY 4.0).

Needs the database file `uspto_rxn_insight.gzip` (parquet; columns REACTION, REF, TAG, TAG2, rxn_str_patt_fp). The index
holds the 2 x 167 MACCS-key fingerprints as a uint8 matrix, the TAG2 groups and an identity key per record (md5 of the
canonical stereo-retained reaction key, atom maps removed, sides sorted). Built once and cached next to the registry.
"""
from __future__ import annotations
import hashlib
import os
from typing import Dict, Optional
import numpy as np

FP_LEN = 334


def side_key(s: str, iso: bool = True) -> Optional[str]:
    from rdkit import Chem
    if s is None:
        return None
    if s == "":
        return ""
    m = Chem.MolFromSmiles(s)
    if m is None:
        return None
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return ".".join(sorted(Chem.MolToSmiles(m, isomericSmiles=iso).split(".")))


def rxn_key(rxn: str, iso: bool = True) -> Optional[str]:
    parts = rxn.split(">")
    kr, kp = side_key(parts[0], iso), side_key(parts[-1], iso)
    return None if kr is None or kp is None else kr + ">>" + kp


def key_hash(k: Optional[str]) -> Optional[str]:
    return None if k is None else hashlib.md5(k.encode()).hexdigest()


def decode_fp_column(strings) -> (np.ndarray, np.ndarray):
    """Digit strings -> (uint8 matrix N x 334, ok mask). Rows that are not 334 digits are zero and ok=False."""
    n = len(strings)
    M = np.zeros((n, FP_LEN), dtype=np.uint8)
    ok = np.zeros(n, dtype=bool)
    for i, s in enumerate(strings):
        if isinstance(s, str) and len(s) == FP_LEN:
            b = np.frombuffer(s.encode(), dtype=np.uint8)
            if ((b >= 48) & (b <= 57)).all():
                M[i] = b - 48
                ok[i] = True
    return M, ok


def _keys_chunk(rxns):
    return [key_hash(rxn_key(r, True)) for r in rxns]


class PrecedentIndex:
    def __init__(self, reaction, ref, tag, tag2, fp, fp_ok, key_stereo, source: str = ""):
        self.reaction, self.ref, self.tag, self.tag2 = reaction, ref, tag, tag2
        self.fp, self.fp_ok, self.key_stereo, self.source = fp, fp_ok, key_stereo, source
        self.groups: Dict[str, np.ndarray] = {}
        order = np.argsort(tag2, kind="stable")
        sorted_tags = tag2[order]
        bounds = np.flatnonzero(np.r_[True, sorted_tags[1:] != sorted_tags[:-1], True])
        for a, b in zip(bounds[:-1], bounds[1:]):
            self.groups[sorted_tags[a]] = np.sort(order[a:b])
        self.by_key: Dict[str, list] = {}
        for i, k in enumerate(key_stereo):
            if isinstance(k, str):
                self.by_key.setdefault(k, []).append(i)

    def __len__(self):
        return len(self.reaction)

    @classmethod
    def from_parquet(cls, path: str, identity_index: Optional[str] = None, n_proc: int = 1, limit: Optional[int] = None,
                     identity_cache: Optional[str] = None):
        import pandas as pd
        df = pd.read_parquet(path, columns=["REACTION", "REF", "TAG", "TAG2", "rxn_str_patt_fp"])
        if limit:
            df = df.iloc[:limit]
        fp, ok = decode_fp_column(df["rxn_str_patt_fp"].values)
        if identity_index:
            ks = pd.read_parquet(identity_index, columns=["key_stereo"]).key_stereo.values[: len(df)]
        elif identity_cache and os.path.exists(identity_cache) and not limit:
            ks = np.array(pd.read_parquet(identity_cache, columns=["key_stereo"]).key_stereo.values, dtype=object)
            if len(ks) != len(df):
                raise ValueError(f"identity cache {identity_cache} has {len(ks)} keys for {len(df)} records; delete it to rebuild")
        else:
            rx = df["REACTION"].tolist()
            if n_proc > 1:
                from multiprocessing import Pool
                chunks = [rx[i:i + 20000] for i in range(0, len(rx), 20000)]
                with Pool(n_proc) as p:
                    ks = np.array([k for c in p.map(_keys_chunk, chunks) for k in c], dtype=object)
            else:
                ks = np.array(_keys_chunk(rx), dtype=object)
            if identity_cache and not limit:
                os.makedirs(os.path.dirname(os.path.abspath(identity_cache)), exist_ok=True)
                pd.DataFrame({"key_stereo": ks}).to_parquet(identity_cache + ".tmp")
                os.replace(identity_cache + ".tmp", identity_cache)
        return cls(df["REACTION"].values, df["REF"].values, df["TAG"].values, df["TAG2"].values, fp, ok, ks, source=os.path.abspath(path))

    def group(self, tag2: str):
        idx = self.groups.get(tag2)
        return idx, (None if idx is None else self.fp[idx]), (None if idx is None else self.fp_ok[idx])
