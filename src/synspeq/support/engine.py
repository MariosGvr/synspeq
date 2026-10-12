"""Query preprocessing and retrieval for Pillar 2 (needs `rxn-insight` and `rxnmapper`, i.e. extras 'support').

prepare(): Rxn-INSIGHT classification of the forward reaction string (atom mapping by RXNMapper) -> TAG, broad TAG2,
sanitised reaction, identity keys. Statuses: ok | unparseable_smiles | preprocess_error | mapping_failed_silent |
reaction_info_error. Failures are reported, never repaired; an unscored step makes a route indeterminate.
search(): similarity of the query's MACCS-key presence pattern to every database record with the same TAG2; the
'analogue' search additionally drops every database record that is identical to the query.
"""
from __future__ import annotations
import json
import time
import warnings
from typing import Dict, List, Optional
import numpy as np
from .index import PrecedentIndex, key_hash, rxn_key


def similarity_presence(q: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Tanimoto of the PRESENCE pattern of the 334 key counts (what the installed SciPy jaccard computes on integer vectors)."""
    pq = q != 0
    pm = M != 0
    inter = (pm & pq).sum(1).astype(float)
    union = (pm | pq).sum(1).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.where(union > 0, 1.0 - inter / union, 0.0)
    return 1.0 - d


class Preprocessor:
    def __init__(self):
        warnings.filterwarnings("ignore")
        try:
            import os
            os.environ.setdefault("TQDM_DISABLE", "1")
            import torch
            # RXNMapper runs on CPU; unrestricted OpenMP threads oversubscribe a scheduler-limited allocation (observed
            # ~100x slowdown on a shared node). One thread per process unless SYNSPEQ_TORCH_THREADS says otherwise.
            torch.set_num_threads(int(os.environ.get("SYNSPEQ_TORCH_THREADS", "1")))
            from importlib.resources import files
            import pandas as pd
            from rxn_insight.utils import curate_smirks
            from rxn_insight.reaction import Reaction
            from rxnmapper import RXNMapper
        except ImportError as e:
            raise RuntimeError("Pillar 2 needs the 'support' extras: pip install rxn-insight rxnmapper (see docs/DATA_RESOURCES_GUIDE.md)") from e
        self._Reaction = Reaction
        self.mapper = RXNMapper()
        self.smirks = curate_smirks(pd.read_json(files("rxn_insight.data").joinpath("smirks.json"), orient="records", lines=True))
        self.fg = pd.read_json(files("rxn_insight.data").joinpath("functional_groups.json"), orient="records", lines=True)

    def prepare(self, rxn: str) -> dict:
        from rdkit import Chem
        rec = dict(rxn=rxn, prep_status="ok", prep_error="")
        parts = rxn.split(">")
        if len(parts) < 3 or Chem.MolFromSmiles(parts[0]) is None or Chem.MolFromSmiles(parts[-1]) is None:
            rec.update(prep_status="unparseable_smiles", prep_error="RDKit cannot parse a side of the query")
            return rec
        try:
            rx = self._Reaction(rxn, rxn_mapper=self.mapper, smirks=self.smirks, fg=self.fg)
        except Exception as e:
            rec.update(prep_status="preprocess_error", prep_error=f"{type(e).__name__}: {str(e)[:300]}")
            return rec
        mapped = rx.classifier.mapped_reaction
        if ":" not in mapped:
            rec.update(prep_status="mapping_failed_silent", prep_error="RXNMapper failed; continuing unmapped is not supported")
            return rec
        try:
            rx.name = "(not computed)"
            info = rx.get_reaction_info()
        except Exception as e:
            rec.update(prep_status="reaction_info_error", prep_error=f"{type(e).__name__}: {str(e)[:300]}")
            return rec
        try:
            tag2 = rx.give_broad_tag()
        except Exception:
            tag2 = ""
        san = rx.classifier.sanitized_reaction
        rec.update(rxni_class=info.get("CLASS"), tag=info.get("TAG"), tag2=tag2, sanitized_reaction=san,
                   raw_reaction_for_fp=rx.reaction, key_asgiven=key_hash(rxn_key(rxn, True)), key_sanitized=key_hash(rxn_key(san, True)))
        return rec


def _query_fp(rxn: str) -> np.ndarray:
    from rxn_insight.utils import get_fp
    return np.asarray(get_fp(rxn, "MACCS", True))


class Searcher:
    def __init__(self, index: PrecedentIndex, k: int = 5):
        self.index, self.k = index, k

    def search(self, rec: dict, exclude_keys: Optional[set] = None) -> dict:
        r = dict(status="", n=0, best=float("nan"), top=[], n_excl=0)
        idx, M, ok = self.index.group(rec["tag2"])
        if idx is None:
            r["status"] = "no_broad_tag_records"
            return r
        r["n"] = int(len(idx))
        sims = np.where(ok, similarity_presence(_query_fp(rec["raw_reaction_for_fp"]), M), np.nan)
        keep = np.ones(len(idx), bool)
        if exclude_keys:
            rows = {x for h in exclude_keys for x in self.index.by_key.get(h, [])}
            if rows:
                ex = np.isin(idx, np.fromiter(rows, int))
                keep &= ~ex
                r["n_excl"] = int(ex.sum())
        valid = keep & ok & (sims > 0.0)
        if not valid.any():
            r["status"] = "no_positive_similarity"
            return r
        cand = np.where(valid)[0]
        order = cand[np.argsort(-sims[cand], kind="stable")][: self.k]
        r["status"] = "retrieved"
        r["best"] = float(sims[order[0]])
        r["top"] = [dict(row=int(idx[o]), sim=float(sims[o]), reaction=str(self.index.reaction[idx[o]]), ref=str(self.index.ref[idx[o]])) for o in order]
        return r


def reaction_features(rxn: str, prep: Preprocessor, searcher: Searcher) -> dict:
    rec = prep.prepare(rxn)
    row = dict(rxn=rxn, prep_status=rec["prep_status"], prep_error=rec.get("prep_error", ""), rxni_class=rec.get("rxni_class"),
               tag=rec.get("tag"), tag2=rec.get("tag2"))
    if rec["prep_status"] != "ok":
        return row
    qk = [h for h in (rec["key_asgiven"], rec["key_sanitized"]) if h]
    idr = sorted({x for h in qk for x in searcher.index.by_key.get(h, [])})
    a = searcher.search(rec)
    b = searcher.search(rec, exclude_keys=set(qk)) if idr else a
    row.update(identity_anywhere=bool(idr), n_identical_records=len(idr),
               broad_status=a["status"], broad_n_candidates=a["n"], broad_best_sim=a["best"], broad_top5=json.dumps(a["top"]),
               broad_analogue_status=b["status"], broad_analogue_best_sim=b["best"], n_excluded_identical=b["n_excl"])
    return row


class FeatureCache:
    """reaction string -> feature row, kept as JSON lines so interrupted runs resume. Rows carry the cache key (`_db`) of the
    database and software that produced them; rows with another key are ignored, never reused.

    A final line without a newline is an interrupted append: it is dropped (with a note) and the file is truncated back to
    the last complete row. A corrupt line anywhere else is an error, because silently skipping it could hide data loss."""

    def __init__(self, path: str, db_tag: str):
        import os
        self.path, self.db_tag, self.rows, self.recovered_truncated_tail = path, db_tag, {}, False
        if not os.path.exists(path):
            return
        with open(path, "rb") as fh:
            data = fh.read()
        good_end = 0
        pos = 0
        lines = data.split(b"\n")
        for i, raw in enumerate(lines):
            end = pos + len(raw) + (1 if i < len(lines) - 1 else 0)
            is_last = i == len(lines) - 1                      # the piece after the final newline (usually empty)
            if raw.strip():
                try:
                    d = json.loads(raw.decode("utf-8"))
                except (ValueError, UnicodeDecodeError):
                    if is_last:                                # no terminating newline: an interrupted append
                        self.recovered_truncated_tail = True
                        break
                    raise ValueError(f"{path}: corrupt feature record on line {i + 1}; delete the file to recompute these features")
                if d.get("_db") == db_tag:
                    self.rows[d["rxn"]] = d
            good_end = end
            pos = end
        if self.recovered_truncated_tail:
            with open(path, "r+b") as fh:
                fh.truncate(good_end)

    def get_or_compute(self, rxn: str, fn):
        if rxn not in self.rows:
            row = fn(rxn)
            row["_db"] = self.db_tag
            self.rows[rxn] = row
            line = json.dumps(row, default=lambda o: None if o != o else str(o)) + "\n"
            from ..locking import file_lock
            with file_lock(self.path + ".lock"):               # evaluators sharing an output folder append one row at a time
                with open(self.path, "ab") as fh:
                    fh.write(line.encode("utf-8"))
                    fh.flush()
        return self.rows[rxn]
