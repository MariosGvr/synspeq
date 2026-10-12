"""One evaluation: population + cached planner records -> the four pillars, honouring the resource policy."""
from __future__ import annotations
import csv
import json
import os
import platform
from typing import Dict, Optional
from . import __version__
from .identity import POLICY_ID
from .metrics import compute_p1, compute_p2, compute_p3, compute_p4, required_reactions
from .populations import Population
from .purchasing import Catalogue, RealLookup, frontier_queries, load_building_blocks
from .resources.policy import Resolution
from .routes import load_outcome


def load_outcomes(pop: Population, routes_dirs) -> Dict[str, object]:
    out = {}
    for m in pop.molecules:
        name = m.source_record.name
        cands = [m.mol_id] + ([name] if name else [])
        path = next((os.path.join(d, f"{c}.json") for d in routes_dirs for c in cands if os.path.exists(os.path.join(d, f"{c}.json"))), None)
        out[m.mol_id] = load_outcome(path)
    return out


def _status(res: Resolution, p: str):
    m = res.metrics[p]
    return m.status


def evaluate(pop: Population, outcomes: Dict[str, object], res: Resolution, out_dir: str, store_entries: dict,
             real_lookup_dir: Optional[str] = None, scaffold_prior: Optional[dict] = None, n_sub: int = 3000,
             cache_dir: Optional[str] = None) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    report: dict = dict(synspeq_version=__version__, identity_policy=POLICY_ID, population=pop.accounting(),
                        resources={k: dict(state=s.state, path=(os.path.basename(s.path.rstrip(os.sep)) if s.path else None), requested_by=s.requested_by, **{x: store_entries.get(k, {}).get(x) for x in ("sha256", "size_bytes", "snapshot_date", "listing_sha256")})
                                   for k, s in res.states.items() if s.state != "not_requested"},
                        status={p: dict(status=m.status, reason=m.reason, sub=m.sub) for p, m in res.metrics.items()}, banner=res.banner)
    mol_rows: Dict[str, dict] = {m.mol_id: dict(mol_id=m.mol_id, pocket=m.pocket, identity=m.identity, n_records=m.n_records) for m in pop.molecules}

    # ---- P1
    if _status(res, "P1") in ("value", "bound"):
        p1 = compute_p1(pop, outcomes)
        report["P1"] = p1
        for mid, o in outcomes.items():
            mol_rows[mid].update(outcome=o.outcome, n_steps=o.n_steps_reported, primary_route_hash=o.primary_route_hash)
    # ---- P2
    if _status(res, "P2") == "value":
        st = res.states["uspto_rxninsight_db"]
        reactions = required_reactions(pop, outcomes)
        feats = _p2_features(reactions, st.path, store_entries.get("uspto_rxninsight_db", {}), out_dir, cache_dir)
        p2 = compute_p2(pop, outcomes, feats)
        for r in p2.pop("molecules"):
            mol_rows[r["mol_id"]].update({f"p2_{k}": v for k, v in r.items() if k not in ("mol_id", "pocket")})
        report["P2"] = p2
    # ---- P3
    p3s = res.metrics["P3"]
    if p3s.status in ("value", "bound") and p3s.sub:
        bb = None
        if res.states["enamine_bb"].state in ("ready", "partial"):
            bb = load_building_blocks(res.states["enamine_bb"].path)
        real = None
        if res.states["enamine_real"].state in ("ready", "partial"):
            real = RealLookup.load(real_lookup_dir) if real_lookup_dir and os.path.isdir(real_lookup_dir) else RealLookup()
        cat = Catalogue(bb, real)
        if real is not None:
            q = frontier_queries(outcomes, cat)
            if q:
                with open(os.path.join(out_dir, "real_queries.smi"), "w") as fh:
                    fh.write("\n".join(sorted(q)) + "\n")
                report["real_queries_pending"] = len(q)
        p3 = compute_p3(pop, outcomes, cat)
        for r in p3.pop("rows"):
            mol_rows[r["mol_id"]].update({f"p3_{k}": v for k, v in r.items() if k not in ("mol_id", "pocket")})
        report["P3"] = p3
    # ---- P4
    if _status(res, "P4") in ("value", "bound"):
        report["P4"] = compute_p4(pop, outcomes, scaffold_prior, n_sub)

    with open(os.path.join(out_dir, "report.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    from .reporting import render_html
    render_html(report, os.path.join(out_dir, "report.html"))
    cols = sorted({k for r in mol_rows.values() for k in r}, key=lambda c: (c not in ("mol_id", "pocket", "identity", "n_records", "outcome"), c))
    with open(os.path.join(out_dir, "molecules.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in mol_rows.values():
            w.writerow(r)
    with open(os.path.join(out_dir, "run_manifest.json"), "w") as fh:
        json.dump(dict(synspeq_version=__version__, python=platform.python_version(), resources=report["resources"], banner=res.banner), fh, indent=1, default=str)
    return report


def _p2_features(reactions, db_path: str, entry: dict, out_dir: str, cache_dir: Optional[str] = None) -> dict:
    from .support.engine import FeatureCache, Preprocessor, Searcher, reaction_features
    from .support.index import PrecedentIndex
    ic = os.path.join(cache_dir, "uspto_identity_keys.parquet") if cache_dir else None
    index = PrecedentIndex.from_parquet(db_path, n_proc=4, identity_cache=ic)
    prep, searcher = Preprocessor(), Searcher(index)
    tag = entry.get("sha256") or entry.get("snapshot_date") or os.path.basename(db_path)
    cache = FeatureCache(os.path.join(out_dir, "p2_reaction_features.jsonl"), tag)
    return {rxn: cache.get_or_compute(rxn, lambda r: reaction_features(r, prep, searcher)) for rxn in reactions}
