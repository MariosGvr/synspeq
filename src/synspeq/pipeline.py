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


def _safe_component(name: Optional[str]) -> bool:
    """A record title may be used as a file name only if it cannot leave the records folder."""
    return bool(name) and len(name) <= 200 and name not in (".", "..") and not any(c in name for c in "/\\\0") and name == name.strip()


def _contained(root: str, candidate: str) -> bool:
    r = os.path.realpath(root)
    c = os.path.realpath(candidate)
    return c == r or c.startswith(r + os.sep)


def load_outcomes(pop: Population, routes_dirs, use_titles: bool = True) -> Dict[str, object]:
    """Planner record for every molecule of the population.

    Lookup order: `<mol_id>.json` (records written by `synspeq plan`), then, when `use_titles` is true, `<title>.json` where
    the title is the source record's name (records from an earlier AiZynthFinder run named after the molecule). A record found
    by title is accepted only if (a) the title is a plain file name inside the records folder and (b) the molecule named in
    the record is the molecule being evaluated; a mismatch is an execution_error, not a result for another molecule.
    """
    out: Dict[str, object] = {}
    for m in pop.molecules:
        cands = [(m.mol_id, False)]
        name = m.source_record.name
        if use_titles and _safe_component(name):
            cands.append((name, True))
        path = None
        by_title = False
        for d in routes_dirs:
            for c, is_title in cands:
                f = os.path.join(d, f"{c}.json")
                if os.path.isfile(f) and _contained(d, f):
                    path, by_title = f, is_title
                    break
            if path:
                break
        out[m.mol_id] = load_outcome(path, expected_identity=m.identity)
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
        feats = _p2_features(reactions, st.path, out_dir, cache_dir)
        p2 = compute_p2(pop, outcomes, feats)
        for r in p2.pop("molecules"):
            mol_rows[r["mol_id"]].update({f"p2_{k}": v for k, v in r.items() if k not in ("mol_id", "pocket")})
        report["P2"] = p2
    # ---- P3
    p3s = res.metrics["P3"]
    if p3s.status in ("value", "bound") and p3s.sub:
        bb = None
        bb_stats: dict = {}
        if res.states["enamine_bb"].state in ("ready", "partial"):
            bb = load_building_blocks(res.states["enamine_bb"].path, stats=bb_stats)
            report["catalogue_stats"] = {"enamine_bb": bb_stats}
        real = None
        if real_lookup_dir and not os.path.isdir(real_lookup_dir):
            raise ValueError(f"--real-lookup {real_lookup_dir!r} is not a folder (create it with `synspeq real-merge`)")
        if res.states["enamine_real"].state in ("ready", "partial"):
            real = RealLookup.load(real_lookup_dir) if real_lookup_dir else RealLookup()
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
        # the measured lookup state refines what the resource check could know in advance
        m = res.metrics["P3"]
        if p3.get("needs_on_demand_status") == "bound":
            m.sub["real_needed"] = dict(status="bound", reason="some intermediates were never looked up in REAL: lower bound")
        if p3.get("frontier_status") == "bound" and m.status == "value":
            m.status = "bound"
            m.reason = "some intermediates were never looked up in REAL: residual steps are upper bounds"
        report["status"]["P3"] = dict(status=m.status, reason=m.reason, sub=m.sub)
        report["banner"] = res.banner
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


def _p2_features(reactions, db_path: str, out_dir: str, cache_dir: Optional[str] = None) -> dict:
    """Retrieval features for every reaction. Cached rows are reused only when the database content (sha256), the RDKit,
    Rxn-INSIGHT and RXNMapper versions and the feature schema are the same; the index and the models are built only if some
    reaction is not cached."""
    import hashlib
    from importlib import metadata
    import rdkit
    from .hashing import content_id
    from .support.engine import FeatureCache
    cache_dir = cache_dir or out_dir
    os.makedirs(cache_dir, exist_ok=True)
    def ver(pkg):
        try:
            return metadata.version(pkg)
        except metadata.PackageNotFoundError:
            return None
    cid = content_id(db_path, os.path.join(cache_dir, "file_hashes.json"))
    meta = dict(db_sha256=cid, rdkit=rdkit.__version__, key_policy="rxn-key-v1", rxn_insight=ver("rxn-insight"), rxnmapper=ver("rxnmapper"),
                features="features-v1")
    key = hashlib.sha256(json.dumps(meta, sort_keys=True).encode()).hexdigest()
    cache = FeatureCache(os.path.join(out_dir, "p2_reaction_features.jsonl"), key)
    todo = [r for r in reactions if r not in cache.rows]
    if todo:
        from .support.engine import Preprocessor, Searcher, reaction_features
        from .support.index import PrecedentIndex
        ic = os.path.join(cache_dir, f"uspto_identity_{cid[:16]}.parquet")
        index = PrecedentIndex.from_parquet(db_path, n_proc=4, identity_cache=ic, cache_meta=dict(db_sha256=cid, rdkit=rdkit.__version__, key_policy="rxn-key-v1"))
        prep, searcher = Preprocessor(), Searcher(index)
        for rxn in todo:
            cache.get_or_compute(rxn, lambda r: reaction_features(r, prep, searcher))
    return {rxn: cache.rows[rxn] for rxn in reactions}
