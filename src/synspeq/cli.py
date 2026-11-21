"""synspeq command line.

  synspeq data list|info|import|check|download|status
  synspeq validate <molecules.sdf|.csv|.smi> [--pocket-field F] [--first N | --last N]
  synspeq findability-replay <molecules> --routes-dir DIR [--routes-dir DIR2 ...]   (cached AiZynthFinder records: <mol_id>.json or <title>.json)
Exit codes: 0 ok; 2 a requested/core resource is missing, an input is unusable, or a measure is blocked (nothing was computed).
"""
from __future__ import annotations
import argparse
import json
import os
import sys
from typing import Optional
from . import __version__
from .config import MEASURE_NAMES
from .io import read_sdf, read_table
from .metrics import compute_findability
from .environment import rdkit_note
from .populations import make_population
from .cli_common import add_measures, config_and_store, nonneg_int, pos_int
from .cli_data import cmd_data
from .resources import load_manifest, resolve


def _scope_name(scope):
    """The measure's name for a message scope such as P2 or P3.real_needed."""
    head, dot, rest = str(scope).partition(".")
    return MEASURE_NAMES.get(head, head) + dot + rest


def _print_messages(res):
    for m in res.messages:
        print(f"[{m['level']}] {_scope_name(m['scope'])}: {m['text']}" + (f"\n        fix: {m['fix']}" if m.get("fix") else ""), file=sys.stderr)


def _note_rdkit():
    note = rdkit_note()
    if note:
        print(f"NOTE {note}", file=sys.stderr)


def _read_any(path: str, pocket_field=None, smiles_col="smiles", name_col="mol_id"):
    low = path.lower()
    if low.endswith((".sdf", ".sdf.gz")):
        return read_sdf(path)
    return read_table(path, smiles_col=smiles_col, name_col=name_col)


def _fmt_pct(x):
    return "n/a" if x is None else f"{100 * x:.1f}%"


# ------------------------------------------------------------------ validate / replay
def _population(a):
    recs = _read_any(a.input, smiles_col=a.smiles_col, name_col=a.name_col)
    return make_population(recs, os.path.basename(a.input), a.pocket_field, first=a.first, last=a.last)


def cmd_validate(a) -> int:
    acc = _population(a).accounting()
    print(json.dumps(acc, indent=1))
    return 0


def cmd_findability_replay(a) -> int:
    from .pipeline import load_outcomes
    pop = _population(a)
    outcomes = load_outcomes(pop, a.routes_dir)
    res = compute_findability(pop, outcomes)
    out = dict(synspeq_version=__version__, population=pop.accounting(), findability=res)
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(out, fh, indent=1)
    c = res["counts"]
    print(f"{pop.label}: {res['n_input_records']} records, {res['n_invalid_records']} invalid, {res['n_distinct']} distinct, duplicate rate {_fmt_pct(res['duplicate_rate'])}")
    print(f"  direct {c['direct_purchase']}  routed {c['routed']}  found-no-saved-route {c['found_no_saved_route']}  unresolved {c['unresolved']}  error {c['execution_error']}  pending {c['pending']}")
    print(f"  Solved {_fmt_pct(res['find'])}  Solved>=1 {_fmt_pct(res['find_ge1'])}  directly purchasable share of solved {_fmt_pct(res['direct_purchase_share_of_found'])}  Solved (all records) {_fmt_pct(res['find_raw'])}  complete={res['complete']}")
    return 0


def cmd_evaluate(a) -> int:
    from .pipeline import evaluate, load_outcomes
    man = load_manifest()
    cfg, store = config_and_store(a)
    res = resolve(cfg, man, store, mode="replay")
    _print_messages(res)
    if res.exit_code:
        print(f"blocked: nothing was computed (see messages above; resource store {store.root})", file=sys.stderr)
        return res.exit_code
    _note_rdkit()
    pf = a.pocket_field or cfg.get("pocket_field")
    recs = _read_any(a.input, smiles_col=a.smiles_col, name_col=a.name_col)
    pop = make_population(recs, os.path.basename(a.input), pf, first=a.first, last=a.last)
    acc = pop.accounting()
    if acc.get("n_valid_records_without_pocket"):
        print(f"NOTE {acc['n_valid_records_without_pocket']} valid record(s) have no value for {pf!r}; they are evaluated as one "
              "group without a pocket", file=sys.stderr)
    outcomes = load_outcomes(pop, a.routes_dir)
    prior = None
    if a.scaffold_prior:
        with open(a.scaffold_prior) as fh:
            prior = {l.split("\t")[0]: int(l.split("\t")[1]) for l in fh if l.strip()}
    rep = evaluate(pop, outcomes, res, a.out, store.all(), real_lookup_dir=a.real_lookup, scaffold_prior=prior, cache_dir=os.path.join(store.root, "cache"))
    for line in res.banner:
        print("NOTE", line)
    p1 = rep.get("findability")
    if p1:
        print(f"Findability  Solved {_fmt_pct(p1['find'])}  Solved>=1 {_fmt_pct(p1['find_ge1'])}  distinct {p1['n_distinct']}  complete={p1['complete']}")
    p2 = rep.get("precedent_support")
    if p2:
        print(f"Precedent support  primary route {_fmt_pct(p2['primary']['conditional'])} of routed, any saved route {_fmt_pct(p2['any_saved']['conditional'])}")
    p3 = rep.get("residual_steps")
    if p3 and p3.get("n_routed"):
        od = p3.get("needs_on_demand")
        od_txt = "n/a (needs both tiers)" if od is None else (("at least " if p3.get("needs_on_demand_status") == "bound" else "") + _fmt_pct(od))
        bnd = " (bound: some intermediates were never looked up in REAL)" if p3.get("frontier_status") == "bound" else ""
        print(f"Residual steps  median {p3['median_residual_steps']}  frontier found {_fmt_pct(p3['frontier_found'])}  on demand {od_txt}{bnd}")
    p4 = rep.get("scaffold_diversity")
    if p4:
        b = p4["routed"]
        if not b["n_molecules"]:
            print("Scaffold diversity  no routed molecules")
        elif b["n50"] is None:
            print(f"Scaffold diversity  N50 n/a  CSR-AUC n/a  (no cyclic molecule among {b['n_molecules']} routed molecules)")
        else:
            print(f"Scaffold diversity  N50 {b['n50']}  unique scaffold ratio {b['unique_scaffold_ratio']:.3f}  CSR-AUC {b['csr_auc']:.3f}")
    print(f"wrote {a.out}/report.html, report.json, molecules.csv, run_manifest.json")
    return 0


def _planner_config(res, man, store, out_dir: str):
    """AiZynthFinder config from the resolved resources: models dir (config.yml) + the planner stock, every member converted
    to an InChIKey stock when it is not one already. Returns (config_path, {resource id: stock file used})."""
    from .planners import build_config
    from .resources.stocks import compiled_planner_stock
    models = res.states["aizynth_models"]
    if models.state != "ready":
        return None, {}
    mc = os.path.join(models.path, "config.yml") if os.path.isdir(models.path) else models.path
    cache = os.path.join(store.root, "cache")
    stocks = {}
    for rid in res.planner_stocks:
        st = res.states[rid]
        stocks[rid] = compiled_planner_stock(rid, st.path, cache, os.path.join(store.root, "file_hashes.json"), man.resources[rid].flags)
    if not stocks:
        return None, {}
    return build_config(mc, stocks, os.path.join(out_dir, "aizynth_config.yml")), stocks


def _check_plan_fingerprint(out_dir: str, res, store, stock_files: dict) -> int:
    """Planner records are reused when a run is repeated. The fingerprint is made of the CONTENT of the model files and of every
    stock file actually given to the planner, the search settings and the planner version; a folder made under another
    fingerprint is refused instead of being mixed."""
    import hashlib
    from .hashing import content_id, file_sha256
    from .planners import PLANNER_SETTINGS
    from importlib import metadata
    hc = os.path.join(store.root, "file_hashes.json")
    models = res.states["aizynth_models"]
    try:
        ver = metadata.version("aizynthfinder")
    except metadata.PackageNotFoundError:
        ver = None
    parts = dict(models=content_id(models.path, hc, exclude=("zinc_stock.hdf5",)),
                 stocks={rid: file_sha256(p, hc) for rid, p in sorted(stock_files.items())},
                 settings=PLANNER_SETTINGS, aizynthfinder=ver)
    fp = hashlib.sha256(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()
    path = os.path.join(out_dir, "plan_fingerprint.json")
    if os.path.exists(path):
        old = json.load(open(path)).get("fingerprint")
        if old != fp:
            print(f"{out_dir} holds planner records made with a different stock, models or settings; use a new --out folder "
                  "(mixing them would make the results incomparable).", file=sys.stderr)
            return 2
    else:
        with open(path, "w") as fh:
            json.dump(dict(fingerprint=fp, **parts), fh, indent=1, default=str)
    return 0


def cmd_plan(a) -> int:
    """Plan every distinct molecule with AiZynthFinder (resumable; one chunk per SLURM array task)."""
    from .planners import plan_population
    man = load_manifest()
    cfg, store = config_and_store(a, dict(pillars=["P1"]))
    res = resolve(cfg, man, store, mode="plan")
    _print_messages(res)
    if res.exit_code:
        print(f"blocked: the planner needs AiZynthFinder models and a stock (see above; resource store {store.root})", file=sys.stderr)
        return res.exit_code
    _note_rdkit()
    os.makedirs(a.out, exist_ok=True)
    conf, stock_files = _planner_config(res, man, store, a.out)
    pf = a.pocket_field or cfg.get("pocket_field")
    pop = make_population(_read_any(a.input, smiles_col=a.smiles_col, name_col=a.name_col), os.path.basename(a.input), pf, first=a.first, last=a.last)
    items = [(m.mol_id, m.identity) for m in pop.molecules]
    rc = _check_plan_fingerprint(a.out, res, store, stock_files)
    if rc:
        return rc
    info = plan_population(items, conf, a.out, workers=a.workers, chunk=a.chunk, n_chunks=a.n_chunks)
    with open(os.path.join(a.out, f"plan_manifest_chunk{a.chunk}.json"), "w") as fh:
        json.dump(dict(synspeq_version=__version__, config=conf, stocks=stock_files, population=pop.accounting(), **info), fh, indent=1)
    print(f"planned {info['planned']} (skipped {info['skipped']}, errors {info['errors']}) of {info['n_molecules']} in {info['seconds']} s -> {a.out}")
    return 0


def cmd_run(a) -> int:
    """plan (resumable) then evaluate, in one command."""
    rc = cmd_plan(argparse.Namespace(**{**vars(a), "out": os.path.join(a.out, "routes"), "chunk": 0, "n_chunks": 1}))
    if rc:
        return rc
    ev = argparse.Namespace(**{**vars(a), "routes_dir": [os.path.join(a.out, "routes")], "out": os.path.join(a.out, "evaluation"),
                               "profile": a.profile if hasattr(a, "profile") else None, "real_lookup": getattr(a, "real_lookup", None),
                               "scaffold_prior": getattr(a, "scaffold_prior", None), "pillars": getattr(a, "pillars", None)})
    return cmd_evaluate(ev)


def cmd_report(a) -> int:
    from .reporting import render_html
    with open(os.path.join(a.dir, "report.json")) as fh:
        rep = json.load(fh)
    out = render_html(rep, os.path.join(a.dir, "report.html"))
    print(f"wrote {out}")
    return 0


def cmd_real_scan(a) -> int:
    from .purchasing.realscan import scan_shard
    m = scan_shard(a.queries, a.real_dir, a.shard, a.n_shards, a.out, smiles_col=a.smiles_col, recanonicalise=a.recanonicalise)
    print(f"shard {m['shard']}/{m['n_shards']}: {m['n_files_in_shard']} files scanned, {m['n_hits']}/{m['n_queries']} queries found -> {a.out}")
    return 0


def cmd_real_merge(a) -> int:
    from .purchasing.realscan import merge_shards
    s = merge_shards(a.queries, a.hits, a.out)
    if s["complete"]:
        print(f"scan complete ({s['n_shards']} of {s['n_shards']} shards): {s['n_hits']} hits; the other {s['n_queries'] - s['n_hits']} "
              f"queries are recorded as non-hits -> {a.out}")
    else:
        missing = sorted(set(range(s["n_shards"])) - set(s["shards_present"]))
        print(f"scan INCOMPLETE ({len(s['shards_present'])} of {s['n_shards']} shards; missing {missing[:20]}): {s['n_hits']} hits kept, "
              f"all other queries stay UNKNOWN -> {a.out}", file=sys.stderr)
    return 0


def cmd_prior(a) -> int:
    """Scaffold prior from a random draw of rows from random catalog files (default: 3000 parts x 100 rows)."""
    import random
    from .coverage import scaffold_prior
    from .purchasing.catalogue import _iter_smiles, catalogue_files
    rng = random.Random(a.seed)
    files = catalogue_files(a.real_dir)
    pick = rng.sample(files, min(a.n_parts, len(files)))
    smis = []
    for f in pick:
        col = list(_iter_smiles(f, a.smiles_col))
        smis += rng.sample(col, min(a.rows_per_part, len(col)))
    pr = scaffold_prior(smis)
    with open(a.out, "w") as fh:
        for k, v in sorted(pr.items(), key=lambda kv: -kv[1]):
            fh.write(f"{k}\t{v}\n")
    print(f"seed {a.seed}: {len(pick)} files, {len(smis)} molecules, {len(pr)} distinct scaffolds -> {a.out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="synspeq", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"synspeq {__version__}")
    ap.add_argument("--data-dir", default=None, help="resource store (default: $SYNSPEQ_DATA_DIR, else ./.synspeq if it exists, else ~/.local/share/synspeq)")
    sp = ap.add_subparsers(dest="cmd", required=True)

    d = sp.add_parser("data", help="resources: list, info, import, check, download")
    ds = d.add_subparsers(dest="data_cmd", required=True)
    ds.add_parser("list")
    ds.add_parser("status")
    x = ds.add_parser("info"); x.add_argument("id")
    x = ds.add_parser("import"); x.add_argument("id"); x.add_argument("path")
    x.add_argument("--no-hash", action="store_true", help="skip the sha256 (very large files)")
    x.add_argument("--snapshot-date", default=None)
    x.add_argument("--accept-terms", action="store_true")
    x = ds.add_parser("download"); x.add_argument("id"); x.add_argument("--dest"); x.add_argument("--accept-terms", action="store_true")
    x = ds.add_parser("build-union", help="merge registered catalogs into one planner stock")
    x.add_argument("--resources", required=True, help="comma list of registered resource ids")
    x.add_argument("--out", required=True)
    x = ds.add_parser("check")
    x.add_argument("--config"); x.add_argument("--profile", default=None); x.add_argument("--on-missing", default=None, choices=["stop", "skip"])
    add_measures(x)
    x.add_argument("--resources", default=None, help="comma list of resource ids to enable")
    x.add_argument("--planner-stock", default=None)
    x.add_argument("--mode", default="plan", choices=["plan", "replay"])
    d.set_defaults(func=cmd_data)

    def common(p):
        p.add_argument("input")
        p.add_argument("--pocket-field", default=None)
        p.add_argument("--smiles-col", default="smiles")
        p.add_argument("--name-col", default="mol_id")
        g = p.add_mutually_exclusive_group()
        g.add_argument("--first", type=pos_int, default=None, help="keep the first N rows (per pocket); N >= 1")
        g.add_argument("--last", type=pos_int, default=None, help="keep the last N rows (per pocket); N >= 1")
    v = sp.add_parser("validate", help="ingestion audit: valid/invalid records, raw and distinct counts")
    common(v); v.set_defaults(func=cmd_validate)
    r = sp.add_parser("findability-replay", aliases=["p1-replay"], help="findability from cached AiZynthFinder records")
    common(r); r.add_argument("--routes-dir", action="append", required=True); r.add_argument("--out", default=None)
    r.set_defaults(func=cmd_findability_replay)
    e = sp.add_parser("evaluate", help="all assessable measures from cached planner records (replay mode)")
    common(e)
    e.add_argument("--routes-dir", action="append", required=True)
    e.add_argument("--out", required=True)
    e.add_argument("--config"); e.add_argument("--profile", default=None); e.add_argument("--on-missing", default=None, choices=["stop", "skip"])
    add_measures(e)
    e.add_argument("--resources", default=None, help="comma list of resource ids to enable, e.g. enamine_bb,enamine_real")
    e.add_argument("--planner-stock", default=None)
    e.add_argument("--real-lookup", default=None, help="folder with hits.txt / queried.txt from `real-scan` + `real-merge`")
    e.add_argument("--scaffold-prior", default=None, help="TSV scaffold<TAB>count from a reference catalog sample (reference affinity of scaffold diversity)")
    e.set_defaults(func=cmd_evaluate)
    pl = sp.add_parser("plan", help="plan with AiZynthFinder (needs the 'planner' extras, models and a stock)")
    common(pl)
    pl.add_argument("--out", required=True)
    pl.add_argument("--config"); pl.add_argument("--on-missing", default=None, choices=["stop", "skip"])
    pl.add_argument("--resources", default=None)
    pl.add_argument("--workers", type=pos_int, default=1)
    pl.add_argument("--chunk", type=nonneg_int, default=0); pl.add_argument("--n-chunks", type=pos_int, default=1)
    pl.add_argument("--planner-stock", default=None, help="comma list of resources that form the planner stock (e.g. enamine_bb)")
    pl.set_defaults(func=cmd_plan)
    rn = sp.add_parser("run", help="plan, then evaluate (SDF to report in one command)")
    common(rn)
    rn.add_argument("--out", required=True)
    rn.add_argument("--config"); rn.add_argument("--profile", default=None); rn.add_argument("--on-missing", default=None, choices=["stop", "skip"])
    add_measures(rn); rn.add_argument("--resources", default=None)
    rn.add_argument("--workers", type=pos_int, default=1)
    rn.add_argument("--planner-stock", default=None, help="comma list of resources that form the planner stock (e.g. enamine_bb)")
    rn.add_argument("--real-lookup", default=None); rn.add_argument("--scaffold-prior", default=None)
    rn.set_defaults(func=cmd_run)
    x = sp.add_parser("report", help="re-render report.html from an evaluation folder's report.json")
    x.add_argument("dir"); x.set_defaults(func=cmd_report)
    x = sp.add_parser("real-scan", help="one shard of the REAL membership scan for the queries written by evaluate")
    x.add_argument("--queries", required=True); x.add_argument("--real-dir", required=True)
    x.add_argument("--shard", type=nonneg_int, default=0); x.add_argument("--n-shards", type=pos_int, default=1)
    x.add_argument("--smiles-col", default="smiles"); x.add_argument("--out", required=True)
    x.add_argument("--recanonicalise", action="store_true",
                   help="canonicalise every catalog SMILES with this RDKit before matching (slow; use when the catalog was "
                        "canonicalised by a different RDKit version)")
    x.set_defaults(func=cmd_real_scan)
    x = sp.add_parser("scaffold-prior", help="reference scaffold prior for scaffold diversity: random draw of rows from random catalog files (REAL)")
    x.add_argument("--real-dir", required=True); x.add_argument("--out", required=True)
    x.add_argument("--seed", type=int, required=True); x.add_argument("--n-parts", type=int, default=3000)
    x.add_argument("--rows-per-part", type=int, default=100); x.add_argument("--smiles-col", default="smiles")
    x.set_defaults(func=cmd_prior)
    x = sp.add_parser("real-merge", help="merge shard hit files into a REAL lookup folder")
    x.add_argument("--queries", required=True); x.add_argument("--hits", nargs="+", required=True); x.add_argument("--out", required=True)
    x.set_defaults(func=cmd_real_merge)
    return ap


def main(argv: Optional[list] = None) -> int:
    a = build_parser().parse_args(argv)
    try:
        return a.func(a)
    except (OSError, EOFError, ValueError) as e:      # unusable input or configuration: a message, not a traceback
        print(f"synspeq: error: {type(e).__name__}: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
