"""synspeq command line.

  synspeq data list|info|import|check|download|status
  synspeq validate <molecules.sdf|.csv|.smi> [--pocket-field F] [--first N | --last N]
  synspeq p1-replay <molecules> --routes-dir DIR [--routes-dir DIR2 ...]   (cached AiZynthFinder records named <record name>.json)
Exit codes: 0 ok; 2 a requested/core resource is missing or a metric is blocked (nothing was computed).
"""
from __future__ import annotations
import argparse
import json
import os
import shutil
import subprocess
import sys
from typing import Optional
from . import __version__
from .config import load_config
from .io import read_sdf, read_table
from .metrics import compute_p1
from .populations import make_population
from .resources import Store, load_manifest, resolve
from .routes import load_outcome


def _read_any(path: str, pocket_field=None, smiles_col="smiles", name_col="mol_id"):
    low = path.lower()
    if low.endswith((".sdf", ".sdf.gz")):
        return read_sdf(path)
    return read_table(path, smiles_col=smiles_col, name_col=name_col)


def _fmt_pct(x):
    return "n/a" if x is None else f"{100 * x:.1f}%"


def _fetch(r, dest: str, store: Store) -> int:
    """Stream r.fetch['url'] to dest, verify md5 (and size when given), then register."""
    import hashlib
    import urllib.request
    os.makedirs(dest, exist_ok=True)
    path = os.path.join(dest, r.fetch["filename"])
    tmp = path + ".part"
    h = hashlib.md5()
    n = last = 0
    with urllib.request.urlopen(r.fetch["url"]) as resp, open(tmp, "wb") as out:
        for chunk in iter(lambda: resp.read(1 << 20), b""):
            out.write(chunk)
            h.update(chunk)
            n += len(chunk)
            if n - last >= 100 << 20:
                last = n
                print(f"{n >> 20} MB", flush=True)
    if r.fetch.get("md5") and h.hexdigest() != r.fetch["md5"]:
        os.remove(tmp)
        print(f"checksum mismatch (got {h.hexdigest()}, expected {r.fetch['md5']}); file removed", file=sys.stderr)
        return 2
    if r.fetch.get("size") and n != r.fetch["size"]:
        os.remove(tmp)
        print(f"size mismatch (got {n}, expected {r.fetch['size']}); file removed", file=sys.stderr)
        return 2
    os.replace(tmp, path)
    e = store.register(r.id, path, fmt=r.format, terms_acknowledged=True)
    print(f"downloaded and registered {r.id} -> {path} (md5 verified)")
    return 0


# ------------------------------------------------------------------ data
def cmd_data(a) -> int:
    man = load_manifest()
    store = Store(a.data_dir)
    sub = a.data_cmd
    if sub == "list":
        reg = store.all()
        print(f"{'id':22} {'role':34} {'default':8} {'access':16} {'fetch':8} registered")
        for r in man.resources.values():
            e = reg.get(r.id)
            print(f"{r.id:22} {','.join(r.roles):34} {'yes' if r.default_enabled else 'no':8} {r.access:16} {'auto' if r.automated else 'guided':8} {e['path'] if e else '-'}")
        return 0
    if sub == "info":
        r = man.resources.get(a.id)
        if r is None:
            print(f"unknown resource {a.id!r}; see `synspeq data list`", file=sys.stderr)
            return 2
        for k, v in (("name", r.name), ("role", ", ".join(r.roles)), ("provides", ", ".join(r.provides)), ("access", r.access),
                     ("how to get it", r.get), ("url", r.url), ("format", r.format), ("cite", r.cite), ("terms", r.terms), ("warning", r.warn),
                     ("if missing", json.dumps(r.if_missing))):
            if v:
                print(f"{k:14}: {v}")
        return 0
    if sub == "import":
        r = man.resources.get(a.id)
        if r is None:
            print(f"unknown resource {a.id!r}", file=sys.stderr)
            return 2
        if r.terms:
            print(f"Terms for {r.id}: {r.terms}")
            if r.warn:
                print(f"Note: {r.warn}")
        e = store.register(a.id, a.path, fmt=r.format, hash_it=not a.no_hash, snapshot_date=a.snapshot_date, terms_acknowledged=a.accept_terms)
        v = e["validation"]
        print(f"registered {a.id} -> {e['path']}  validation: {v.get('status')}" + (f" ({v.get('reason')})" if v.get("reason") else ""))
        return 0 if v.get("status") in ("ok", "unchecked") else 2
    if sub == "status":
        print(json.dumps(store.all(), indent=1))
        return 0
    if sub == "check":
        cfg = load_config(a.config, dict(on_missing=a.on_missing, profile=a.profile, pillars=a.pillars.split(",") if a.pillars else None,
                                        resources=a.resources.split(",") if a.resources else None))
        res = resolve(cfg, man, store, mode=a.mode)
        print("resources")
        for s in res.states.values():
            if s.state != "not_requested":
                print(f"  {s.id:22} {s.state:10} requested by {s.requested_by:8} {s.path or ''} {('- ' + s.reason) if s.reason else ''}")
        print("metrics")
        for p, m in res.metrics.items():
            print(f"  {p}: {m.status}" + (f"  ({m.reason})" if m.reason else ""))
            for k, d in m.sub.items():
                print(f"      {k}: {d['status']}" + (f"  ({d['reason']})" if d.get("reason") else ""))
        for m in res.messages:
            print(f"[{m['level']}] {m['scope']}: {m['text']}" + (f"\n        fix: {m['fix']}" if m.get("fix") else ""))
        return res.exit_code
    if sub == "build-union":
        from .resources.union import build_union
        ids = [x for x in a.resources.split(",") if x]
        members = {}
        for rid in ids:
            e = store.get(rid)
            path = (a.paths or {}).get(rid) if hasattr(a, "paths") else None
            if e is None:
                print(f"{rid} is not registered; `synspeq data import {rid} <path>` first", file=sys.stderr)
                return 2
            members[rid] = e["path"]
        flags = {rid: man.resources[rid].flags for rid in ids if rid in man.resources}
        out = build_union(members, a.out, flags)
        store.register("user_union_stock", out["file"], fmt="aizynthfinder_hdf5", terms_acknowledged=True)
        print(f"union of {len(ids)} catalogues: {out['n_unique_inchikeys']:,} unique InChIKeys -> {out['file']}")
        for rid, info in out["members"].items():
            print(f"  {rid:18} records {info['n_records']:>10,}  unparseable {info['n_unparseable']:>7,}  unique {info['n_unique_keys']:>10,}  new {info['new_keys_added']:>10,}")
        if out["noncommercial"] or out["sharealike"]:
            print("  note: the union inherits a non-commercial / share-alike condition from its members")
        print("  registered as user_union_stock; enable it with --resources user_union_stock")
        return 0
    if sub == "download":
        r = man.resources.get(a.id)
        if r is None:
            print(f"unknown resource {a.id!r}", file=sys.stderr)
            return 2
        print(f"{r.name}\n  source: {r.url}\n  terms : {r.terms}\n  cite  : {r.cite}")
        if r.warn:
            print(f"  note  : {r.warn}")
        if r.id in ("aizynth_zinc", "aizynth_models"):
            if not a.accept_terms:
                print("Re-run with --accept-terms to confirm you have read the terms above.", file=sys.stderr)
                return 2
            exe = shutil.which("download_public_data")
            if exe is None or not a.dest:
                print("Needs AiZynthFinder installed (`pip install aizynthfinder`) and --dest <dir>; then it runs `download_public_data <dir>`.", file=sys.stderr)
                return 2
            os.makedirs(a.dest, exist_ok=True)
            rc = subprocess.call([exe, a.dest])
            if rc == 0:
                print(f"Downloaded to {a.dest}. Register the files with `synspeq data import {a.id} <file>` (stock: zinc_stock.hdf5; models: the .onnx/.hdf5 files listed in config.yml).")
            return rc
        if r.fetch:
            if not a.accept_terms or not a.dest:
                print("Re-run with --dest <dir> and --accept-terms (confirms you have read the terms above).", file=sys.stderr)
                return 2
            return _fetch(r, a.dest, store)
        print(f"\nAutomatic download is not implemented for {r.id} in this version. How to get it:\n  {r.get}\nThen: synspeq data import {r.id} <path>")
        return 0
    return 2


# ------------------------------------------------------------------ validate / replay
def _population(a):
    recs = _read_any(a.input, smiles_col=a.smiles_col, name_col=a.name_col)
    return make_population(recs, os.path.basename(a.input), a.pocket_field, first=a.first, last=a.last)


def cmd_validate(a) -> int:
    acc = _population(a).accounting()
    print(json.dumps(acc, indent=1))
    return 0


def cmd_p1(a) -> int:
    pop = _population(a)
    dirs = a.routes_dir
    outcomes = {}
    for m in pop.molecules:
        name = m.source_record.name
        path = next((os.path.join(d, f"{name}.json") for d in dirs if name and os.path.exists(os.path.join(d, f"{name}.json"))), None)
        outcomes[m.mol_id] = load_outcome(path)
    res = compute_p1(pop, outcomes)
    out = dict(synspeq_version=__version__, population=pop.accounting(), p1=res)
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(out, fh, indent=1)
    c = res["counts"]
    print(f"{pop.label}: {res['n_input_records']} records, {res['n_invalid_records']} invalid, {res['n_distinct']} distinct, duplicate rate {_fmt_pct(res['duplicate_rate'])}")
    print(f"  direct {c['direct_purchase']}  routed {c['routed']}  found-no-saved-route {c['found_no_saved_route']}  unresolved {c['unresolved']}  error {c['execution_error']}  pending {c['pending']}")
    print(f"  Find {_fmt_pct(res['find'])}  Find>=1 {_fmt_pct(res['find_ge1'])}  direct share of found {_fmt_pct(res['direct_purchase_share_of_found'])}  Find(raw) {_fmt_pct(res['find_raw'])}  complete={res['complete']}")
    return 0


def cmd_evaluate(a) -> int:
    from .pipeline import evaluate, load_outcomes
    man = load_manifest()
    store = Store(a.data_dir)
    cfg = load_config(a.config, dict(on_missing=a.on_missing, profile=a.profile, pillars=a.pillars.split(",") if a.pillars else None, pocket_field=a.pocket_field,
                                    resources=a.resources.split(",") if a.resources else None))
    res = resolve(cfg, man, store, mode="replay")
    for m in res.messages:
        print(f"[{m['level']}] {m['scope']}: {m['text']}" + (f"\n        fix: {m['fix']}" if m.get("fix") else ""), file=sys.stderr)
    if res.exit_code:
        print("blocked: nothing was computed (see messages above)", file=sys.stderr)
        return res.exit_code
    pf = a.pocket_field or cfg.get("pocket_field")
    recs = _read_any(a.input, smiles_col=a.smiles_col, name_col=a.name_col)
    pop = make_population(recs, os.path.basename(a.input), pf, first=a.first, last=a.last)
    outcomes = load_outcomes(pop, a.routes_dir)
    prior = None
    if a.scaffold_prior:
        with open(a.scaffold_prior) as fh:
            prior = {l.split("\t")[0]: int(l.split("\t")[1]) for l in fh if l.strip()}
    rep = evaluate(pop, outcomes, res, a.out, store.all(), real_lookup_dir=a.real_lookup, scaffold_prior=prior, cache_dir=os.path.join(store.root, "cache"))
    for line in res.banner:
        print("NOTE", line)
    p1 = rep.get("P1")
    if p1:
        print(f"P1  Find {_fmt_pct(p1['find'])}  Find>=1 {_fmt_pct(p1['find_ge1'])}  distinct {p1['n_distinct']}  complete={p1['complete']}")
    p2 = rep.get("P2")
    if p2:
        print(f"P2  supported (primary) {_fmt_pct(p2['primary']['conditional'])} of routed, (any saved) {_fmt_pct(p2['any_saved']['conditional'])}")
    p3 = rep.get("P3")
    if p3 and p3.get("n_routed"):
        print(f"P3  median residual steps {p3['median_residual_steps']}  frontier found {_fmt_pct(p3['frontier_found'])}  needs on-demand {_fmt_pct(p3['needs_on_demand'])}")
    p4 = rep.get("P4")
    if p4:
        b = p4["routed"]
        print(f"P4  N50 {b['n50']}  unique scaffold ratio {b['unique_scaffold_ratio']:.3f}  CSR-AUC {b['csr_auc']:.3f}" if b["n_molecules"] else "P4  no routed molecules")
    print(f"wrote {a.out}/report.html, report.json, molecules.csv, run_manifest.json")
    return 0


def _planner_config(res, man, store, out_dir: str) -> Optional[str]:
    """AiZynthFinder config from the resolved resources: models dir (config.yml) + every enabled planner-stock file."""
    from .planners import build_config
    models = res.states["aizynth_models"]
    if models.state != "ready":
        return None
    mc = os.path.join(models.path, "config.yml") if os.path.isdir(models.path) else models.path
    stocks = {}
    for rid, st in res.states.items():
        if st.state == "ready" and "planner_stock" in man.resources[rid].provides:
            p = st.path
            if os.path.isdir(p):                       # the AiZynthFinder data folder holds zinc_stock.hdf5
                p = os.path.join(p, "zinc_stock.hdf5")
            stocks[rid] = p
    return build_config(mc, stocks, os.path.join(out_dir, "aizynth_config.yml")) if stocks else None


def cmd_plan(a) -> int:
    """Plan every distinct molecule with AiZynthFinder (resumable; one chunk per SLURM array task)."""
    from .planners import plan_population
    man = load_manifest()
    store = Store(a.data_dir)
    cfg = load_config(a.config, dict(on_missing=a.on_missing, pillars=["P1"], resources=a.resources.split(",") if a.resources else None,
                                    pocket_field=a.pocket_field))
    res = resolve(cfg, man, store, mode="plan")
    for m in res.messages:
        print(f"[{m['level']}] {m['scope']}: {m['text']}" + (f"\n        fix: {m['fix']}" if m.get("fix") else ""), file=sys.stderr)
    if res.exit_code:
        print("blocked: the planner needs AiZynthFinder models and a stock (see above)", file=sys.stderr)
        return res.exit_code
    os.makedirs(a.out, exist_ok=True)
    conf = _planner_config(res, man, store, a.out)
    pf = a.pocket_field or cfg.get("pocket_field")
    pop = make_population(_read_any(a.input, smiles_col=a.smiles_col, name_col=a.name_col), os.path.basename(a.input), pf, first=a.first, last=a.last)
    items = [(m.mol_id, m.identity) for m in pop.molecules]
    info = plan_population(items, conf, a.out, workers=a.workers, chunk=a.chunk, n_chunks=a.n_chunks)
    with open(os.path.join(a.out, f"plan_manifest_chunk{a.chunk}.json"), "w") as fh:
        json.dump(dict(synspeq_version=__version__, config=conf, population=pop.accounting(), **info), fh, indent=1)
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
    from .purchasing import scan_membership
    from .purchasing.catalogue import catalogue_files
    files = catalogue_files(a.real_dir)
    files = files[a.shard::a.n_shards]
    queries = {l.strip() for l in open(a.queries) if l.strip()}
    hits = scan_membership(queries, files, smiles_col=a.smiles_col, already_canonical=not a.recanonicalise)
    with open(a.out, "w") as fh:
        for k, f in sorted(hits.items()):
            fh.write(f"{k}\t{f}\n")
    print(f"shard {a.shard}/{a.n_shards}: {len(files)} files scanned, {len(hits)}/{len(queries)} queries found -> {a.out}")
    return 0


def cmd_real_merge(a) -> int:
    from .purchasing import RealLookup
    queries = {l.strip() for l in open(a.queries) if l.strip()}
    hits = set()
    for f in a.hits:
        hits |= {l.split("\t")[0].strip() for l in open(f) if l.strip()}
    lk = RealLookup.load(a.out) if os.path.isdir(a.out) else RealLookup()
    lk.add_scan(queries, hits)
    lk.save(a.out)
    import rdkit
    with open(os.path.join(a.out, "meta.json"), "w") as fh:
        json.dump(dict(rdkit_version=rdkit.__version__, n_queried=len(lk.queried), n_hits=len(lk.hits)), fh, indent=1)
    print(f"lookup now holds {len(lk.queried)} queried molecules, {len(lk.hits)} hits -> {a.out}")
    return 0


def cmd_prior(a) -> int:
    """Scaffold prior from a random draw of rows from random catalogue files (default: 3000 parts x 100 rows)."""
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


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(prog="synspeq", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"synspeq {__version__}")
    ap.add_argument("--data-dir", default=None, help="resource store (default: $SYNSPEQ_DATA_DIR, else ./.synspeq in the current directory)")
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
    x = ds.add_parser("build-union", help="merge registered catalogues into one planner stock")
    x.add_argument("--resources", required=True, help="comma list of registered resource ids")
    x.add_argument("--out", required=True)
    x = ds.add_parser("check")
    x.add_argument("--config"); x.add_argument("--profile", default=None); x.add_argument("--on-missing", default=None, choices=["stop", "skip"])
    x.add_argument("--pillars", default=None, help="comma list, e.g. P1,P2")
    x.add_argument("--resources", default=None, help="comma list of resource ids to enable")
    x.add_argument("--mode", default="plan", choices=["plan", "replay"])
    d.set_defaults(func=cmd_data)

    def common(p):
        p.add_argument("input")
        p.add_argument("--pocket-field", default=None)
        p.add_argument("--smiles-col", default="smiles")
        p.add_argument("--name-col", default="mol_id")
        g = p.add_mutually_exclusive_group()
        g.add_argument("--first", type=int, default=None)
        g.add_argument("--last", type=int, default=None)
    v = sp.add_parser("validate", help="ingestion audit: valid/invalid records, raw and distinct counts")
    common(v); v.set_defaults(func=cmd_validate)
    r = sp.add_parser("p1-replay", help="Pillar 1 from cached AiZynthFinder records")
    common(r); r.add_argument("--routes-dir", action="append", required=True); r.add_argument("--out", default=None)
    r.set_defaults(func=cmd_p1)
    e = sp.add_parser("evaluate", help="all assessable pillars from cached planner records (replay mode)")
    common(e)
    e.add_argument("--routes-dir", action="append", required=True)
    e.add_argument("--out", required=True)
    e.add_argument("--config"); e.add_argument("--profile", default=None); e.add_argument("--on-missing", default=None, choices=["stop", "skip"])
    e.add_argument("--pillars", default=None)
    e.add_argument("--resources", default=None, help="comma list of resource ids to enable, e.g. enamine_bb,enamine_real")
    e.add_argument("--real-lookup", default=None, help="folder with hits.txt / queried.txt from `real-scan` + `real-merge`")
    e.add_argument("--scaffold-prior", default=None, help="TSV scaffold<TAB>count from a reference catalogue sample (P4 reference affinity)")
    e.set_defaults(func=cmd_evaluate)
    pl = sp.add_parser("plan", help="plan with AiZynthFinder (needs the 'planner' extras, models and a stock)")
    common(pl)
    pl.add_argument("--out", required=True)
    pl.add_argument("--config"); pl.add_argument("--on-missing", default=None, choices=["stop", "skip"])
    pl.add_argument("--resources", default=None)
    pl.add_argument("--workers", type=int, default=1)
    pl.add_argument("--chunk", type=int, default=0); pl.add_argument("--n-chunks", type=int, default=1)
    pl.set_defaults(func=cmd_plan)
    rn = sp.add_parser("run", help="plan, then evaluate (SDF to report in one command)")
    common(rn)
    rn.add_argument("--out", required=True)
    rn.add_argument("--config"); rn.add_argument("--profile", default=None); rn.add_argument("--on-missing", default=None, choices=["stop", "skip"])
    rn.add_argument("--pillars", default=None); rn.add_argument("--resources", default=None)
    rn.add_argument("--workers", type=int, default=1)
    rn.add_argument("--real-lookup", default=None); rn.add_argument("--scaffold-prior", default=None)
    rn.set_defaults(func=cmd_run)
    x = sp.add_parser("report", help="re-render report.html from an evaluation folder's report.json")
    x.add_argument("dir"); x.set_defaults(func=cmd_report)
    x = sp.add_parser("real-scan", help="one shard of the REAL membership scan for the queries written by evaluate")
    x.add_argument("--queries", required=True); x.add_argument("--real-dir", required=True)
    x.add_argument("--shard", type=int, default=0); x.add_argument("--n-shards", type=int, default=1)
    x.add_argument("--smiles-col", default="smiles"); x.add_argument("--out", required=True)
    x.add_argument("--recanonicalise", action="store_true",
                   help="canonicalise every catalogue SMILES with this RDKit before matching (slow; use when the catalogue was "
                        "canonicalised by a different RDKit version)")
    x.set_defaults(func=cmd_real_scan)
    x = sp.add_parser("scaffold-prior", help="P4 reference scaffold prior: random draw of rows from random catalogue files (REAL)")
    x.add_argument("--real-dir", required=True); x.add_argument("--out", required=True)
    x.add_argument("--seed", type=int, required=True); x.add_argument("--n-parts", type=int, default=3000)
    x.add_argument("--rows-per-part", type=int, default=100); x.add_argument("--smiles-col", default="smiles")
    x.set_defaults(func=cmd_prior)
    x = sp.add_parser("real-merge", help="merge shard hit files into a REAL lookup folder")
    x.add_argument("--queries", required=True); x.add_argument("--hits", nargs="+", required=True); x.add_argument("--out", required=True)
    x.set_defaults(func=cmd_real_merge)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
