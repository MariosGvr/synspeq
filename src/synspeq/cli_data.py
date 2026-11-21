"""`synspeq data ...`: list, info, import, check, download, build-union."""
from __future__ import annotations
import json
import os
import shutil
import subprocess
import sys
from .cli_common import config_and_store
from .config import MEASURE_NAMES
from .environment import rdkit_note
from .resources import Store, load_manifest, resolve


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
    store.register(r.id, path, fmt=r.format, terms_acknowledged=True, check=r.check or None)
    print(f"downloaded and registered {r.id} -> {path} (md5 verified)")
    return 0


# ------------------------------------------------------------------ data
def cmd_data(a) -> int:
    man = load_manifest()
    if a.data_cmd == "check":
        cfg, store = config_and_store(a)
    else:
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
        e = store.register(a.id, a.path, fmt=r.format, hash_it=not a.no_hash, snapshot_date=a.snapshot_date, terms_acknowledged=a.accept_terms,
                           check=r.check or None)
        v = e["validation"]
        print(f"registered {a.id} -> {e['path']}  validation: {v.get('status')}" + (f" ({v.get('reason')})" if v.get("reason") else ""))
        return 0 if v.get("status") in ("ok", "unchecked") else 2
    if sub == "status":
        print(json.dumps(store.all(), indent=1))
        return 0
    if sub == "check":
        res = resolve(cfg, man, store, mode=a.mode)
        print(f"resource store: {store.root}")
        note = rdkit_note()
        if note:
            print(f"NOTE {note}")
        print("resources")
        for s in res.states.values():
            if s.state != "not_requested":
                print(f"  {s.id:22} {s.state:10} requested by {s.requested_by:8} {s.path or ''} {('- ' + s.reason) if s.reason else ''}")
        print("metrics")
        for p, m in res.metrics.items():
            print(f"  {MEASURE_NAMES.get(p, p)}: {m.status}" + (f"  ({m.reason})" if m.reason else ""))
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
            if e is None:
                print(f"{rid} is not registered; `synspeq data import {rid} <path>` first", file=sys.stderr)
                return 2
            members[rid] = e["path"]
        flags = {rid: man.resources[rid].flags for rid in ids if rid in man.resources}
        out = build_union(members, a.out, flags)
        store.register("user_union_stock", out["file"], fmt=None, terms_acknowledged=True)
        print(f"union of {len(ids)} catalogs: {out['n_unique_inchikeys']:,} unique InChIKeys -> {out['file']}")
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
                # download_public_data fetches models, templates and the stock together: register both resources
                dest = os.path.abspath(a.dest)
                stock = os.path.join(dest, "zinc_stock.hdf5")
                if os.path.exists(os.path.join(dest, "config.yml")):
                    e = store.register("aizynth_models", dest, fmt=man.resources["aizynth_models"].format, terms_acknowledged=True,
                                       check=man.resources["aizynth_models"].check or None)
                    print(f"registered aizynth_models -> {dest}  validation: {e['validation'].get('status')}")
                if os.path.exists(stock):
                    e = store.register("aizynth_zinc", stock, fmt=man.resources["aizynth_zinc"].format, terms_acknowledged=True)
                    print(f"registered aizynth_zinc -> {stock}  validation: {e['validation'].get('status')}")
                else:
                    print(f"{stock} not found; register the stock with `synspeq data import aizynth_zinc <file>`", file=sys.stderr)
            return rc
        if r.fetch:
            if not a.accept_terms or not a.dest:
                print("Re-run with --dest <dir> and --accept-terms (confirms you have read the terms above).", file=sys.stderr)
                return 2
            return _fetch(r, a.dest, store)
        print(f"\nAutomatic download is not implemented for {r.id} in this version. How to get it:\n  {r.get}\nThen: synspeq data import {r.id} <path>")
        return 0
    return 2


