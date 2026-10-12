"""Live planning with AiZynthFinder (needs `aizynthfinder`; extras 'planner').

Search settings are fixed: USPTO expansion policy plus ring-breaker (both selected), filter policy on,
maximum depth 10 transforms, 500 iterations, 60 s per molecule, return_first False, target not excluded from the stock,
default cutoffs, one search per molecule. Output per molecule is the same JSON record the replay reads:
{mol_id, smiles, solvable, n_steps, trees, error}; trees are AiZynthFinder's route dictionaries (partial trees included;
a returned route is a tree with every leaf in stock, see routes/aizynth.py).

The stock is whatever the user enabled: AiZynthFinder's own ZINC stock by default, or any set of stock files
(HDF5 or plain SMILES converted by `stock_from_smiles`). Several stocks are used as one union (select_all).
"""
from __future__ import annotations
import json
import multiprocessing as mp
import os
import time
from typing import Dict, Iterable, List, Optional

PLANNER_SETTINGS = dict(max_transforms=10, iteration_limit=500, time_limit=60, return_first=False,
                        exclude_target_from_stock=False, seed=42)


def count_reactions(node) -> int:
    if not isinstance(node, dict):
        return 0
    n = 1 if node.get("type") == "reaction" else 0
    for c in node.get("children", []):
        n += count_reactions(c)
    return n


def build_config(models_config: str, stock_files: Dict[str, str], out_path: str) -> str:
    """Write an AiZynthFinder config: expansion/filter policies from the models config.yml (as produced by
    `download_public_data`), stocks replaced by the given {name: file} map."""
    import yaml
    with open(models_config) as fh:
        cfg = yaml.safe_load(fh)
    base = os.path.dirname(os.path.abspath(models_config))

    def absolutise(v):
        if isinstance(v, str):
            return v if os.path.isabs(v) else os.path.join(base, v)
        if isinstance(v, list):
            return [absolutise(x) for x in v]
        if isinstance(v, dict):
            return {k: absolutise(x) for k, x in v.items()}
        return v
    for key in ("expansion", "filter"):
        if key in cfg:
            cfg[key] = absolutise(cfg[key])
    cfg["stock"] = {k: os.path.abspath(v) for k, v in stock_files.items()}
    # search settings are applied through PLANNER_SETTINGS at run time; keep the file free of them
    cfg.pop("search", None)
    with open(out_path, "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    return out_path


def stock_from_smiles(smiles_files: List[str], out_hdf5: str, source: str = "plain") -> str:
    """Convert plain SMILES catalogues (one SMILES per line) to an AiZynthFinder HDF5 stock (the `smiles2stock` tool)."""
    from aizynthfinder.tools import make_stock
    make_stock.main(["--files", *smiles_files, "--output", out_hdf5, "--source", source])
    return out_hdf5


class AiZynthPlanner:
    def __init__(self, config_path: str, settings: Optional[dict] = None):
        from aizynthfinder.aizynthfinder import AiZynthFinder
        self.settings = dict(PLANNER_SETTINGS, **(settings or {}))
        finder = None
        for attempt in range(5):                      # the large HDF5 stock occasionally fails under read contention
            try:
                finder = AiZynthFinder(configfile=config_path)
                break
            except Exception:
                time.sleep(10 * (attempt + 1))
        if finder is None:
            finder = AiZynthFinder(configfile=config_path)
        finder.expansion_policy.select_all()
        finder.filter_policy.select_all()
        finder.stock.select_all()
        s = finder.config.search
        s.max_transforms = self.settings["max_transforms"]
        s.iteration_limit = self.settings["iteration_limit"]
        s.time_limit = self.settings["time_limit"]
        s.return_first = self.settings["return_first"]
        if hasattr(s, "exclude_target_from_stock"):
            s.exclude_target_from_stock = self.settings["exclude_target_from_stock"]
        for attr in ("random_seed", "seed"):
            if hasattr(s, attr):
                setattr(s, attr, self.settings["seed"])
        self.finder = finder
        self.version = _version()

    def plan(self, mol_id: str, smiles: str) -> dict:
        r = {"mol_id": mol_id, "smiles": smiles, "solvable": False, "n_steps": None, "trees": [], "error": None}
        try:
            self.finder.target_smiles = smiles
            self.finder.tree_search()
            self.finder.build_routes()
            stats = self.finder.extract_statistics()
            r["solvable"] = bool(stats.get("is_solved", False))
            r["trees"] = self.finder.routes.dicts
            if r["solvable"] and r["trees"]:
                r["n_steps"] = count_reactions(r["trees"][0])
        except Exception as exc:                     # recorded, never swallowed: the outcome is `execution_error`
            r["error"] = f"{type(exc).__name__}: {exc}"
        return r


def _version() -> Optional[str]:
    try:
        import aizynthfinder
        return getattr(aizynthfinder, "__version__", None)
    except Exception:
        return None


def _worker(args):
    config, settings, items, out_dir, skip = args
    pl = AiZynthPlanner(config, settings)
    n_ok = n_err = n_skip = 0
    for mol_id, smi in items:
        path = os.path.join(out_dir, f"{mol_id}.json")
        if skip and os.path.exists(path):
            try:
                old = json.load(open(path))
                from ..routes.aizynth import validate_record
                # reuse only a complete record that is for this very molecule
                if old.get("error") is None and validate_record(old) is None and old.get("mol_id") == mol_id and old.get("smiles") == smi:
                    n_skip += 1
                    continue
            except Exception:
                pass
        res = pl.plan(mol_id, smi)
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(res, fh)
        os.replace(tmp, path)                        # atomic: an interrupted run never leaves a half-written record
        n_ok += 1
        n_err += bool(res["error"])
    return n_ok, n_err, n_skip


def plan_population(items: List[tuple], config_path: str, out_dir: str, settings: Optional[dict] = None, workers: int = 1,
                    chunk: int = 0, n_chunks: int = 1, skip_existing: bool = True) -> dict:
    """items: [(mol_id, smiles)]. Plans chunk `chunk` of `n_chunks` (for SLURM arrays) with `workers` processes, resumable."""
    for name, v, lo in (("workers", workers, 1), ("n_chunks", n_chunks, 1)):
        if isinstance(v, bool) or not isinstance(v, int) or v < lo:
            raise ValueError(f"{name} must be an integer >= {lo}, got {v!r}")
    if isinstance(chunk, bool) or not isinstance(chunk, int) or not (0 <= chunk < n_chunks):
        raise ValueError(f"chunk must satisfy 0 <= chunk < n_chunks, got chunk={chunk!r}, n_chunks={n_chunks!r}")
    os.makedirs(out_dir, exist_ok=True)
    mine = items[chunk::n_chunks]
    parts = [mine[i::workers] for i in range(workers)]
    jobs = [(config_path, settings, p, out_dir, skip_existing) for p in parts if p]
    t0 = time.time()
    if len(jobs) <= 1:
        res = [_worker(jobs[0])] if jobs else []
    else:
        with mp.get_context("spawn").Pool(len(jobs)) as pool:
            res = pool.map(_worker, jobs)
    return dict(n_molecules=len(mine), planned=sum(r[0] for r in res), errors=sum(r[1] for r in res), skipped=sum(r[2] for r in res),
                seconds=round(time.time() - t0, 1), settings=dict(PLANNER_SETTINGS, **(settings or {})), chunk=chunk, n_chunks=n_chunks)
