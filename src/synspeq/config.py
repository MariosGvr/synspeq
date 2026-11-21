"""Run configuration (config.yaml). Unknown keys are rejected so a typo cannot silently change a run."""
from __future__ import annotations
import yaml

ALLOWED = {"profile", "resources", "planner_stock", "paths", "on_missing", "measures", "pillars", "data_dir", "pocket_field"}
PILLARS = ("P1", "P2", "P3", "P4")      # internal identifiers of the four measures; the output uses MEASURE_NAMES and MEASURE_KEYS
MEASURE_NAMES = {"P1": "findability", "P2": "precedent support", "P3": "residual steps", "P4": "scaffold diversity"}
MEASURE_KEYS = {"P1": "findability", "P2": "precedent_support", "P3": "residual_steps", "P4": "scaffold_diversity"}   # keys of report.json
_BY_NAME = {n.replace(" ", "-"): k for k, n in MEASURE_NAMES.items()}


def measure_id(x: str) -> str:
    """A measure given by name (findability, precedent-support, residual-steps, scaffold-diversity) or by its identifier (P1-P4)."""
    key = x.strip().lower().replace("_", "-").replace(" ", "-")
    if key.upper() in PILLARS:
        return key.upper()
    if key in _BY_NAME:
        return _BY_NAME[key]
    raise ValueError(f"unknown measure {x!r} (use findability, precedent-support, residual-steps or scaffold-diversity)")


def _str_list(cfg, key):
    v = cfg.get(key)
    if v is None:
        return []
    if not isinstance(v, list) or not all(isinstance(x, str) for x in v):
        raise ValueError(f"{key} must be a list of strings, got {v!r}")
    return v


def load_config(path: str | None = None, overrides: dict | None = None) -> dict:
    """Precedence: command-line overrides > config file > defaults. `data_dir` follows the same rule
    (the command line's --data-dir is applied by the caller, then this file's data_dir, then $SYNSPEQ_DATA_DIR)."""
    cfg = {}
    if path:
        with open(path) as fh:
            cfg = yaml.safe_load(fh) or {}
        if not isinstance(cfg, dict):
            raise ValueError(f"{path}: the configuration must be a mapping")
    cfg.update({k: v for k, v in (overrides or {}).items() if v is not None})
    bad = set(cfg) - ALLOWED
    if bad:
        raise ValueError(f"unknown config keys: {sorted(bad)} (allowed: {sorted(ALLOWED)})")
    cfg["resources"] = _str_list(cfg, "resources")
    cfg["planner_stock"] = _str_list(cfg, "planner_stock")
    cfg.setdefault("profile", "default")
    cfg.setdefault("paths", {})
    if not isinstance(cfg["paths"], dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in cfg["paths"].items()):
        raise ValueError("paths must map resource ids to path strings")
    cfg.setdefault("on_missing", "stop")
    if "measures" in cfg:               # "measures" is the documented key; "pillars" is the older name for the same setting
        if "pillars" in cfg:
            raise ValueError("give either measures or pillars, not both")
        cfg["pillars"] = cfg.pop("measures")
    cfg.setdefault("pillars", list(PILLARS))
    cfg["pillars"] = [measure_id(x) for x in _str_list(cfg, "pillars")]
    if cfg["on_missing"] not in ("stop", "skip"):
        raise ValueError("on_missing must be 'stop' or 'skip'")
    if cfg["profile"] not in ("default", "custom"):
        raise ValueError("profile must be default | custom")
    return cfg
