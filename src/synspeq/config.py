"""Run configuration (config.yaml). Unknown keys are rejected so a typo cannot silently change a run."""
from __future__ import annotations
import yaml

ALLOWED = {"profile", "resources", "paths", "on_missing", "auto_download", "pillars", "data_dir", "pocket_field"}
PILLARS = ("P1", "P2", "P3", "P4")


def load_config(path: str | None = None, overrides: dict | None = None) -> dict:
    cfg = {}
    if path:
        with open(path) as fh:
            cfg = yaml.safe_load(fh) or {}
    cfg.update({k: v for k, v in (overrides or {}).items() if v is not None})
    bad = set(cfg) - ALLOWED
    if bad:
        raise ValueError(f"unknown config keys: {sorted(bad)} (allowed: {sorted(ALLOWED)})")
    cfg.setdefault("profile", "default")
    cfg.setdefault("resources", [])
    cfg.setdefault("paths", {})
    cfg.setdefault("on_missing", "stop")
    cfg.setdefault("auto_download", False)
    cfg.setdefault("pillars", list(PILLARS))
    if cfg["on_missing"] not in ("stop", "skip"):
        raise ValueError("on_missing must be 'stop' or 'skip'")
    if cfg["profile"] not in ("default", "custom"):
        raise ValueError("profile must be default | custom")
    for p in cfg["pillars"]:
        if p not in PILLARS:
            raise ValueError(f"unknown pillar {p!r}")
    return cfg
