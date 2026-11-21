"""Helpers shared by the command modules."""
from __future__ import annotations
import argparse
from .config import load_config
from .resources import Store


def pos_int(v):
    try:
        n = int(v)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{v!r} is not an integer")
    if n < 1:
        raise argparse.ArgumentTypeError(f"must be >= 1, got {n}")
    return n


def nonneg_int(v):
    n = int(v)
    if n < 0:
        raise argparse.ArgumentTypeError(f"must be >= 0, got {n}")
    return n


def config_and_store(a, extra=None):
    """Configuration first, then the store: command line > config file > $SYNSPEQ_DATA_DIR > ./.synspeq."""
    ov = dict(on_missing=getattr(a, "on_missing", None), profile=getattr(a, "profile", None),
              pillars=a.pillars.split(",") if getattr(a, "pillars", None) else None,
              resources=a.resources.split(",") if getattr(a, "resources", None) else None,
              planner_stock=a.planner_stock.split(",") if getattr(a, "planner_stock", None) else None,
              pocket_field=getattr(a, "pocket_field", None))
    ov.update(extra or {})
    cfg = load_config(getattr(a, "config", None), ov)
    return cfg, Store(a.data_dir or cfg.get("data_dir"))


def add_measures(parser):
    parser.add_argument("--measures", dest="pillars", default=None, metavar="MEASURES",
                        help="comma list of findability, precedent-support, residual-steps, scaffold-diversity (default: all)")
    parser.add_argument("--pillars", dest="pillars", default=None, help=argparse.SUPPRESS)    # older name of --measures
