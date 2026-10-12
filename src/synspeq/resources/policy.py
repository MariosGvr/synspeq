"""Missing-resource policy (resource_manifest.yaml: missing_resource_policy, metric_requirements, use_cases).

resolve(config, manifest, store, mode) -> Resolution with
  * a state per resource: ready | partial | not_requested | missing | invalid
  * a status per pillar:  value | bound | not_assessed | unavailable | blocked   (+ sub-metrics for P3, P4)
  * exit_code 2 when anything is blocked, else 0, and the messages/banner a report must print.
Rules implemented: R1 unspecified is not an error; R2 requested-but-absent is an error with the exact fix (skip
turns it into 'unavailable'); R3 a user-named stock is never replaced by the default; R4 absence is never zero;
R5 partial resources give bounds; core planner resources are never skippable.
"""
from __future__ import annotations
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional
from .manifest import Manifest
from .store import Store, validate

STOCK_CAP, MODELS_CAP, DB_CAP = "planner_stock", "planner_models", "reaction_database"
STOCKED, ONDEMAND, SCAFFOLD_REF = "purchasing_stocked_tier", "purchasing_on_demand_tier", "reference_scaffold_profile"


@dataclass
class ResourceState:
    id: str
    state: str                         # ready | partial | not_requested | missing | invalid
    requested_by: Optional[str] = None  # user | profile | default | core
    path: Optional[str] = None
    reason: Optional[str] = None
    coverage: Optional[dict] = None
    fix: Optional[str] = None


@dataclass
class MetricStatus:
    pillar: str
    status: str                        # value | bound | not_assessed | unavailable | blocked
    reason: str = ""
    sub: Dict[str, dict] = field(default_factory=dict)


@dataclass
class Resolution:
    states: Dict[str, ResourceState]
    metrics: Dict[str, MetricStatus]
    messages: List[dict]
    exit_code: int

    @property
    def banner(self) -> List[str]:
        out = []
        for p, m in self.metrics.items():
            if m.status in ("not_assessed", "unavailable", "blocked", "bound"):
                out.append(f"{p}: {m.status}" + (f" ({m.reason})" if m.reason else ""))
            for s, d in m.sub.items():
                if d["status"] != "value":
                    out.append(f"{p}.{s}: {d['status']}" + (f" ({d['reason']})" if d.get("reason") else ""))
        return out


def _fix(manifest: Manifest, rid: str) -> str:
    r = manifest.resources[rid]
    how = f"synspeq data download {rid}" if r.automated else f"synspeq data import {rid} <path>"
    return f"{how}   ({r.get})"


def _locate(rid: str, cfg: dict, store: Store, fmt: Optional[str]):
    """-> (path | None, coverage | None, validation | None, source)"""
    p = (cfg.get("paths") or {}).get(rid)
    if p:
        p = os.path.abspath(os.path.expanduser(p))
        return p, None, validate(p, fmt), "config"
    e = store.get(rid)
    if e:
        v = validate(e["path"], fmt) if not os.path.exists(e["path"]) else e.get("validation") or validate(e["path"], fmt)
        return e["path"], e.get("coverage"), v, "registry"
    return None, None, None, None


def resolve(cfg: dict, manifest: Manifest, store: Store, mode: str = "plan") -> Resolution:
    """mode: 'plan' (planner runs: stock + models needed for P1) or 'replay' (cached planner records supplied)."""
    pillars = set(cfg["pillars"])
    skip = cfg.get("on_missing") == "skip"
    R = manifest.resources
    unknown = [r for r in cfg.get("resources", []) if r not in R] + [r for r in (cfg.get("paths") or {}) if r not in R]
    if unknown:
        raise ValueError(f"unknown resource ids in config: {sorted(set(unknown))} (see `synspeq data list`)")

    requested: Dict[str, str] = {}
    for r in cfg.get("resources", []):
        requested[r] = "user"
    for r in (cfg.get("paths") or {}):
        requested.setdefault(r, "user")
    # R1/R3: a planner stock is added by default only when the user named none.
    if mode == "plan" and "P1" in pillars and not any(STOCK_CAP in R[r].provides for r in requested):
        requested.setdefault("aizynth_zinc", "default")
    if mode == "plan" and "P1" in pillars:
        requested.setdefault("aizynth_models", "core")
    if "P2" in pillars:
        requested.setdefault("uspto_rxninsight_db", "core")

    states: Dict[str, ResourceState] = {}
    for rid, res in R.items():
        if rid not in requested:
            states[rid] = ResourceState(rid, "not_requested")
            continue
        path, cov, v, src = _locate(rid, cfg, store, res.format)
        fix = _fix(manifest, rid)
        if path is None:
            states[rid] = ResourceState(rid, "missing", requested[rid], None, "not configured or imported", fix=fix)
        elif v is None or v.get("status") == "invalid":
            why = (v or {}).get("reason", "invalid")
            states[rid] = ResourceState(rid, "missing" if why in ("path_not_found",) else "invalid", requested[rid], path, why, fix=fix)
        elif cov and cov.get("fraction", 1.0) < 1.0:
            states[rid] = ResourceState(rid, "partial", requested[rid], path, f"coverage {cov['fraction']:.3f}", cov)
        else:
            states[rid] = ResourceState(rid, "ready", requested[rid], path)

    have = set()
    partial_caps = set()
    for s in states.values():
        if s.state == "ready":
            have.update(R[s.id].provides)
        elif s.state == "partial":
            have.update(R[s.id].provides)
            partial_caps.update(R[s.id].provides)

    def broken(cap: str) -> List[ResourceState]:
        return [s for s in states.values() if s.state in ("missing", "invalid") and cap in R[s.id].provides]

    msgs: List[dict] = []
    metrics: Dict[str, MetricStatus] = {}

    def note(level, scope, text, fix=None):
        msgs.append(dict(level=level, scope=scope, text=text, fix=fix))

    # ---- P1
    if "P1" not in pillars:
        metrics["P1"] = MetricStatus("P1", "not_assessed", "pillar not selected")
    elif mode == "replay":
        metrics["P1"] = MetricStatus("P1", "value", "cached planner records; no stock or models needed")
    else:
        miss = []
        if MODELS_CAP not in have:
            miss.append("planner models")
            for s in broken(MODELS_CAP) or [states["aizynth_models"]]:
                note("error", "P1", f"{s.id}: {s.reason}", s.fix)
        if STOCK_CAP not in have:
            miss.append("planner stock")
            for s in broken(STOCK_CAP):
                note("error", "P1", f"{s.id}: {s.reason} (requested by {s.requested_by}; the tool will not fall back to another stock)", s.fix)
        metrics["P1"] = MetricStatus("P1", "blocked", "missing: " + ", ".join(miss)) if miss else MetricStatus("P1", "value")

    # ---- P2
    if "P2" not in pillars:
        metrics["P2"] = MetricStatus("P2", "not_assessed", "pillar not selected")
    elif DB_CAP in have:
        metrics["P2"] = MetricStatus("P2", "value")
    else:
        for s in broken(DB_CAP) or [states["uspto_rxninsight_db"]]:
            note("error", "P2", f"{s.id}: {s.reason}", s.fix)
        metrics["P2"] = MetricStatus("P2", "unavailable" if skip else "blocked", "reaction database missing")

    # ---- P3
    if "P3" not in pillars:
        metrics["P3"] = MetricStatus("P3", "not_assessed", "pillar not selected")
    else:
        bs, od = STOCKED in have, ONDEMAND in have
        asked = [s for s in states.values() if s.requested_by and {STOCKED, ONDEMAND} & set(R[s.id].provides) and s.state in ("missing", "invalid")]
        asked_purchasing = [s for s in asked if s.requested_by in ("user", "profile")]
        sub = {}
        if bs or od:
            which = ("stocked catalogue" if bs else "") + (" + " if bs and od else "") + ("on-demand catalogue" if od else "")
            top = "bound" if ({STOCKED, ONDEMAND} & partial_caps) else "value"
            sub["residual_steps_and_frontier"] = dict(status=top, reason=f"frontier defined by: {which}" + ("; partial catalogue coverage -> bounds" if top == "bound" else ""))
            if bs and od:
                sub["real_needed"] = dict(status="bound" if ONDEMAND in partial_caps else "value", reason="partial on-demand coverage: lower bound" if ONDEMAND in partial_caps else "")
            else:
                sub["real_needed"] = dict(status="not_assessed", reason="needs both the stocked and the on-demand tier")
            for s in asked_purchasing:
                note("warning", "P3", f"{s.id}: {s.reason}; the affected tier is {'unavailable' if skip else 'blocked'}", s.fix)
                if not skip:
                    top = "blocked"
            metrics["P3"] = MetricStatus("P3", top, sub["residual_steps_and_frontier"]["reason"], sub)
            if top == "blocked":
                metrics["P3"].reason = "a requested catalogue is missing"
        elif asked_purchasing:
            for s in asked_purchasing:
                note("error", "P3", f"{s.id}: {s.reason}", s.fix)
            metrics["P3"] = MetricStatus("P3", "unavailable" if skip else "blocked", "requested purchasing catalogue missing")
        else:
            metrics["P3"] = MetricStatus("P3", "not_assessed", "no purchasing catalogue configured (REAL / Enamine building blocks); planner-stock hits are already leaves, so the default stock alone adds no frontier information")

    # ---- P4
    if "P4" not in pillars:
        metrics["P4"] = MetricStatus("P4", "not_assessed", "pillar not selected")
    else:
        sub = {}
        if SCAFFOLD_REF in have:
            sub["enamine_affinity"] = dict(status="bound" if SCAFFOLD_REF in partial_caps else "value", reason="partial reference catalogue" if SCAFFOLD_REF in partial_caps else "")
        else:
            asked = [s for s in states.values() if s.requested_by in ("user", "profile") and s.state in ("missing", "invalid") and SCAFFOLD_REF in R[s.id].provides]
            for s in asked:
                note("warning", "P4", f"{s.id}: {s.reason}", s.fix)
            sub["enamine_affinity"] = dict(status="unavailable" if (asked and skip) else ("blocked" if asked else "not_assessed"),
                                           reason="requested catalogue missing" if asked else "no reference catalogue configured")
        top = "blocked" if sub["enamine_affinity"]["status"] == "blocked" else "value"
        metrics["P4"] = MetricStatus("P4", top, "", sub)

    # hint (never silent): a resource that is registered but not enabled is not used.
    for rid, st in states.items():
        if st.state == "not_requested" and store.get(rid) is not None and rid not in ("aizynth_zinc", "aizynth_models", "uspto_rxninsight_db"):
            note("info", rid, "registered but not enabled, so it is not used; add it to `resources:` in config.yaml or pass --resources", None)
    # R2: any resource the user or profile asked for must be present, whatever capability it provides.
    unmet = [s for s in states.values() if s.state in ("missing", "invalid") and s.requested_by in ("user", "profile")]
    hard = False
    for s in unmet:
        if skip:
            note("warning", s.id, f"{s.reason}; continuing without it (on_missing: skip). Results exclude it and are not comparable with a run that has it.", s.fix)
            if STOCK_CAP in R[s.id].provides and metrics["P1"].status == "value":
                metrics["P1"].reason = (metrics["P1"].reason + "; " if metrics["P1"].reason else "") + f"declared stock incomplete: {s.id} missing"
        else:
            hard = True
            if not any(m["scope"] == s.id for m in msgs):
                note("error", s.id, f"{s.reason}; the run stops before any computation (on_missing: stop)", s.fix)
    code = 2 if hard or any(m.status == "blocked" or any(d["status"] == "blocked" for d in m.sub.values()) for m in metrics.values()) else 0
    return Resolution(states, metrics, msgs, code)
