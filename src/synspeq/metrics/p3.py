"""Pillar 3: deliverability on the primary route of routed molecules (purchasing frontier, strict matching)."""
from __future__ import annotations
from statistics import median as _median


def median(values):
    m = _median(list(values))
    return int(m) if float(m).is_integer() else m
from typing import Dict, List, Optional
from ..populations import Population
from ..purchasing import Catalogue, frontier_for_outcome


def compute_p3(pop: Population, outcomes: Dict[str, object], cat: Catalogue) -> dict:
    rows: List[dict] = []
    for m in pop.molecules:
        o = outcomes.get(m.mol_id)
        if o is None or o.outcome != "routed":
            continue
        f = frontier_for_outcome(o, cat)
        if f is None:
            continue
        rows.append(dict(mol_id=m.mol_id, pocket=m.pocket, full_steps=f.full_steps, residual_steps=f.residual_steps,
                         residual_depth=f.residual_depth, tier=f.tier, n_unknown=f.n_unknown, n_frontier=f.n_frontier))
    n = len(rows)
    if n == 0:
        return dict(label=pop.label, n_routed=0, rows=rows)
    n_unknown_routes = sum(1 for r in rows if r["n_unknown"] > 0)
    out = dict(
        label=pop.label, n_routed=n, rows=rows,
        median_route_steps=median(r["full_steps"] for r in rows),
        median_residual_steps=median(r["residual_steps"] for r in rows),
        frontier_found=sum(1 for r in rows if r["tier"] != "none") / n,
        needs_on_demand=sum(1 for r in rows if r["tier"] == "on_demand") / n,
        n_routes_with_unrecorded_real_state=n_unknown_routes,
        residual_steps_counts={str(k): sum(1 for r in rows if r["residual_steps"] == k) for k in sorted({r["residual_steps"] for r in rows})},
        tiers_available=dict(stocked=cat.has_stocked, on_demand=cat.has_on_demand),
        bound=("residual steps are upper bounds and the on-demand share is a lower bound for "
               f"{n_unknown_routes} routes with intermediates never looked up in REAL") if n_unknown_routes else None,
        denominators=dict(all="routed molecules (primary route)"),
    )
    return out
