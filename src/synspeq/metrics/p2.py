"""Pillar 2: chemical support = retrieved reaction-precedent support of the saved routes.

Statement the report must carry: retrieved precedent evidence is not experimental feasibility; no qualifying
precedent means support was not retrieved (not that the reaction cannot work); 0.6412 is a similarity cutoff, not a
probability; the reaction database is shared with the planner's training data.
"""
from __future__ import annotations
from typing import Callable, Dict, List, Optional
from ..populations import Population
from ..support.aggregate import THRESHOLD, molecule_state, route_forward_reactions, route_state, step_state


def p2_molecules(pop: Population, outcomes: Dict[str, object], features: Dict[str, dict], threshold: float = THRESHOLD) -> List[dict]:
    """One record per distinct molecule that has saved routes. features: forward reaction string -> feature row."""
    out = []
    for m in pop.molecules:
        o = outcomes.get(m.mol_id)
        if o is None or o.outcome != "routed":
            continue
        rec = dict(mol_id=m.mol_id, pocket=m.pocket, primary_len=len(o.routes[0]), n_saved_routes=len(o.routes))
        for view in ("all", "analogue"):
            rs = [route_state([step_state(features.get(fr), view, threshold) for fr in route_forward_reactions(rt)]) for rt in o.routes]
            rec[f"primary_{view}"] = rs[0]
            rec[f"any_{view}"] = molecule_state(rs)
        rec["any_via_alternative"] = int(rec["any_all"] == 1 and rec["primary_all"] != 1)
        rec["n_unscored_steps_primary"] = sum(1 for fr in route_forward_reactions(o.routes[0]) if step_state(features.get(fr), "all", threshold) == -1)
        out.append(rec)
    return out


def required_reactions(pop: Population, outcomes: Dict[str, object]) -> List[str]:
    seen, out = set(), []
    for m in pop.molecules:
        o = outcomes.get(m.mol_id)
        if o is None or o.outcome != "routed":
            continue
        for rt in o.routes:
            for fr in route_forward_reactions(rt):
                if fr not in seen:
                    seen.add(fr)
                    out.append(fr)
    return out


def compute_p2(pop: Population, outcomes: Dict[str, object], features: Dict[str, dict], n_distinct: Optional[int] = None,
               threshold: float = THRESHOLD) -> dict:
    mols = p2_molecules(pop, outcomes, features, threshold)
    n_routed = len(mols)
    n_distinct = n_distinct if n_distinct is not None else len(pop.molecules)
    n_direct = sum(1 for m in pop.molecules if outcomes.get(m.mol_id) is not None and outcomes[m.mol_id].outcome == "direct_purchase")

    def block(scope: str, view: str):
        k = f"{scope}_{view}"
        sup = sum(1 for r in mols if r[k] == 1)
        ind = sum(1 for r in mols if r[k] == -1)
        return dict(supported=sup, indeterminate=ind, unsupported=n_routed - sup - ind,
                    conditional=(sup / n_routed) if n_routed else None, unconditional=(sup / n_distinct) if n_distinct else None)
    by_len = {}
    for b in (1, 2, 3, 4):
        sel = [r for r in mols if min(r["primary_len"], 4) == b]
        by_len[("4+" if b == 4 else str(b))] = dict(n=len(sel), primary=(sum(1 for r in sel if r["primary_all"] == 1) / len(sel)) if sel else None,
                                                   any_saved=(sum(1 for r in sel if r["any_all"] == 1) / len(sel)) if sel else None)
    return dict(label=pop.label, threshold=threshold, n_routed=n_routed, n_distinct=n_distinct, n_direct_purchase=n_direct,
                primary=block("primary", "all"), any_saved=block("any", "all"),
                primary_analogue=block("primary", "analogue"), any_saved_analogue=block("any", "analogue"),
                by_primary_route_length=by_len, molecules=mols,
                notes=["retrieved precedent evidence is not experimental feasibility",
                       "no qualifying precedent = support not retrieved",
                       f"{threshold} is a similarity cutoff, not a probability",
                       "the USPTO database is shared with the planner's training data"])
