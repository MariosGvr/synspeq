"""Pillar 1: findability (route recovery). Exact counts and denominators; nothing is rounded here.

All shares are over distinct molecules (the evaluation unit) unless stated. A population is `complete` only when
no molecule is pending and none is an execution error; a partial population still reports its counts but its
three-way partition does not sum to 1 and `complete` is False.
"""
from __future__ import annotations
from typing import Dict
from ..populations import Population
from ..routes.aizynth import OUTCOMES, RouteOutcome


def _share(n, d):
    return None if not d else n / d


def compute_p1(pop: Population, outcomes: Dict[str, RouteOutcome]) -> dict:
    """outcomes: mol_id -> RouteOutcome (missing mol_id counts as pending)."""
    cnt = {o: 0 for o in OUTCOMES}
    raw = {o: 0 for o in OUTCOMES}
    for m in pop.molecules:
        o = outcomes.get(m.mol_id)
        k = o.outcome if o is not None else "pending"
        cnt[k] += 1
        raw[k] += m.n_records
    n = len(pop.molecules)
    n_raw = sum(m.n_records for m in pop.molecules)
    found = cnt["direct_purchase"] + cnt["routed"] + cnt["found_no_saved_route"]
    found_ge1 = cnt["routed"] + cnt["found_no_saved_route"]
    found_raw = raw["direct_purchase"] + raw["routed"] + raw["found_no_saved_route"]
    acc = pop.accounting()
    evaluated = n - cnt["pending"] - cnt["execution_error"]
    return dict(
        label=pop.label,
        complete=(cnt["pending"] == 0 and cnt["execution_error"] == 0),
        n_distinct=n, n_raw_valid=n_raw, n_input_records=acc["n_input_records"], n_invalid_records=acc["n_invalid_records"],
        counts=cnt, raw_counts=raw,
        duplicate_rate=acc["duplicate_rate"],
        find=_share(found, n),                         # found / distinct
        find_raw=_share(found_raw, n_raw),             # found over raw valid rows (duplicates inherit their molecule's outcome)
        direct_purchase_share_of_found=_share(cnt["direct_purchase"], found),
        find_ge1=_share(found_ge1, n),                 # found and needing at least one reaction / distinct
        routed=_share(cnt["routed"], n),
        direct=_share(cnt["direct_purchase"], n),
        unresolved=_share(cnt["unresolved"], n),
        non_stock_route_rate=_share(found_ge1, found_ge1 + cnt["unresolved"]),
        routed_per_1000_raw=_share(cnt["routed"] * 1000, n_raw),
        n_evaluated=evaluated,
        denominators=dict(find="distinct molecules", find_raw="raw valid rows", find_ge1="distinct molecules",
                          non_stock_route_rate="found-with-reactions + unresolved"),
    )
