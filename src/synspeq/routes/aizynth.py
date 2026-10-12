"""Cached AiZynthFinder route records -> one route representation.

Record fields used: smiles, solvable, n_steps, trees (molecule/reaction trees), error.
* A "returned route" is a tree whose leaves are ALL in stock (`leaves_ok`); the planner's `trees` list also
  holds partial trees with non-stock leaves, which are not routes.
* The primary route is the first all-in-stock tree in the returned order (verified against the planner's top tree).
* Reactions are rebuilt at molecule level as 'parent>>child1.child2' (retrosynthetic direction), duplicates
  across trees removed, order kept.
Outcome vocabulary (exactly one per molecule):
  direct_purchase   found, zero reactions (the molecule itself is in stock)
  routed            found, at least one reaction and an all-in-stock tree is saved
  found_no_saved_route  planner reports solved with >=1 reaction but no all-in-stock tree was saved (counted as found, excluded from route-based pillars)
  unresolved        planner ran to its budget and found no route
  execution_error   the planner record carries an error
  pending           no planner record for this molecule yet (partial run)
"""
from __future__ import annotations
import hashlib
import json
import os
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

OUTCOMES = ("direct_purchase", "routed", "found_no_saved_route", "unresolved", "execution_error", "pending")


def leaves_ok(n: dict) -> bool:
    ch = n.get("children") or []
    if not ch:
        return bool(n.get("in_stock"))
    return all(leaves_ok(m) for r in ch for m in (r.get("children") or []))


def reactions_of(tree: dict) -> List[str]:
    out: List[str] = []

    def walk(n):
        for r in n.get("children") or []:
            if r.get("type") == "reaction":
                kids = [c for c in (r.get("children") or []) if c.get("type") == "mol"]
                out.append(n["smiles"] + ">>" + ".".join(c["smiles"] for c in kids))
                for c in kids:
                    walk(c)
    walk(tree)
    return out


@dataclass
class RouteOutcome:
    outcome: str
    n_steps_reported: Optional[int] = None
    n_trees: int = 0
    n_instock_trees: int = 0
    routes: List[Tuple[str, ...]] = field(default_factory=list)   # distinct all-in-stock reaction tuples, primary first
    primary_route_hash: Optional[str] = None
    error: Optional[str] = None
    source: Optional[str] = None
    smiles_in_record: Optional[str] = None
    primary_tree: Optional[dict] = None          # first all-in-stock tree (the planner's top returned route)
    instock_trees: List[dict] = field(default_factory=list, repr=False)

    @property
    def found(self) -> bool:
        return self.outcome in ("direct_purchase", "routed", "found_no_saved_route")


def classify_planner_json(j: dict, source: Optional[str] = None) -> RouteOutcome:
    if j.get("error"):
        return RouteOutcome("execution_error", error=str(j["error"]), source=source, smiles_in_record=j.get("smiles"))
    trees = j.get("trees") or []
    ok = [t for t in trees if leaves_ok(t)]
    routes, seen = [], set()
    for t in ok:
        rx = tuple(reactions_of(t))
        if rx not in seen:
            seen.add(rx)
            routes.append(rx)
    solvable = bool(j.get("solvable"))
    ns = j.get("n_steps")
    ns = int(ns) if ns is not None else (len(reactions_of(trees[0])) if trees else None)
    base = dict(n_steps_reported=ns, n_trees=len(trees), n_instock_trees=len(ok), source=source, smiles_in_record=j.get("smiles"))
    if not solvable:
        return RouteOutcome("unresolved", **base)
    if ns == 0 or (routes and len(routes[0]) == 0):
        return RouteOutcome("direct_purchase", routes=routes[:1] if routes else [()], **base)
    if not routes:
        return RouteOutcome("found_no_saved_route", **base)
    h = hashlib.sha1("|".join(routes[0]).encode()).hexdigest()[:16]
    return RouteOutcome("routed", routes=routes, primary_route_hash=h, primary_tree=ok[0], instock_trees=ok, **base)


def load_outcome(path: Optional[str]) -> RouteOutcome:
    if path is None or not os.path.exists(path):
        return RouteOutcome("pending", source=path)
    try:
        with open(path) as fh:
            return classify_planner_json(json.load(fh), source=path)
    except (OSError, ValueError) as e:
        return RouteOutcome("execution_error", error=f"unreadable planner record: {e}", source=path)
