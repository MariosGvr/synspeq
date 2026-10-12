"""Cached AiZynthFinder route records -> one route representation.

A planner record is {mol_id, smiles, solvable, n_steps, trees, error}. It is validated before it is classified:
a malformed record is an `execution_error` with a reason, never a chemical outcome.

* A "returned route" is a tree whose leaves are ALL in stock (`leaves_ok`); the planner's `trees` list also holds
  partial trees with non-stock leaves, which are not routes.
* The primary route is the first all-in-stock tree in the returned order.
* Reactions are rebuilt at molecule level as 'parent>>child1.child2' (retrosynthetic direction), duplicates across
  trees removed, order kept.
* If the caller states which molecule the record must describe (`expected_identity`), the record's own target is
  compared with it after canonicalisation; a record for another molecule is rejected.

Outcome vocabulary (exactly one per molecule):
  direct_purchase       found, zero reactions (a single in-stock molecule)
  routed                found, at least one reaction and an all-in-stock tree is saved
  found_no_saved_route  planner reports solved with >=1 reaction but no all-in-stock tree was saved
  unresolved            the planner finished and found no route
  execution_error       the record carries an error, is malformed, or describes another molecule
  pending               no planner record for this molecule yet
"""
from __future__ import annotations
import hashlib
import json
import os
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

OUTCOMES = ("direct_purchase", "routed", "found_no_saved_route", "unresolved", "execution_error", "pending")
MAX_DEPTH = 100


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


def _check_tree(n, depth=0) -> Optional[str]:
    """Molecule nodes alternate with reaction nodes; reactions have at least one reactant; leaves say whether they are in stock."""
    if depth > MAX_DEPTH:
        return f"tree deeper than {MAX_DEPTH} levels"
    if not isinstance(n, dict):
        return "tree node is not an object"
    if n.get("type") != "mol":
        return f"expected a molecule node, got type={n.get('type')!r}"
    if not isinstance(n.get("smiles"), str) or not n["smiles"]:
        return "molecule node without a SMILES"
    ch = n.get("children")
    if ch is None or ch == []:
        if not isinstance(n.get("in_stock"), bool):
            return "leaf molecule without a boolean in_stock"
        return None
    if not isinstance(ch, list):
        return "children is not a list"
    for r in ch:
        if not isinstance(r, dict) or r.get("type") != "reaction":
            return "a molecule's child is not a reaction node"
        kids = r.get("children")
        if not isinstance(kids, list) or not kids:
            return "reaction node without reactants"
        for k in kids:
            e = _check_tree(k, depth + 1)
            if e:
                return e
    return None


def validate_record(j) -> Optional[str]:
    """None if the record is well formed, else a short reason."""
    if not isinstance(j, dict):
        return f"record is {type(j).__name__}, not an object"
    if "solvable" not in j or not isinstance(j["solvable"], bool):
        return "missing or non-boolean 'solvable' (planner did not record a completed search)"
    trees = j.get("trees", [])
    if trees is None:
        trees = []
    if not isinstance(trees, list):
        return "'trees' is not a list"
    ns = j.get("n_steps")
    if ns is not None and (isinstance(ns, bool) or not isinstance(ns, int) or ns < 0):
        return "'n_steps' is not a non-negative integer"
    for t in trees:
        e = _check_tree(t)
        if e:
            return e
    if j["solvable"] and not trees:
        return "solvable record without any tree"
    return None


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


def _canon(smiles) -> Optional[str]:
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) and smiles else None
    return Chem.MolToSmiles(m) if m is not None else None


def classify_planner_json(j: dict, source: Optional[str] = None, expected_identity: Optional[str] = None) -> RouteOutcome:
    if isinstance(j, dict) and j.get("error"):
        return RouteOutcome("execution_error", error=str(j["error"]), source=source, smiles_in_record=j.get("smiles"))
    why = validate_record(j)
    if why:
        return RouteOutcome("execution_error", error=f"invalid planner record: {why}", source=source,
                            smiles_in_record=j.get("smiles") if isinstance(j, dict) else None)
    if j.get("error"):
        return RouteOutcome("execution_error", error=str(j["error"]), source=source, smiles_in_record=j.get("smiles"))
    trees = j.get("trees") or []
    if expected_identity is not None:
        targets = {_canon(j.get("smiles"))} if j.get("smiles") else set()
        targets |= {_canon(t.get("smiles")) for t in trees}
        targets.discard(None)
        if not targets:
            return RouteOutcome("execution_error", error="planner record does not state its target molecule", source=source)
        if targets != {expected_identity}:
            return RouteOutcome("execution_error", source=source, smiles_in_record=j.get("smiles"),
                                error=f"planner record is for another molecule ({sorted(targets)[0]} instead of {expected_identity})")
    ok = [t for t in trees if leaves_ok(t)]
    routes, seen = [], set()
    for t in ok:
        rx = tuple(reactions_of(t))
        if rx not in seen:
            seen.add(rx)
            routes.append(rx)
    solvable = bool(j.get("solvable"))
    ns = j.get("n_steps")
    base = dict(n_steps_reported=ns, n_trees=len(trees), n_instock_trees=len(ok), source=source, smiles_in_record=j.get("smiles"))
    if not solvable:
        return RouteOutcome("unresolved", **base)
    if routes and len(routes[0]) == 0:
        return RouteOutcome("direct_purchase", routes=routes[:1], primary_tree=ok[0], instock_trees=ok, **base)
    if not routes:
        return RouteOutcome("found_no_saved_route", **base)
    h = hashlib.sha1("|".join(routes[0]).encode()).hexdigest()[:16]
    return RouteOutcome("routed", routes=routes, primary_route_hash=h, primary_tree=ok[0], instock_trees=ok, **base)


def load_outcome(path: Optional[str], expected_identity: Optional[str] = None) -> RouteOutcome:
    if path is None or not os.path.exists(path):
        return RouteOutcome("pending", source=path)
    try:
        with open(path) as fh:
            j = json.load(fh)
        return classify_planner_json(j, source=path, expected_identity=expected_identity)
    except (OSError, ValueError, RecursionError) as e:
        return RouteOutcome("execution_error", error=f"unreadable planner record: {type(e).__name__}: {e}", source=path)
