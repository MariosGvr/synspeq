"""Pillar 3 purchasing frontier on the primary route (the planner's first all-in-stock tree).

Walk from the target: at each molecule, if the catalogue lists it, stop (its subtree is bought) and record the tier;
otherwise expand its reaction and walk its precursors. Residual steps = reactions left above the frontier.
Frontier tier of the route: 'on_demand' if any frontier molecule is make-on-demand only, else 'in_stock' if any is
stocked, else 'none'. Counts and depths use the retained reaction nodes of the pruned tree.

Unknown REAL state: an intermediate that is not stocked and whose REAL state was never recorded may or may not be
buyable. It is counted (`n_unknown`); when n_unknown > 0 the residual steps are an UPPER bound and 'needs REAL'
a LOWER bound for that route.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Optional, Set
from .catalogue import Catalogue, canonical


@dataclass
class Frontier:
    full_steps: int
    full_depth: int
    residual_steps: int
    residual_depth: int
    tier: str                        # in_stock | on_demand | none
    n_unknown: int = 0               # non-stocked intermediates inside the walked region with unrecorded REAL state
    n_frontier: int = 0

    @property
    def complete(self) -> bool:
        return self.n_unknown == 0


def _count(node: dict):
    n, d = (1 if node.get("type") == "reaction" else 0), 0
    for c in node.get("children") or []:
        a, b = _count(c)
        n += a
        d = max(d, b)
    return n, d + (1 if node.get("type") == "reaction" else 0)


def walk_frontier(tree: dict, cat: Catalogue) -> Frontier:
    full, fdep = _count(tree)
    unknown: Set[str] = set()
    tiers: List[str] = []

    def walk(node):
        tier, st = cat.tier_of(node.get("smiles", ""))
        if tier is not None:
            tiers.append(tier)
            return 0, 0
        kids = node.get("children")
        if not kids:
            return 0, 0                      # terminal stock leaf of the planner, not listed in the purchasing catalogue
        if st == "unknown":
            unknown.add(node.get("smiles", ""))
        tot, dep = 1, 0
        for pre in kids[0].get("children") or []:
            a, b = walk(pre)
            tot += a
            dep = max(dep, b)
        return tot, dep + 1

    res, rdep = walk(tree)
    tier = "on_demand" if "on_demand" in tiers else ("in_stock" if "in_stock" in tiers else "none")
    return Frontier(full, fdep, res, rdep, tier, len(unknown), len(tiers))


def frontier_for_outcome(outcome, cat: Catalogue) -> Optional[Frontier]:
    return walk_frontier(outcome.primary_tree, cat) if getattr(outcome, "primary_tree", None) is not None else None


def frontier_queries(outcomes: Dict[str, object], cat: Catalogue) -> Set[str]:
    """Canonical SMILES of non-stocked, non-leaf intermediates on primary trees that REAL has not been asked about."""
    q: Set[str] = set()

    def visit(node):
        c = canonical(node.get("smiles", ""))
        if c is None:
            return
        if cat.bb is not None and c in cat.bb:
            return
        kids = node.get("children")
        if not kids:
            return
        if cat.real is None or cat.real.state(c) == "unknown":
            q.add(c)
        for pre in kids[0].get("children") or []:
            visit(pre)

    for o in outcomes.values():
        t = getattr(o, "primary_tree", None)
        if t is not None:
            visit(t)
    return q
