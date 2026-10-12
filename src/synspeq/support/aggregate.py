"""Pillar 2 aggregation (pure functions; no chemistry dependencies).

Reaction supported  : the best Tanimoto of MACCS-key PRESENCE among database records sharing the query's broad
                      Rxn-INSIGHT tag (TAG2) is >= THRESHOLD. THRESHOLD = 0.6412 is a similarity cutoff chosen on the
                      FREA calibration half, not a probability. No qualifying record = support not retrieved.
Reaction unscored   : query preprocessing failed (unparseable, mapping failed, ...): state -1.
Route supported     : every step supported; unsupported if any scored step is unsupported; indeterminate (-1) if all
                      scored steps pass but one is unscored.
Molecule supported  : any evaluated route supported ('any_saved' scope) or the primary route supported ('primary').
Views               : 'all' keeps records identical to the query; 'analogue' excludes every identical record.
"""
from __future__ import annotations
import math
from typing import Dict, List, Optional, Sequence

THRESHOLD = 0.6412


def forward(retro: str) -> str:
    """'parent>>child1.child2' (retrosynthetic) -> 'child1.child2>>parent' (forward, the string retrieval is run on)."""
    left, right = retro.split(">>")
    return right + ">>" + left


def route_forward_reactions(route) -> List[str]:
    return [forward(r) for r in route]


def step_state(feat: Optional[dict], view: str = "all", threshold: float = THRESHOLD) -> int:
    """1 supported, 0 unsupported (including 'no qualifying record'), -1 unscored."""
    if feat is None or feat.get("prep_status") != "ok":
        return -1
    sim = feat.get("broad_analogue_best_sim") if (view == "analogue" and "broad_analogue_best_sim" in feat) else feat.get("broad_best_sim")
    if sim is None or (isinstance(sim, float) and math.isnan(sim)):
        return 0
    return 1 if sim >= threshold else 0


def route_state(step_states: Sequence[int]) -> int:
    if any(s == 0 for s in step_states):
        return 0
    return -1 if any(s == -1 for s in step_states) else 1


def molecule_state(route_states: Sequence[int]) -> int:
    if any(s == 1 for s in route_states):
        return 1
    return -1 if any(s == -1 for s in route_states) else 0
