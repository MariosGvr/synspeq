"""Pillar 4 building blocks: generic Bemis-Murcko scaffolds, N50, CSR-AUC, fixed-size subsample N50, scaffold prior.

Acyclic molecules have an empty scaffold: they are excluded from scaffold statistics and counted separately.
N50 = fewest scaffolds that cover at least half of the (cyclic) molecules. CSR-AUC = mean cumulative coverage over
scaffold ranks (1 = a few scaffolds dominate). The subsample N50 is the mean over seeds of N50 on a random
fixed-size sample (Python random.Random(seed).sample over the molecules in input order).
"""
from __future__ import annotations
import math
import random
from collections import Counter
from typing import Dict, Iterable, List, Optional, Sequence
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def generic_scaffold(smiles_or_mol) -> Optional[str]:
    """Generic Bemis-Murcko scaffold SMILES; '' for acyclic molecules; None if it cannot be computed."""
    m = Chem.MolFromSmiles(smiles_or_mol) if isinstance(smiles_or_mol, str) else smiles_or_mol
    if m is None:
        return None
    try:
        return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m)))
    except Exception:
        return None


def n50_csr(scaffolds: Iterable[Optional[str]]):
    """-> (n50, n_distinct_scaffolds, csr_auc, n_cyclic). Empty/None scaffolds are excluded."""
    c = Counter(s for s in scaffolds if s)
    total = sum(c.values())
    if total == 0:
        return 0, 0, 0.0, 0
    freqs = sorted(c.values(), reverse=True)
    cum, n50 = 0, 0
    for i, f in enumerate(freqs, 1):
        cum += f
        if cum >= 0.5 * total:
            n50 = i
            break
    cum, ys = 0, []
    for f in freqs:
        cum += f
        ys.append(cum / total)
    return n50, len(c), sum(ys) / len(ys), total


def subsample_n50(scaffolds: Sequence[Optional[str]], n: int = 3000, seeds: Sequence[int] = (42, 43, 44, 45, 46)):
    """Mean N50 over seeded subsamples of n molecules; None when fewer than n molecules (not an n-sample estimate)."""
    sc = list(scaffolds)
    if len(sc) < n:
        return None
    vals = [n50_csr(random.Random(s).sample(sc, n))[0] for s in seeds]
    return sum(vals) / len(vals)


def scaffold_prior(smiles: Iterable[str]) -> Dict[str, int]:
    """Scaffold frequency prior from a sample of reference-catalogue SMILES (e.g. a random draw from REAL parts)."""
    c: Counter = Counter()
    for s in smiles:
        g = generic_scaffold(s)
        if g is not None:
            c[g] += 1
    return dict(c)


def cosine_to_prior(scaffolds: Iterable[Optional[str]], prior: Dict[str, int]) -> float:
    m = Counter(s for s in scaffolds if s)
    keys = set(m) | set(prior)
    dot = sum(m.get(k, 0) * prior.get(k, 0) for k in keys)
    nm = math.sqrt(sum(v * v for v in m.values()))
    npr = math.sqrt(sum(v * v for v in prior.values()))
    return dot / (nm * npr) if nm and npr else 0.0
