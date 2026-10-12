"""Pillar 4: coverage (model level), on the routed molecules and, beside it, on all distinct valid molecules."""
from __future__ import annotations
from statistics import median
from typing import Dict, Optional
from ..coverage import generic_scaffold, n50_csr, subsample_n50, cosine_to_prior
from ..populations import Population


def _block(mols, prior: Optional[dict], n_sub: int, pocket_min: int = 20):
    sc = [generic_scaffold(m.identity) for m in mols]
    n = len(sc)
    n50, nd, csr, n_cyc = n50_csr(sc)
    by: Dict[Optional[str], list] = {}
    for m, s in zip(mols, sc):
        by.setdefault(m.pocket, []).append(s)
    per_pocket = [n50_csr(v)[0] for k, v in by.items() if k is not None and len(v) >= pocket_min]
    return dict(
        n_molecules=n, n_cyclic=n_cyc, n_acyclic=n - n_cyc - sum(1 for s in sc if s is None), n_scaffold_failures=sum(1 for s in sc if s is None),
        n_distinct_scaffolds=nd,
        unique_scaffold_ratio=(nd / n) if n else None,                       # distinct scaffolds / all molecules in the set
        unique_scaffold_ratio_cyclic=(nd / n_cyc) if n_cyc else None,         # excluding acyclic molecules
        n50=n50, csr_auc=csr,
        n50_subsample_mean=subsample_n50(sc, n_sub), subsample_n=n_sub,
        n50_per_pocket_median=median(per_pocket) if per_pocket else None, n_pockets_with_min=len(per_pocket), pocket_min=pocket_min,
        reference_affinity=(cosine_to_prior(sc, prior) if prior else None),
    )


def compute_p4(pop: Population, outcomes: Dict[str, object], prior: Optional[dict] = None, n_sub: int = 3000) -> dict:
    routed = [m for m in pop.molecules if outcomes.get(m.mol_id) is not None and outcomes[m.mol_id].outcome in ("routed", "found_no_saved_route")]
    return dict(label=pop.label,
                routed=_block(routed, prior, n_sub),
                all_distinct=_block(pop.molecules, prior, n_sub),
                scope="routed = found and needing at least one reaction",
                note="subsample N50 is only reported when the set has at least the subsample size; 'reference_affinity' needs a scaffold prior from a reference catalogue")
