<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/synspeq_logo_dark.svg">
    <img src="docs/assets/synspeq_logo.svg" alt="SYNSPEQ" width="520">
  </picture>
</p>

<p align="center">
  <b>SYN</b>thesizability through <b>S</b>tandardised <b>P</b>lanner <b>E</b>valuation and <b>Q</b>uantification<br>
  SYNSPEQ reports route recovery, chemical support, deliverability and scaffold coverage.
</p>

<p align="center">
  <a href="https://github.com/MariosGvr/synspeq/actions/workflows/tests.yml"><img src="https://github.com/MariosGvr/synspeq/actions/workflows/tests.yml/badge.svg" alt="tests"></a>
  <img src="https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue" alt="python">
  <img src="https://img.shields.io/badge/licence-Apache--2.0-green" alt="licence">
  <img src="https://img.shields.io/badge/status-early%20development-orange" alt="status">
</p>

---

## Overview

SYNSPEQ evaluates the synthesizability of molecule sets supplied as SDF files. It reports route recovery, retrieved
precedent support, remaining synthesis after purchasing and scaffold coverage, with counts and denominators for each
metric. Results are comparable only when the stock and planner settings are the same.

| Pillar | The question it answers | How |
|---|---|---|
| **P1 Route recovery** | Can a retrosynthesis planner find a route that ends in purchasable building blocks? | AiZynthFinder, fixed search settings |
| **P2 Chemical support** | Does each reaction in those routes have precedent in the reaction literature? | Rxn-INSIGHT retrieval against the USPTO database |
| **P3 Deliverability** | How much of the route would still have to be run, if we buy what we can? | Purchasing frontier on Enamine building blocks and REAL |
| **P4 Coverage** | Is the set diverse, or collapsed onto a few scaffolds? | Generic Bemis–Murcko scaffolds, N50, CSR-AUC |

## Quick start

```bash
git clone https://github.com/MariosGvr/synspeq && cd synspeq
conda env create -f environment.yml        # everything, including AiZynthFinder and Rxn-INSIGHT
conda activate synspeq

synspeq data check                         # what is installed, and what each pillar still needs
synspeq data download uspto_rxninsight_db --dest data/ --accept-terms
synspeq run molecules.sdf --out results/   # plan, then evaluate
open results/evaluation/report.html
```

To evaluate saved JSON from a previous AiZynthFinder run without rerunning planning:

```bash
synspeq evaluate molecules.sdf --routes-dir aizynth_json/ --pocket-field pocket_id --out results/
```

### Example output
```
P1  Find 24.4%  Find>=1 22.7%  distinct 19,970  complete=True
P2  supported (primary) …% of routed, (any saved) …%
P3  median residual steps 4  frontier found 92.0%  needs on-demand 15.0%
P4  N50 38  unique scaffold ratio 0.197  CSR-AUC 0.840
```
This output is for a set of 20,000 raw molecules.

Outputs include the HTML report, `report.json`, per-molecule results in `molecules.csv` and resource metadata in
`run_manifest.json`. The manifest records the resources used, including their sizes and checksums.

## Reporting conventions
- **Missing resources.** Missing resources are reported as *not assessed*, *unavailable*, *bound* or *blocked*. These
  states are not reported as 0%.
- **Stock selection.** If the configured stock is missing, SYNSPEQ stops without substituting another stock. With
  `on_missing: skip` it continues with the remaining configured stocks and marks the stock set as incomplete.
- **Denominators.** Each metric identifies the population used as its denominator: raw rows, distinct molecules or
  routed molecules.
- **Retrieved precedent is not experimental feasibility.** Pillar 2 reports retrieved precedent support. "No qualifying
  precedent" means that support was not retrieved. The 0.6412 threshold is a similarity cut-off, not a probability.
- **No third-party data is shipped.** Users supply catalogues or download them after reviewing the source terms.

## Resources
By default only AiZynthFinder's public ZINC stock is used. Additional resources must be configured:

| Resource | Used for | How you get it |
|---|---|---|
| AiZynthFinder ZINC stock (default) and models | P1 | `synspeq data download aizynth_zinc` |
| USPTO reaction database (Zenodo, CC BY 4.0) | P2 | `synspeq data download uspto_rxninsight_db` (checksum verified) |
| ZINC20 tranches, ASKCOS buyables | alternative planner stocks | guided import |
| Enamine building blocks | P3 stocked tier | download from Enamine, then `synspeq data import` |
| Enamine REAL | P3 on-demand tier, P4 reference | request from Enamine, then `synspeq data import` |
| eMolecules, Molport | alternative planner stocks | your own account or agreement |

`synspeq data info <id>` prints the source, terms and citation for each resource. See
[docs/DATA_RESOURCES_GUIDE.md](docs/DATA_RESOURCES_GUIDE.md) for resource requirements and file formats. The guide also
describes missing-resource handling.

| If you… | SYNSPEQ… |
|---|---|
| never mention REAL | reports the REAL-dependent numbers as *not assessed* and runs |
| list REAL but the path is wrong | stops with the exact fix (or continues with *unavailable* under `on_missing: skip`) |
| have only part of REAL scanned | reports *bounds* and says how much was covered |
| have the default stock and the USPTO database | evaluates P1, P2 and P4; P3 needs a catalogue beyond the planner stock |

## How it works
<p align="center"><img src="docs/assets/how_it_works.svg" alt="SYNSPEQ pipeline" width="100%"></p>
P1 requires planner models and a stock. P2 requires the USPTO database. P3 requires a purchasing catalogue. P4's
scaffold metrics require no additional catalogue.

- **Identity** is the RDKit canonical isomeric SMILES as read. No cleanup is applied. Invalid records are listed with a reason.
- **Routes**: analysis uses returned routes whose leaves are all in stock. The primary route is the first such route in planner order.
- **Pillar 3** stops each branch at its first catalogue-listed molecule and counts the remaining reactions.

## Commands
| Command | What it does |
|---|---|
| `synspeq run` | plan with AiZynthFinder, then evaluate (SDF → report) |
| `synspeq plan` | planning only; resumable, splits into chunks for SLURM arrays |
| `synspeq evaluate` | all assessable pillars from cached planner records |
| `synspeq validate` | input audit: valid / invalid, raw and distinct counts |
| `synspeq data list \| info \| import \| check \| download \| build-union` | manage resources; `build-union` merges your catalogues into one planner stock |
| `synspeq real-scan` / `real-merge` | one-pass REAL lookup of route intermediates |
| `synspeq scaffold-prior` | P4 reference scaffold prior from a REAL sample |
| `synspeq report` | re-render `report.html` from `report.json` |

## Installation options
```bash
conda env create -f environment.yml        # full: planner + Pillar 2 (recommended)
conda env create -f environment-core.yml   # core: ingestion, P1 replay, P3, P4
pip install -e .                           # into an environment you already manage
```
`environment.yml` installs AiZynthFinder and Rxn-INSIGHT with pip inside the conda environment. The environment pins
RDKit 2023.09, numpy < 2 and pandas < 3 for compatibility with AiZynthFinder 4.4.1. A conda-build recipe for the
core package is in [`conda/`](conda).

## Status
Early development. All four pillars and live planning are implemented, and the package produces an HTML report.
Multi-population comparison plots and a PyPI / conda-forge release are planned.
v1 does not include analysis of generator-supplied routes or Chemformer.

## Contributing and testing
```bash
python -m unittest discover -s tests/unit -t .     # 45+ tests, no network, no third-party data
```
Issues and pull requests are welcome.

## Acknowledgements
SYNSPEQ builds on [AiZynthFinder](https://github.com/MolecularAI/aizynthfinder) (MIT),
[Rxn-INSIGHT](https://github.com/mrodobbe/Rxn-INSIGHT) (MIT), [RXNMapper](https://github.com/rxn4chemistry/rxnmapper) (MIT),
[RDKit](https://www.rdkit.org) (BSD), the USPTO reaction data distributed with Rxn-INSIGHT (CC BY 4.0) and the ZINC database.
Cite the works corresponding to the components you use. `synspeq data info <id>` lists their references.

## Licence
Copyright 2026 The SYNSPEQ authors. Released under the [Apache License 2.0](LICENSE).
The licence covers SYNSPEQ's code only. Third-party software and data (AiZynthFinder, Rxn-INSIGHT, ZINC, Enamine,
eMolecules, Molport, REAL, the USPTO database) remain under their own terms, which you must follow.
