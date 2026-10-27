# Minimum setup

You can run SYNSPEQ on your own molecules with two free downloads and no vendor agreement. What you get is route recovery,
precedent support and scaffold coverage. Deliverability (Pillar 3) needs a purchasing catalogue that you obtain yourself.

## What each pillar needs
| pillar | needs | where it comes from | in the minimum setup |
|---|---|---|---|
| P1 Route recovery | AiZynthFinder, its planner models and a planner stock | `synspeq data download aizynth_zinc` (public; the default stock is AiZynthFinder's ZINC stock) | yes |
| P2 Chemical support | the USPTO reaction database used by Rxn-INSIGHT (and RXNMapper, installed with the environment) | `synspeq data download uspto_rxninsight_db` (Zenodo, CC BY 4.0) | yes |
| P3 Deliverability | a purchasing catalogue beyond the planner stock: Enamine building blocks, and optionally Enamine REAL for the on-demand tier | from the vendor, under its terms; register with `synspeq data import` | no: reported as *not assessed* |
| P4 Coverage | nothing for N50, CSR-AUC and the unique-scaffold ratio; the optional scaffold-profile metric needs a REAL or Enamine copy | n/a | yes, without the optional metric |

## Steps
```bash
git clone https://github.com/MariosGvr/synspeq && cd synspeq
conda env create -f environment.yml && conda activate synspeq        # or: pip install -e ".[planner,support]" -c constraints.txt
synspeq data download aizynth_zinc --dest data/aizynth --accept-terms       # planner models and default stock
synspeq data download uspto_rxninsight_db --dest data/uspto --accept-terms  # Pillar 2 reaction database
synspeq data check                                                          # what is ready, what each pillar still needs
synspeq run molecules.sdf --pocket-field pocket_id --out results/
```
About 2 GB in total. Planning takes up to 60 s per distinct molecule by default, so split large sets over processes or
array tasks (see the README). Open `results/evaluation/report.html`.

## What the report says when something is missing
Pillar 3 reads *not assessed* and the REAL-dependent values read `n/a`; nothing is shown as 0%. Add a catalogue later and
re-run `synspeq evaluate` on the saved planner records; the records are reused.

## Comparing with other work
Numbers are comparable only with runs that used the same planner stock and settings. The default stock above has a different
composition from the stock behind the paper's results; to reproduce that stock see [`PAPER_RESOURCES.md`](PAPER_RESOURCES.md). For how to
compare tools fairly see [`COMPARING_TOOLS.md`](COMPARING_TOOLS.md).
