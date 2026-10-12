# SYNSPEQ data resources: where to get them and how the tool helps (user-facing draft)

SYNSPEQ ships **no** third-party data. It ships a manifest (`resource_manifest.yaml`), importers, validators and download help. You choose which catalogues to use; the report always states which ones were used.

## Quick start
The resource registry lives in `./.synspeq` (the current directory, like `.git`); set `SYNSPEQ_DATA_DIR` or pass `--data-dir` to share one store between projects.
```bash
synspeq data list                      # every supported resource, whether it is enabled, and where it comes from
synspeq data info zinc20_tranches      # source URL, access type, file format, citation, terms, warnings
synspeq data download aizynth_zinc     # shows the terms and citation first; fetches only the AiZynthFinder items (needs --dest and --accept-terms)
synspeq data import enamine_bb data/Enamine_Building_Blocks_Stock.sdf   # registers a file or directory you obtained yourself (nothing is copied or shared)
synspeq data check                     # counts, invalid rows, duplicates, format, snapshot dates, hashes
synspeq data build-union --resources enamine_bb,zinc20_tranches --out my_union/   # merge registered catalogues into one planner stock + manifest
```
Compute nodes without internet: run `data download` and `data import` on a machine with network access, then copy the resulting store.

## Which catalogue can the tool fetch for you?
| resource | role | how you get it | `data download` |
|---|---|---|---|
| AiZynthFinder ZINC stock (**default**) | planner stock | `download_public_data` | yes |
| ZINC20 tranches | planner stock | files.docking.org/2D, purchasability A and B | guided (automatic download planned) |
| ASKCOS buyables | planner stock | askcos-data repository | guided (CC BY-NC-SA: non-commercial) |
| Enamine building blocks | planner stock, purchasing | Enamine catalogue page, 'Global Stock' SDF | no: you download, then `import` |
| eMolecules building blocks | planner stock | your eMolecules account/agreement | no |
| Molport | planner stock | your Molport registration | no |
| Enamine REAL | purchasing only | your REAL access (request to Enamine) | no; the tool runs the membership scan on your copy |
| USPTO reaction database | Chemical Support | Zenodo 10.5281/zenodo.10171745 | yes (`--dest`, `--accept-terms`; size and md5 verified) |
| AiZynthFinder models | route recovery | Zenodo 7797465 / 7341155 | yes, via AiZynthFinder's `download_public_data` |

## Profiles
- `default`: AiZynthFinder's ZINC stock only.
- `custom`: any subset you enable (`resources:` in `config.yaml`).

## What you must know before you use a resource
- Terms differ per provider and change over time. `data info <id>` prints what we found and the link; you are responsible for complying with them. Enamine's terms in particular restrict redistribution and some uses of its data; SYNSPEQ never uploads or shares catalogue files or indexes.
- Results are not comparable across stocks. Compare tools only within one stock and one set of planner settings.
- Catalogue membership is a vendor-listed snapshot, not a live quote or a guarantee of availability.
- REAL: if you do not provide it, REAL-dependent metrics are "unavailable", not zero.

## What happens when a resource is missing
The rules are defined in `resource_manifest.yaml` (`missing_resource_policy`, `metric_requirements`, `use_cases`). In short:

| situation | what SYNSPEQ does | exit |
|---|---|---|
| You never mention a resource | not requested: dependent metrics read **not assessed**, the report lists them and how to enable them | 0 |
| You list a resource but the file/path is missing or invalid | error at `data check` and at run start, with the fix; `on_missing: skip` (or `--allow-missing`) continues and marks the metrics **unavailable** | 2 (0 with skip) |
| The stock you named is missing | stops. It does **not** fall back to the default stock, because that would change every P1 number | 2 |
| REAL is only partly available (shards missing, scan interrupted) | REAL hits are a **lower bound**, unqueried molecules are *unknown* (never non-hit), residual steps an **upper bound**; coverage is printed | 0 |
| Nothing but the default stock | P1 and P2 run; P3 (needs a catalogue beyond the planner stock) and Enamine affinity (P4) are not assessed | 0 |
| Enamine BB but no REAL | P3 uses the stocked tier only, labelled as such; "REAL needed" not assessed | 0 |
| REAL but no Enamine BB | frontier on the on-demand tier only; no stocked/on-demand split | 0 |
| USPTO database missing | Pillar 2 blocked/unavailable; P1, P3, P4 unaffected | 2 / 0 |
| Offline compute node | the run uses the local store only; missing items give the exact command to run on a connected machine | 2 |
| You add REAL after a first run | only the REAL scan and the REAL-dependent P3 parts re-run; earlier results are kept | 0 |
| You compare reports from different stocks | refused unless `--force`, with a warning | 2 |

Absence is never printed as zero or left blank: a cell reads `n/a (REAL not configured)` or `n/a (REAL missing: path not found)`.

Example `config.yaml`:
```yaml
profile: default                 # default | custom
resources: [aizynth_zinc, enamine_bb, enamine_real]   # omit anything you cannot get
paths:
  enamine_bb:   data/Enamine_Building_Blocks_Stock.sdf
  enamine_real: data/REAL/
on_missing: stop                 # stop (default) | skip
auto_download: false             # true only for resources with automated downloads, after you accepted their terms
```

## Using REAL (purchasing lookups and the scaffold prior)
REAL is too large to load, so SYNSPEQ uses it in two one-pass scans over *your* copy; nothing is uploaded.
```bash
synspeq data import enamine_real data/REAL/                     # a folder of .parquet / .smi / .csv(.gz) parts
synspeq evaluate molecules.sdf --routes-dir json/ --resources enamine_bb,enamine_real --out results/
#   -> results/real_queries.smi : intermediates REAL has not been asked about yet (P3 reports bounds until you scan)
synspeq real-scan --queries results/real_queries.smi --real-dir data/REAL --shard 0 --n-shards 50 --out hits_0.txt   # one task per shard (SLURM array)
synspeq real-merge --queries results/real_queries.smi --hits hits_*.txt --out real_lookup/
synspeq evaluate molecules.sdf --routes-dir json/ --resources enamine_bb,enamine_real --real-lookup real_lookup/ --out results/
synspeq scaffold-prior --real-dir data/REAL --seed 1 --out prior.tsv   # P4 reference affinity; pass --scaffold-prior prior.tsv
```
Every queried molecule is recorded as hit or non-hit; a molecule that was never queried is *unknown*, never a non-hit.
A pillar-3 result that depends on unknown molecules is reported as a bound (residual steps upper, REAL-needed lower).
`.parquet` parts need `pip install pyarrow`.
