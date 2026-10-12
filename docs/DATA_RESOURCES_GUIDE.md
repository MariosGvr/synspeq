# SYNSPEQ data resources: where to get them and how the tool helps

SYNSPEQ ships **no** third-party data. It ships a manifest (`src/synspeq/data/resource_manifest.yaml`), importers, validators
and download help. You choose which catalogues to use; the report lists the resources used.

## Commands
The resource registry lives in `./.synspeq` (the current directory, like `.git`). Set `SYNSPEQ_DATA_DIR`, put `data_dir:` in
`config.yaml`, or pass `--data-dir` to use another store (command line, then config file, then environment, then default).
```bash
synspeq data list                      # every supported resource, how it is obtained, and where your copy is registered
synspeq data info zinc20_tranches      # source URL, access type, file format, citation, terms, warnings
synspeq data download aizynth_zinc --dest data/aizynth --accept-terms   # terms and citation are printed first
synspeq data import enamine_bb data/Enamine_Building_Blocks_Stock.sdf   # registers a file or folder you obtained yourself (nothing is copied)
synspeq data check                     # per resource: ready / missing / invalid / partial, and what each pillar can do with it
synspeq data build-union --resources enamine_bb,zinc20_tranches --out my_union/   # merge registered catalogues into one planner stock
```
`import` validates the content: a file must contain at least one usable record (SMILES and SDF are parsed with the same
readers the pipeline uses; HDF5 and Parquet must have the expected dataset or columns). A registered file that changed after
import is reported as `invalid` until it is imported again.

The registry stores absolute paths. To use a store on another machine, copy the data files and import them there.

## Which resource can the tool fetch for you?
| resource | role | where it comes from | `data download` |
|---|---|---|---|
| AiZynthFinder ZINC stock (**default**) and models | planner stock, planner models | AiZynthFinder's public data | yes, through AiZynthFinder's `download_public_data` (needs `--dest`, `--accept-terms`) |
| USPTO reaction database | Pillar 2 | Zenodo 10.5281/zenodo.10171745 | yes (`--dest`, `--accept-terms`; size and md5 verified) |
| ZINC20 tranches | planner stock | files.docking.org/2D, purchasability A and B | no: download, then `import` |
| ASKCOS buyables | planner stock | askcos-data repository (CC BY-NC-SA: non-commercial) | no |
| Enamine building blocks | purchasing catalogue (Pillar 3); planner stock if listed under `planner_stock` | Enamine catalogue page, 'Global Stock' SDF | no |
| eMolecules building blocks | planner stock | your eMolecules account or agreement | no |
| Molport | planner stock | your Molport registration | no |
| Enamine REAL | purchasing catalogue (on-demand tier) and scaffold reference | your REAL access (request to Enamine) | no; `real-scan` runs on your copy |

A catalogue that is not already an InChIKey stock is converted to one before planning (cached under `<store>/cache/stocks`).

## Profiles and planner stock
- `default`: AiZynthFinder's ZINC stock only.
- `custom`: the resources you list under `resources:`. Any planner-stock resource you list (ZINC20, ASKCOS, eMolecules,
  Molport, a union stock) forms the planner stock together, and the default stock is then **not** added.
- Enamine building blocks are a purchasing catalogue by default. They become part of the planner stock only if you list them
  under `planner_stock:` (or pass `--planner-stock enamine_bb`). Using one catalogue both as planner stock and as the
  purchasing catalogue is allowed but makes Pillar 3 find nothing new, because stock hits are already leaves.

## What you must know before you use a resource
- Terms differ per provider and change over time. `data info <id>` prints what we found and the link; you are responsible for
  complying with them. Enamine's terms in particular restrict redistribution and some uses of its data. SYNSPEQ never uploads
  or shares catalogue files or indexes.
- Results are comparable only within one stock and one set of planner settings.
- Catalogue membership is a vendor-listed snapshot, not a live quote or a guarantee of availability.
- Without REAL, REAL-dependent metrics are *not assessed* (null in `report.json`, `n/a` in the report), not zero.

## What happens when a resource is missing
The rules are in `resource_manifest.yaml` (`missing_resource_policy`, `metric_requirements`, `use_cases`).

| situation | what SYNSPEQ does | exit |
|---|---|---|
| You never mention a resource | not requested: dependent metrics read **not assessed**; the report lists them | 0 |
| You list a resource but the file or path is missing, invalid, or changed since import | stops before any computation and prints the fix; with `--on-missing skip` (or `on_missing: skip`) it continues and marks the affected metrics **unavailable** | 2 (0 with skip) |
| The planner stock or planner models are missing | stops; `skip` does not apply. It never substitutes the default stock for a stock you named | 2 |
| A planner stock is only partly available | stops: planning with it would change every P1 number | 2 |
| REAL scanned in part (some shards missing) | hits found so far are kept; every other query stays **unknown**, never a non-hit; P3 residual steps are an upper bound and the on-demand share a lower bound | 0 |
| Only the default stock | P1, P4 (and P2 with the USPTO database) run; P3 needs a catalogue beyond the planner stock | 0 |
| Enamine BB but no REAL | P3 uses the stocked tier only; "needs on-demand" is not assessed | 0 |
| REAL but no Enamine BB | frontier on the on-demand tier only; "needs on-demand" is not assessed | 0 |
| USPTO database missing | Pillar 2 blocked (unavailable with `skip`); P1, P3, P4 unaffected | 2 / 0 |
| Compute node without internet | nothing is downloaded during a run; a missing item gives the command to run on a connected machine | 2 |
| You re-run with a different stock, model set or settings into the same output folder | refused: planner records are reused only when stock content, model content and settings are identical | 2 |

A metric that cannot be assessed is `null` in `report.json` with a status, and `n/a` in the HTML report and on the console.

## Caches
- **Planner records** (`<out>/routes`) are reused when a run is repeated with the same planner settings, model files and stock
  file contents (checked by sha256). Otherwise SYNSPEQ refuses to reuse the folder; choose a new `--out`.
- **Reaction features** for Pillar 2 are reused when the database content, the RDKit, Rxn-INSIGHT and RXNMapper versions and the
  feature schema are the same. The reaction-database index is built only if a reaction is not cached.
- **Summaries** are always recomputed. Adding REAL after a first run means re-running `evaluate` with the REAL lookup.

Example `config.yaml`:
```yaml
profile: custom                  # default | custom
resources: [enamine_bb, enamine_real]    # omit anything you cannot get
planner_stock: []                # resources that form the planner stock; empty means AiZynthFinder's ZINC stock
paths:                           # optional: use a file or folder directly instead of the registry
  enamine_bb:   data/Enamine_Building_Blocks_Stock.sdf
  enamine_real: data/REAL/
on_missing: stop                 # stop (default) | skip
data_dir: .synspeq               # optional store location
```

## Using REAL (purchasing lookups and the scaffold prior)
REAL is read in one pass per shard and is not loaded as a whole; nothing is uploaded.
```bash
synspeq data import enamine_real data/REAL/                     # a folder of .smi / .csv / .sdf / .parquet parts (optionally .gz)
synspeq evaluate molecules.sdf --routes-dir json/ --resources enamine_bb,enamine_real --out results/
#   -> results/real_queries.smi : intermediates REAL has not been asked about yet (P3 is a bound until they are scanned)
synspeq real-scan --queries results/real_queries.smi --real-dir data/REAL --shard 0 --n-shards 50 --out hits_0.txt   # one task per shard (SLURM array)
synspeq real-merge --queries results/real_queries.smi --hits hits_*.txt --out real_lookup/
synspeq evaluate molecules.sdf --routes-dir json/ --resources enamine_bb,enamine_real --real-lookup real_lookup/ --out results/
synspeq scaffold-prior --real-dir data/REAL --seed 1 --out prior.tsv   # P4 reference affinity; pass --scaffold-prior prior.tsv
```
Each shard file records the query set, the catalogue (file names and sizes), the shard number and the number of shards.
`real-merge` records queried molecules as non-hits only when every shard `0..N-1` of one scan is present and the headers
agree. With fewer shards it keeps the hits and leaves everything else unknown, and prints which shards are missing. Shard files
without a header (from older versions) are refused. A lookup folder is tied to one catalogue; merging shards of a different
catalogue into it is refused.

Memory: a `.parquet` catalogue part is read as a whole column, so peak memory follows the size of the largest part. Choose
the shard count so that each shard stays within your memory budget.

### Canonical SMILES depend on the RDKit version
Catalogue matching compares canonical SMILES. Catalogues you import (`enamine_bb`, ...) are canonicalised by your own RDKit when
they are loaded, so they agree with the molecules being tested. REAL is not loaded; `real-scan` compares its strings as they are,
and they were canonicalised by whichever RDKit version produced the file. A different RDKit version can write a different
(equally valid) string for a small share of molecules, which then count as non-hits. If your RDKit differs from the one that
wrote the file, add `--recanonicalise` to `real-scan`; it parses every catalogue line and is much slower. Shard headers and
`real_lookup/meta.json` record the RDKit version and whether this option was used.
