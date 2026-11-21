# Purchasing catalogs and residual steps

Residual steps count the reactions of a route that are still left after everything purchasable is bought. They need a purchasing
catalog beyond the planner stock, for example Enamine building blocks. Without one, residual steps are reported as *not assessed*
and nothing is shown as 0%.

## Enamine building blocks (stocked tier)
Download the catalog from Enamine under its terms, then register it. The building blocks are used as a purchasing catalog only. The
planner stock stays AiZynthFinder's ZINC stock unless you list them under `--planner-stock`.

```bash
synspeq data import enamine_bb data/Enamine_Building_Blocks_Stock.sdf
synspeq run molecules.sdf --pocket-field pocket_id --resources enamine_bb --out results/
```

## Enamine REAL (on-demand tier)
REAL is too large to load into memory, so SYNSPEQ looks route intermediates up in one pass over its files, split into shards that
can run as a SLURM array. Residual steps that depend on the on-demand tier are reported as a bound until every shard has been
scanned. The commands are in [`DATA_RESOURCES_GUIDE.md`](DATA_RESOURCES_GUIDE.md#using-real-purchasing-lookups-and-the-scaffold-prior).

## Reusing saved planner records
Planner records in `results/routes` and the reaction features used for precedent support are reused when a run is repeated with the
same stock, model files and settings. The summaries are recomputed. If the stock or the models change, SYNSPEQ refuses to reuse the
folder. Give it a new `--out`.

To evaluate saved JSON from a previous AiZynthFinder run without planning again:

```bash
synspeq evaluate molecules.sdf --routes-dir aizynth_json/ --pocket-field pocket_id --out results/
```
