# Large sets

Planning is the slow step: AiZynthFinder searches for up to 60 s per distinct molecule. Precedent support adds a fixed cost on the
first run, because it builds an index of the 1.5 million database reactions and caches it in the resource store. Later runs reuse it and need about half a minute for the lookup setup, whatever the size of the set.

## One machine
`--workers N` plans N molecules at a time. Planning stops where it is interrupted and continues when the same command is run again.

```bash
synspeq run molecules.sdf --pocket-field pocket_id --workers 8 --out results/
```

## A SLURM array
Split the planning over array tasks with `--chunk` and `--n-chunks`, then evaluate once on the finished planner records.

```bash
#!/bin/bash
#SBATCH --array=0-49
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH -t 04:00:00
export SYNSPEQ_DATA_DIR=/path/to/store          # the store where `synspeq data download` registered the resources
synspeq plan molecules.sdf --pocket-field pocket_id --out results/routes \
    --workers 4 --chunk $SLURM_ARRAY_TASK_ID --n-chunks 50
```

```bash
synspeq evaluate molecules.sdf --pocket-field pocket_id --routes-dir results/routes --out results/evaluation
```

Evaluation reads the planner records of every chunk. A molecule whose record is missing is reported as *pending*, and the report
says the population is incomplete.

## Comparing several runs
Use the same planner stock and settings for every run, and keep each run in its own `--out` folder. SYNSPEQ refuses to mix planner
records made under another stock, model set or settings in one folder.
