# Output files

`synspeq run` writes `results/routes/` (the planner records) and `results/evaluation/`. `synspeq evaluate` writes the evaluation
folder to the `--out` you give it.

| File | Content |
|---|---|
| `report.html` | The report: a summary of the four measures, then one section per measure, the resources used and how to read the numbers. One self-contained file, no external assets. |
| `report.json` | Every number with its counts and denominator. Keys: `findability`, `precedent_support`, `residual_steps`, `scaffold_diversity`, plus `population`, `status`, `resources`, `banner` and `environment`. |
| `molecules.csv` | One row per distinct molecule. |
| `run_manifest.json` | The resources used (size, checksum, snapshot date) and the versions of SYNSPEQ, Python, RDKit and the planner. |
| `precedent_support_features.jsonl` | Cache of the reaction lookups, reused when the same run is repeated. |

## `status`
Each measure has a status in `report.json`, never a bare 0.

| Status | Meaning |
|---|---|
| `value` | Computed from complete inputs. |
| `bound` | Computed from partial inputs, for example REAL scanned in part. The number is a lower or an upper bound. |
| `not_assessed` | The measure was not selected or its resource was never named. |
| `unavailable` | A requested resource is missing and the run was told to continue (`--on-missing skip`). |
| `blocked` | The run refused to start, with exit code 2. |

## `molecules.csv` columns
`mol_id`, `pocket`, `identity` (RDKit canonical isomeric SMILES), `n_records` and `outcome` (`direct_purchase`, `routed`,
`found_no_saved_route`, `unresolved`, `execution_error` or `pending`) come first. The rest are named after the measure that
produced them: `n_steps` and `primary_route_hash` for findability, `precedent_support_*` and `residual_steps_*` for the other two.

## Versions
`environment` in `report.json` and `run_manifest.json` records the versions of SYNSPEQ, Python, RDKit, numpy, pandas, AiZynthFinder,
Rxn-INSIGHT and RXNMapper. Molecule identity is the RDKit canonical SMILES, so results made with another RDKit release can differ.
SYNSPEQ prints a note when RDKit is not the 2023.09 series behind the reported results.
