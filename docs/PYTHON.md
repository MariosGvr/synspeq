# Using the results from Python

SYNSPEQ is a command-line tool, and its output files are the interface. The files are plain JSON and CSV, so they load directly.

```python
import json
import pandas as pd

report = json.load(open("results/evaluation/report.json"))
print(report["findability"]["find"], report["findability"]["find_ge1"])     # Solved, Solved>=1 (fractions)
print(report["status"]["residual_steps"]["status"])                          # value, bound or not_assessed
molecules = pd.read_csv("results/evaluation/molecules.csv")
unresolved = molecules[molecules.outcome == "unresolved"]
```

To run a command from a script, call the entry point. It returns the exit code (0 ok, 2 when something is missing or blocked).

```python
from synspeq.cli import main

rc = main(["evaluate", "molecules.sdf", "--routes-dir", "results/routes", "--pocket-field", "pocket_id", "--out", "results/evaluation"])
```

The functions inside the package (`synspeq.metrics`, `synspeq.pipeline`) are not a stable interface yet and can change before a release.
See [`OUTPUTS.md`](OUTPUTS.md) for every file and column.
