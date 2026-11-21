# Changelog

## Unreleased
- The output uses the measure names throughout: `report.json` keys, `molecules.csv` columns, console output, messages and the
  module names. The older identifiers P1 to P4 no longer appear.
- The report opens with a summary of the four measures. The versions of RDKit, AiZynthFinder, Rxn-INSIGHT and the other packages are
  recorded in `report.json` and `run_manifest.json`, and a note is printed when RDKit is not the 2023.09 series.
- Dependencies are pinned to the versions behind the reported results (RDKit 2023.09, numpy below 2, pandas below 3, pyarrow below
  21), so `pip install` and the conda environments give the same results. Python 3.10 to 3.12.
- The resource store defaults to `./.synspeq` when that folder exists and otherwise to a per-user folder, so a command no longer
  reports every resource as missing when it is run from another directory.
- A downloaded stock or a log file next to the planner models no longer makes the registered models invalid.
- A missing resource is reported once. The results folder no longer holds a stray `.lock` file.
- First-run progress messages for the precedent index.
- New guides for large sets, purchasing catalogs, output files and use from Python.
- Table-format HDF5 stocks are counted correctly.
