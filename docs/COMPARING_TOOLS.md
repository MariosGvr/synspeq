# Comparing tools with SYNSPEQ

SYNSPEQ scores any file of molecules, however the molecules were produced. It does not know how a set was generated,
filtered or selected, so a fair comparison between tools depends on choices made before and after the file reaches it.
This page lists those choices. Results are comparable only when the points below hold.

## 1. Fix the resources
Use the same planner stock, planner settings and purchasing catalogue for every tool. SYNSPEQ records them with every run
and refuses to reuse a results folder whose stock, model files or settings have changed, so two runs that differ in them
cannot be mixed by accident. Give each configuration its own `--out`.

## 2. Compare sets of the same declared size, and say how each was made
- State how many molecules each tool generated per target, and report both the raw and the distinct count. SYNSPEQ reports
  raw rows, distinct molecules and routed molecules for every run.
- State any filtering or selection applied before scoring (a docking cutoff, a synthesizability filter, taking the first
  or best N). SYNSPEQ cannot see this, so it must come from the authors. A rate on a filtered set is not comparable with a
  rate on an unfiltered one; in our benchmark, keeping the best quarter by docking score moved findability by -15.3
  percentage points for one tool and +12.4 for another.
- Declare the budget before looking at results. We do not prescribe a number; report an interval with every rate.

## 3. Report denominators and the purchasable share
Every rate names its denominator (raw rows, distinct molecules or routed molecules). Report separately the molecules that
are already purchasable and need no route, because including them can reorder tools.

## 4. Treat the planner as noisy
AiZynthFinder is not seeded. Re-planning the same molecules changes findability slightly (at most 1.7 percentage points in
our test on 20 populations), so differences of that size should not be interpreted. Re-evaluating saved planner records is
deterministic: SYNSPEQ reuses the records and recomputes the summaries.

## 5. Disclose
- Training data that overlaps the evaluation molecules or targets, and any tuning on the evaluation targets.
- Shared corpora between a tool and the evaluation: a generator trained on the reactions or building blocks that the planner
  and the support database also use will look better than one that was not.
- Every departure from the tool's documented settings, tabulated per row, and which rows could not be matched. Do not use
  unmatched rows in ranking statements.

## 6. Report the pillars together, and beside the other objectives
Do not collapse the four pillars into one score: they rank tools differently. Report them next to the objective the tool
was built for and other properties, so that a gain on one is visible against its cost on the others.

## 7. Missing resources are not zero
If a resource is missing, SYNSPEQ reports an explicit state (not assessed, unavailable, bound or blocked), never 0%.
Report the same in comparisons.

## Evaluating your own tool
1. Write the molecules to an SDF (or give saved AiZynthFinder JSON records with `--routes-dir`). Add a property that names
   the target if you want per-target results (`--pocket-field`).
2. `synspeq run molecules.sdf --pocket-field pocket_id --out results/<tool>/`
3. Repeat with the same resources for every tool, each in its own folder, and compare the reports.
