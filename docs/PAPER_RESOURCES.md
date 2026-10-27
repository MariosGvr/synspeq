# Resources behind the paper's numbers, and how to rebuild the planner stock

Every headline result in the paper was computed with one planner stock (**Union2026**), one set of planner models and one
reaction database. None of the vendor catalogues is redistributed with SYNSPEQ or with the paper: terms differ by provider and
several restrict redistribution. This page records what was used so that you can rebuild an equivalent stock from your own
copies. You are responsible for following each provider's terms.

## Resources
| resource | role | version or date | size | how to obtain |
|---|---|---|---|---|
| AiZynthFinder 4.4.1 USPTO expansion, ring-breaker and filter policies | planner models | public release (Zenodo 7797465, 7341155; CC BY 4.0) | n/a | `synspeq data download aizynth_zinc` |
| Enamine building blocks (`enamine_bb`) | part of Union2026; purchasing catalogue for Pillar 3 | Global Stock list dated 2026-06-02 (`Enamine_Building_Blocks_Stock_1100064cmpd_20260602.sdf`) | 1,100,064 records; 1,098,867 InChIKeys | Enamine website, under its terms |
| ASKCOS buyables (`askcos_buyables`) | part of Union2026 | static snapshot | 193,321 | askcos-data repository (CC BY-NC-SA 4.0) |
| ZINC20 tranches (`zinc20_tranches`) | part of Union2026 | purchasability letters A and B, all six reactivity classes, downloaded 2026-09-21 (server files dated 2025-01-28) | 10,976,655 | files.docking.org/2D |
| eMolecules building blocks (`emolecules_bb`) | part of Union2026 | order-building-block parent structures (salts removed), file dated 2026-09-06 | 5,172,254 | provider |
| Molport (`molport`) | part of Union2026 | full catalogue, file dated 2026-06-02 | 5,827,234 | provider (registered users) |
| **Union2026** | planner stock for every headline row | built 2026-09-21 | 18,683,322 compounds in the AiZynthFinder HDF5 | rebuild, below |
| Enamine REAL (`enamine_real`) | on-demand tier of Pillar 3 only; never a planner stock | version not recorded here | about 9.6 billion | Enamine, by request |
| USPTO reaction database for Rxn-INSIGHT (`uspto_rxninsight_db`) | Pillar 2 | Zenodo 10.5281/zenodo.10171745 | 1,779,551 records | `synspeq data download uspto_rxninsight_db` |

The stock counter-screen in the paper replaces Union2026 by the ASKCOS buyables alone.

## Rebuilding Union2026
Register your copies, merge them, and list the result as the planner stock.
```bash
synspeq data import enamine_bb      path/to/Enamine_Building_Blocks_Stock.sdf
synspeq data import askcos_buyables path/to/askcos_buyables
synspeq data import zinc20_tranches path/to/zinc20_tranches
synspeq data import emolecules_bb   path/to/emolecules_orderbb
synspeq data import molport         path/to/molport
synspeq data build-union --resources enamine_bb,askcos_buyables,zinc20_tranches,emolecules_bb,molport --out union2026/
```
Then list `user_union_stock` under `planner_stock:` in `config.yaml` (see the data guide), and keep the planner settings of the
paper: 500 iterations, depth limit 10, 60 s per molecule, USPTO and ring-breaker expansion policies, USPTO filter policy.

## Checking that your stock matches
Do not expect a byte-identical file: vendor lists change, and the HDF5 layout depends on the build. Compare instead:
- the number of unique InChIKeys per constituent against the table above;
- the number of compounds in the merged stock, 18,683,322 for ours (the five lists contain 18,683,351 unique InChIKeys; the
  29-key difference between the two counts was not investigated);
- the SHA-256 that SYNSPEQ records for your planner stock with every run.

SHA-256 of our files:
- `union_2026_stock.hdf5`: `f853e37cde2854d466a67572306bd127f66c6fb6e0571da6496b63e80e6066f1`
- `Enamine_Building_Blocks_Stock_1100064cmpd_20260602.sdf`: `7aad6f30fd505d32510be31987be59f04e72c6efd4fd5d6e0001a64e03823891`
- ASKCOS buyables key list (193,321 InChIKeys): `ea4814975dba44acb1bf49897264f9dd05656e7b056ef1fe5303afd435b04125`

File-level hashes of the ZINC, eMolecules and Molport source files were not recorded.

## What the stock is and is not
Membership means a vendor listed the compound in the stated file. It is not a live quote, a stock check or a guarantee of
availability, and the lists differ in how well they document stock status. The terms of the merged stock follow its most
restrictive member (for example, share-alike and non-commercial terms propagate), so do not share a rebuilt union without
checking each source.
