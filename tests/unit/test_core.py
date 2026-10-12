import gzip, json, os, tempfile, unittest
from rdkit import Chem
from synspeq.identity import identity_from_smiles
from synspeq.io import read_sdf, read_table
from synspeq.populations import make_population
from synspeq.routes.aizynth import classify_planner_json
from synspeq.metrics import compute_p1


def mol_block(smiles, name):
    m = Chem.MolFromSmiles(smiles)
    m.SetProp("_Name", name)
    return Chem.MolToMolBlock(m)


def write_sdf(path, entries):
    """entries: list of (name, smiles | raw molblock, {props})."""
    txt = ""
    for name, smi, props in entries:
        blk = smi if smi.startswith("\n") or "M  END" in smi else mol_block(smi, name)
        txt += blk
        for k, v in props.items():
            txt += f">  <{k}>\n{v}\n\n"
        txt += "$$$$\n"
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "wt") as fh:
        fh.write(txt)


def _bad_valence_block():
    m = Chem.MolFromSmiles("C[N](C)(C)(C)C", sanitize=False)      # pentavalent neutral nitrogen: RDKit refuses to sanitise it
    m.SetProp("_Name", "bad")
    return Chem.MolToMolBlock(m, kekulize=False)


BAD_VALENCE = _bad_valence_block()


class Ingestion(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.sdf = os.path.join(self.d, "m.sdf")
        write_sdf(self.sdf, [
            ("a", "CCO", {"pocket": "P1"}),
            ("a", "OCC", {"pocket": "P1"}),            # same molecule, same name (names may repeat), duplicate in P1
            ("b", "OCC", {"pocket": "P2"}),            # same molecule, different pocket: kept
            ("c", "C[C@H](N)C(=O)O", {"pocket": "P1"}),
            ("d", "C[C@@H](N)C(=O)O", {"pocket": "P1"}),  # enantiomer: different identity
            ("e", "CC(=O)[O-].[Na+]", {"pocket": "P1"}),  # salt kept as given
            ("bad", BAD_VALENCE, {"pocket": "P1"}),
        ])

    def test_records_and_reasons(self):
        recs = list(read_sdf(self.sdf))
        self.assertEqual(len(recs), 7)
        bad = recs[-1]
        self.assertEqual(bad.status, "invalid")
        self.assertTrue(bad.reason.startswith("sanitize_failed") or bad.reason in ("no_atoms", "molblock_unreadable"), bad.reason)
        self.assertEqual(len({r.record_id for r in recs}), 7)
        self.assertEqual(recs[0].props["pocket"], "P1")

    def test_pocket_grouping(self):
        pop = make_population(read_sdf(self.sdf), "x", pocket_field="pocket")
        acc = pop.accounting()
        self.assertEqual(acc["n_input_records"], 7)
        self.assertEqual(acc["n_invalid_records"], 1)
        self.assertEqual(acc["n_valid_records"], 6)
        self.assertEqual(acc["n_distinct_molecules"], 5)       # ethanol P1, ethanol P2, L-Ala, D-Ala, sodium acetate
        self.assertEqual(acc["n_duplicate_records"], 1)
        self.assertAlmostEqual(acc["duplicate_rate"], 1 / 6)
        self.assertEqual(acc["n_multi_fragment"], 1)
        self.assertEqual(acc["n_pockets"], 2)

    def test_global_dedup_without_pocket(self):
        pop = make_population(read_sdf(self.sdf), "x")
        self.assertEqual(pop.accounting()["n_distinct_molecules"], 4)   # ethanol collapses across pockets

    def test_first_record_is_source(self):
        pop = make_population(read_sdf(self.sdf), "x", pocket_field="pocket")
        eth = [m for m in pop.molecules if m.identity == "CCO" and m.pocket == "P1"][0]
        self.assertEqual(eth.n_records, 2)
        self.assertEqual(eth.source_record.index, 0)

    def test_gzip_and_input_order_invariance_of_counts(self):
        gz = os.path.join(self.d, "m.sdf.gz")
        entries = [("a", "CCO", {"pocket": "P1"}), ("b", "CCN", {"pocket": "P1"}), ("c", "OCC", {"pocket": "P1"})]
        write_sdf(gz, entries)
        a = make_population(read_sdf(gz), "g", "pocket").accounting()
        write_sdf(gz, entries[::-1])
        b = make_population(read_sdf(gz), "g", "pocket").accounting()
        for k in ("n_distinct_molecules", "n_duplicate_records", "n_valid_records"):
            self.assertEqual(a[k], b[k])

    def test_window_first_applied_before_dedup(self):
        pop = make_population(read_sdf(self.sdf), "w", pocket_field="pocket", first=2)
        # P1 rows: a, a, c, d, e, bad -> first 2 = a, a ; P2 rows: b
        acc = pop.accounting()
        self.assertEqual(acc["n_input_records"], 3)
        self.assertEqual(acc["n_distinct_molecules"], 2)
        self.assertIn("before deduplication", acc["window"])

    def test_smiles_table(self):
        p = os.path.join(self.d, "t.csv")
        open(p, "w").write("mol_id,smiles,pocket_id\nx1,CCO,A\nx2,not_a_smiles,A\n")
        recs = list(read_table(p))
        self.assertEqual([r.status for r in recs], ["ok", "invalid"])
        self.assertEqual(recs[1].reason, "unparsable_smiles")

    def test_identity_keeps_stereo_and_salts(self):
        self.assertNotEqual(identity_from_smiles("C[C@H](N)C(=O)O")[0], identity_from_smiles("C[C@@H](N)C(=O)O")[0])
        self.assertIn(".", identity_from_smiles("CC(=O)[O-].[Na+]")[0])


def mol(smi, children=None, stock=False):
    n = {"type": "mol", "smiles": smi, "in_stock": stock}
    if children is not None:
        n["children"] = children
    return n


def rxn(*kids):
    return {"type": "reaction", "children": list(kids)}


class Routes(unittest.TestCase):
    def test_direct(self):
        o = classify_planner_json({"smiles": "CCO", "solvable": True, "n_steps": 0, "trees": [mol("CCO", stock=True)]})
        self.assertEqual(o.outcome, "direct_purchase")

    def test_routed_primary_is_first_instock_tree_not_first_tree(self):
        partial = mol("CCCC", [rxn(mol("CCO", stock=True), mol("CCN", stock=False))])             # non-stock leaf: not a route
        good = mol("CCCC", [rxn(mol("CCO", stock=True), mol("CCC", stock=True))])
        o = classify_planner_json({"smiles": "CCCC", "solvable": True, "n_steps": 1, "trees": [partial, good]})
        self.assertEqual(o.outcome, "routed")
        self.assertEqual(o.n_trees, 2)
        self.assertEqual(o.n_instock_trees, 1)
        self.assertEqual(o.routes[0], ("CCCC>>CCO.CCC",))

    def test_unresolved_error_and_inconsistent(self):
        self.assertEqual(classify_planner_json({"solvable": False, "trees": []}).outcome, "unresolved")
        self.assertEqual(classify_planner_json({"error": "boom"}).outcome, "execution_error")
        t = mol("CCCC", [rxn(mol("CCN", stock=False))])
        self.assertEqual(classify_planner_json({"solvable": True, "n_steps": 1, "trees": [t]}).outcome, "found_no_saved_route")

    def test_duplicate_trees_collapse(self):
        t = mol("CCCC", [rxn(mol("CCO", stock=True))])
        o = classify_planner_json({"solvable": True, "n_steps": 1, "trees": [t, t]})
        self.assertEqual(len(o.routes), 1)


class P1(unittest.TestCase):
    def test_partition_and_raw_view(self):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "t.csv")
        rows = [("r0", "CCO"), ("r1", "OCC"), ("r2", "CCN"), ("r3", "CCC"), ("r4", "CCCC"), ("r5", "CCCCC")]
        open(p, "w").write("mol_id,smiles\n" + "\n".join(f"{a},{b}" for a, b in rows) + "\n")
        pop = make_population(read_table(p), "t")
        by = {m.source_record.name: m.mol_id for m in pop.molecules}
        from synspeq.routes.aizynth import RouteOutcome
        out = {by["r0"]: RouteOutcome("direct_purchase"), by["r2"]: RouteOutcome("routed"), by["r3"]: RouteOutcome("unresolved"),
               by["r4"]: RouteOutcome("routed"), by["r5"]: RouteOutcome("execution_error")}
        r = compute_p1(pop, out)
        self.assertEqual(r["n_distinct"], 5)                  # CCO and OCC collapse
        self.assertEqual(r["n_raw_valid"], 6)
        self.assertAlmostEqual(r["find"], 3 / 5)
        self.assertAlmostEqual(r["find_ge1"], 2 / 5)
        self.assertAlmostEqual(r["find_raw"], 4 / 6)         # the duplicate inherits the direct-purchase outcome
        self.assertFalse(r["complete"])                       # one execution error
        self.assertEqual(r["counts"]["pending"], 0)
        self.assertAlmostEqual(r["non_stock_route_rate"], 2 / 3)
        self.assertLessEqual(r["find_ge1"], r["find"])


if __name__ == "__main__":
    unittest.main()
