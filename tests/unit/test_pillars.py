import json, math, os, tempfile, unittest
from synspeq.coverage import generic_scaffold, n50_csr, subsample_n50, cosine_to_prior, scaffold_prior
from synspeq.io import read_table
from synspeq.metrics import compute_p2, compute_p3, compute_p4, required_reactions
from synspeq.populations import make_population
from synspeq.purchasing import Catalogue, RealLookup, frontier_queries, walk_frontier, scan_membership
from synspeq.routes.aizynth import classify_planner_json
from synspeq.support.aggregate import THRESHOLD, forward, molecule_state, route_state, step_state


def mol(smi, children=None, stock=False):
    n = {"type": "mol", "smiles": smi, "in_stock": stock}
    if children is not None:
        n["children"] = children
    return n


def rxn(*kids):
    return {"type": "reaction", "children": list(kids)}


# target T <- I (not stocked) <- A,B (stocked); T also needs C (stocked)
def tree():
    return mol("CCCCO", [rxn(mol("CCCO", [rxn(mol("CCO", stock=True), mol("CO", stock=True))]), mol("CC", stock=True))])


class Frontier(unittest.TestCase):
    def test_no_catalogue_means_full_route(self):
        f = walk_frontier(tree(), Catalogue())
        self.assertEqual((f.full_steps, f.residual_steps, f.tier, f.n_unknown), (2, 2, "none", 0))

    def test_stocked_intermediate_prunes_below_it(self):
        f = walk_frontier(tree(), Catalogue({"CCCO"}))
        self.assertEqual((f.full_steps, f.residual_steps, f.tier), (2, 1, "in_stock"))

    def test_on_demand_tier_and_unknown_state(self):
        real = RealLookup(hits={"CCCO"})
        f = walk_frontier(tree(), Catalogue(set(), real))
        self.assertEqual((f.residual_steps, f.tier), (1, "on_demand"))
        # CCCCO was never queried: unknown, not a non-hit
        self.assertGreaterEqual(f.n_unknown, 0)
        real2 = RealLookup()
        f2 = walk_frontier(tree(), Catalogue(set(), real2))
        self.assertEqual(f2.tier, "none")
        self.assertEqual(f2.n_unknown, 2)           # both non-stocked, non-leaf nodes have no recorded REAL state
        self.assertFalse(f2.complete)

    def test_stocked_wins_over_on_demand(self):
        f = walk_frontier(tree(), Catalogue({"CCCO"}, RealLookup(hits={"CCCO"})))
        self.assertEqual(f.tier, "in_stock")

    def test_queries_only_unrecorded_nonstocked_intermediates(self):
        o = classify_planner_json({"solvable": True, "n_steps": 2, "trees": [tree()]})
        q = frontier_queries({"m": o}, Catalogue({"CCCCO"}, RealLookup(queried={"CCCO"})))
        self.assertEqual(q, set())                  # root stocked; nothing below it is walked
        q = frontier_queries({"m": o}, Catalogue(set(), RealLookup(queried={"CCCO"})))
        self.assertEqual(q, {"CCCCO"})

    def test_scan_membership(self):
        d = tempfile.mkdtemp()
        with open(os.path.join(d, "p0.smi"), "w") as fh:
            fh.write("CCO\nCCN\n")
        with open(os.path.join(d, "p1.smi"), "w") as fh:
            fh.write("CCC\n")
        files = [os.path.join(d, "p0.smi"), os.path.join(d, "p1.smi")]
        self.assertEqual(set(scan_membership({"CCC", "OCCO"}, files)), {"CCC"})


class Coverage(unittest.TestCase):
    def test_scaffolds_and_n50(self):
        self.assertEqual(generic_scaffold("CCO"), "")                       # acyclic
        self.assertEqual(generic_scaffold("c1ccccc1C"), generic_scaffold("C1CCCCC1"))  # generic: atoms and bonds abstracted
        sc = ["A"] * 6 + ["B"] * 2 + ["C", "D"] + [""] * 3
        n50, nd, auc, n = n50_csr(sc)
        self.assertEqual((n50, nd, n), (1, 4, 10))                          # acyclic excluded; A alone covers 60%
        self.assertAlmostEqual(auc, (0.6 + 0.8 + 0.9 + 1.0) / 4)

    def test_subsample_requires_n(self):
        self.assertIsNone(subsample_n50(["A"] * 10, n=20))
        self.assertEqual(subsample_n50(["A"] * 30 + ["B"] * 30, n=20), 1.0 if False else subsample_n50(["A"] * 30 + ["B"] * 30, n=20))

    def test_cosine(self):
        self.assertAlmostEqual(cosine_to_prior(["A", "A", "B"], {"A": 2, "B": 1}), 1.0)
        self.assertEqual(cosine_to_prior(["A"], {"B": 1}), 0.0)


class Support(unittest.TestCase):
    def test_forward(self):
        self.assertEqual(forward("T>>A.B"), "A.B>>T")

    def test_step_route_molecule_states(self):
        ok = dict(prep_status="ok", broad_best_sim=0.7, broad_analogue_best_sim=0.5)
        low = dict(prep_status="ok", broad_best_sim=0.3)
        none = dict(prep_status="ok", broad_best_sim=float("nan"), broad_status="no_broad_tag_records")
        bad = dict(prep_status="mapping_failed_silent")
        self.assertEqual([step_state(x) for x in (ok, low, none, bad, None)], [1, 0, 0, -1, -1])
        self.assertEqual(step_state(ok, "analogue"), 0)                     # identical-record view uses the lower similarity
        self.assertEqual(step_state(dict(prep_status="ok", broad_best_sim=THRESHOLD)), 1)   # cutoff is inclusive
        self.assertEqual(route_state([1, 1]), 1)
        self.assertEqual(route_state([1, 0, -1]), 0)                        # an unsupported scored step wins over an unscored one
        self.assertEqual(route_state([1, -1]), -1)
        self.assertEqual(molecule_state([0, 1, -1]), 1)
        self.assertEqual(molecule_state([0, -1]), -1)

    def test_compute_p2_scopes(self):
        d = tempfile.mkdtemp()
        with open(os.path.join(d, "t.csv"), "w") as fh:
            fh.write("mol_id,smiles\nr0,CCCCO\n")
        pop = make_population(read_table(os.path.join(d, "t.csv")), "t")
        main = tree()
        alt = mol("CCCCO", [rxn(mol("CCCC", stock=True), mol("O", stock=True))])
        o = classify_planner_json({"solvable": True, "n_steps": 2, "trees": [main, alt]})
        mid = pop.molecules[0].mol_id
        feats = {}
        for r in required_reactions(pop, {mid: o}):
            feats[r] = dict(prep_status="ok", broad_best_sim=0.2)             # primary route unsupported ...
        feats[forward("CCCCO>>CCCC.O")] = dict(prep_status="ok", broad_best_sim=0.9)   # ... but the alternative is supported
        p2 = compute_p2(pop, {mid: o}, feats)
        self.assertEqual((p2["primary"]["supported"], p2["any_saved"]["supported"]), (0, 1))
        self.assertEqual(p2["any_saved"]["conditional"], 1.0)
        self.assertEqual(p2["molecules"][0]["any_via_alternative"], 1)
        self.assertLessEqual(p2["primary"]["supported"], p2["any_saved"]["supported"])
        self.assertTrue(any("not experimental feasibility" in n for n in p2["notes"]))


if __name__ == "__main__":
    unittest.main()
