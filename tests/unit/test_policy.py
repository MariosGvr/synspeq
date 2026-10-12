import os, tempfile, unittest
from synspeq.config import load_config
from synspeq.resources import Store, load_manifest, resolve


class PolicyBase(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.store = Store(os.path.join(self.d, "store"))
        self.man = load_manifest()

    def make(self, rid, coverage=None, fmt=None):
        p = os.path.join(self.d, f"{rid}.smi")
        with open(p, "w") as fh:
            fh.write("CCO\nCCN\nCCC\n")
        self.store.register(rid, p, fmt=fmt or "smiles", coverage=coverage)
        return p

    def core(self):
        for r in ("aizynth_zinc", "aizynth_models", "uspto_rxninsight_db"):
            self.make(r)

    def run_(self, mode="plan", **cfg):
        return resolve(load_config(None, cfg), self.man, self.store, mode=mode)


class UseCases(PolicyBase):
    def test_uc1_nothing_installed_blocks_with_the_exact_fix(self):
        r = self.run_()
        self.assertEqual(r.exit_code, 2)
        self.assertEqual(r.metrics["P1"].status, "blocked")
        self.assertTrue(any("synspeq data download aizynth_zinc --dest" in (m["fix"] or "") for m in r.messages))

    def test_uc1_default_core_present_runs_and_reports_not_assessed(self):
        self.core()
        r = self.run_()
        self.assertEqual(r.exit_code, 0)
        self.assertEqual(r.metrics["P1"].status, "value")
        self.assertEqual(r.metrics["P2"].status, "value")
        self.assertEqual(r.metrics["P3"].status, "not_assessed")          # R1: unspecified is not an error
        self.assertEqual(r.metrics["P4"].sub["enamine_affinity"]["status"], "not_assessed")
        self.assertTrue(any("P3: not_assessed" in b for b in r.banner))     # never silent

    def test_uc2_real_not_mentioned_but_bb_present(self):
        self.core(); self.make("enamine_bb", fmt="smiles")
        r = self.run_(resources=["enamine_bb"])
        self.assertEqual(r.exit_code, 0)
        self.assertEqual(r.metrics["P3"].status, "value")
        self.assertEqual(r.metrics["P3"].sub["real_needed"]["status"], "not_assessed")      # UC5

    def test_uc3_real_listed_but_missing(self):
        self.core()
        r = self.run_(resources=["enamine_real"], paths={"enamine_real": os.path.join(self.d, "nope")})
        self.assertEqual(r.exit_code, 2)
        self.assertEqual(r.metrics["P3"].status, "blocked")
        r = self.run_(resources=["enamine_real"], paths={"enamine_real": os.path.join(self.d, "nope")}, on_missing="skip")
        self.assertEqual(r.exit_code, 0)
        self.assertEqual(r.metrics["P3"].status, "unavailable")
        self.assertEqual(r.metrics["P1"].status, "value")                                   # other pillars unaffected

    def test_uc4_partial_real_gives_bounds_not_zero(self):
        self.core(); self.make("enamine_bb")
        os.makedirs(os.path.join(self.d, "REAL"), exist_ok=True)
        open(os.path.join(self.d, "REAL", "part0.smi"), "w").write("CCO\n")
        self.store.register("enamine_real", os.path.join(self.d, "REAL"), coverage={"fraction": 0.6})
        r = self.run_(resources=["enamine_bb", "enamine_real"])
        self.assertEqual(r.exit_code, 0)
        self.assertEqual(r.metrics["P3"].sub["real_needed"]["status"], "bound")
        self.assertEqual(r.states["enamine_real"].state, "partial")

    def test_uc6_real_only(self):
        self.core()
        os.makedirs(os.path.join(self.d, "REAL"), exist_ok=True)
        open(os.path.join(self.d, "REAL", "part0.smi"), "w").write("CCO\n")
        self.store.register("enamine_real", os.path.join(self.d, "REAL"))
        r = self.run_(resources=["enamine_real"])
        self.assertEqual(r.metrics["P3"].status, "value")
        self.assertEqual(r.metrics["P3"].sub["real_needed"]["status"], "not_assessed")

    def test_uc8_named_stock_missing_never_falls_back_to_default(self):
        self.core()                                        # aizynth_zinc present and valid
        r = self.run_(resources=["zinc20_tranches"])       # user named a stock that is not installed
        self.assertEqual(r.exit_code, 2)
        self.assertEqual(r.metrics["P1"].status, "blocked")
        self.assertEqual(r.states["aizynth_zinc"].state, "not_requested")
        r = self.run_(resources=["zinc20_tranches"], on_missing="skip")
        self.assertEqual(r.metrics["P1"].status, "blocked")  # a core resource is not skippable

    def test_uc9_reaction_database_missing(self):
        self.make("aizynth_zinc"); self.make("aizynth_models")
        r = self.run_()
        self.assertEqual((r.exit_code, r.metrics["P2"].status, r.metrics["P1"].status), (2, "blocked", "value"))
        r = self.run_(on_missing="skip")
        self.assertEqual((r.exit_code, r.metrics["P2"].status), (0, "unavailable"))

    def test_replay_mode_needs_no_stock_or_models(self):
        self.make("uspto_rxninsight_db")
        r = self.run_(mode="replay")
        self.assertEqual(r.exit_code, 0)
        self.assertEqual(r.metrics["P1"].status, "value")

    def test_listed_catalogue_set_with_one_missing_stops(self):
        self.core()
        for rid in ("enamine_bb", "askcos_buyables", "zinc20_tranches", "emolecules_bb"):   # molport missing
            self.make(rid)
        r = self.run_(resources=["enamine_bb", "askcos_buyables", "zinc20_tranches", "emolecules_bb", "molport"])
        self.assertEqual(r.exit_code, 2)
        self.assertTrue(any(m["scope"] == "molport" for m in r.messages))
        r = self.run_(resources=["enamine_bb", "askcos_buyables", "zinc20_tranches", "emolecules_bb", "molport"], on_missing="skip")
        self.assertEqual(r.exit_code, 0)
        self.assertIn("declared stock incomplete", r.metrics["P1"].reason)

    def test_unknown_resource_id_is_an_error(self):
        with self.assertRaises(ValueError):
            self.run_(resources=["enamine_rael"])

    def test_unknown_config_key_rejected(self):
        with self.assertRaises(ValueError):
            load_config(None, {"on_mising": "skip"})


class Union(PolicyBase):
    def test_build_union_counts_overlap_and_flags(self):
        from synspeq.resources.union import build_union
        a = os.path.join(self.d, "a.smi"); b = os.path.join(self.d, "b.smi")
        open(a, "w").write("CCO\nCCN\nnot_smiles\n")
        open(b, "w").write("OCC\nCCC\n")                      # OCC is ethanol again
        out = build_union({"a": a, "b": b}, os.path.join(self.d, "u"), flags={"a": {"redistributable": True}, "b": {"noncommercial": True}})
        self.assertEqual(out["n_unique_inchikeys"], 3)
        self.assertEqual(out["members"]["a"]["n_unparseable"], 1)
        self.assertEqual(out["members"]["b"]["overlap_with_previous"], 1)
        self.assertTrue(out["noncommercial"])
        self.assertFalse(out["redistributable"])                 # b is not redistributable
        e = self.store.register("user_union_stock", out["file"], fmt="aizynthfinder_hdf5")
        self.assertEqual(e["validation"]["status"], "ok")


class StoreTests(PolicyBase):
    def test_register_records_hash_size_and_validation(self):
        p = self.make("molport", fmt="smiles")
        e = self.store.get("molport")
        self.assertEqual(e["path"], p)
        self.assertEqual(len(e["sha256"]), 64)
        self.assertEqual(e["validation"]["status"], "ok")
        self.assertEqual(e["validation"]["n_records"], 3)

    def test_invalid_file_is_flagged(self):
        p = os.path.join(self.d, "junk.smi")
        open(p, "w").write("not a smiles\n" * 50)
        e = self.store.register("molport", p, fmt="smiles")
        self.assertEqual(e["validation"]["status"], "invalid")
        r = resolve(load_config(None, {"resources": ["molport"]}), self.man, self.store)
        self.assertEqual(r.states["molport"].state, "invalid")


if __name__ == "__main__":
    unittest.main()
