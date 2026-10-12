"""Regression tests for the independent pre-release review (findings F01-F22). Each test reproduces the review's input and
asserts the value written down BEFORE looking at the code's answer (hand-computed or taken from the stated rule)."""
import gzip, json, multiprocessing, os, re, shlex, tempfile, unittest
from unittest import mock
from rdkit import Chem
from synspeq.cli import build_parser, main
from synspeq.config import load_config
from synspeq.coverage import cosine_to_prior, n50_csr, scaffold_prior, subsample_n50
from synspeq.io import read_sdf, read_table
from synspeq.metrics import compute_p1, compute_p3, compute_p4
from synspeq.pipeline import load_outcomes
from synspeq.populations import make_population
from synspeq.purchasing import Catalogue, RealLookup, load_building_blocks
from synspeq.purchasing.realscan import merge_shards, scan_shard
from synspeq.resources import Store, load_manifest, resolve
from synspeq.routes.aizynth import classify_planner_json, load_outcome
from tests.unit.test_core import write_sdf

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def mol(smi, children=None, stock=False):
    n = {"type": "mol", "smiles": smi, "in_stock": stock}
    if children is not None:
        n["children"] = children
    return n


def rxn(*kids):
    return {"type": "reaction", "children": list(kids)}


def direct(smi, name=None):
    return {"mol_id": name or smi, "smiles": smi, "solvable": True, "n_steps": 0, "trees": [mol(smi, stock=True)], "error": None}


class Tmp(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def p(self, *a):
        return os.path.join(self.d, *a)

    def write(self, name, text, mode="w"):
        path = self.p(name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, mode) as fh:
            fh.write(text)
        return path


class PlannerRecordBinding(Tmp):
    def test_F01_title_collision_does_not_bind_the_wrong_molecule(self):
        sdf = self.p("m.sdf")
        write_sdf(sdf, [("same", "CCO", {}), ("same", "CCN", {})])
        os.makedirs(self.p("r"))
        json.dump(direct("CCO", "same"), open(self.p("r", "same.json"), "w"))
        pop = make_population(read_sdf(sdf), "t")
        out = load_outcomes(pop, [self.p("r")])
        by = {m.identity: out[m.mol_id] for m in pop.molecules}
        self.assertEqual(by["CCO"].outcome, "direct_purchase")
        self.assertEqual(by["CCN"].outcome, "execution_error")           # the record describes another molecule
        self.assertIn("another molecule", by["CCN"].error)
        r = compute_p1(pop, out)
        self.assertFalse(r["complete"])
        self.assertAlmostEqual(r["find"], 0.5)

    def test_F09_missing_pocket_and_literal_none_are_two_units(self):
        sdf = self.p("m.sdf")
        write_sdf(sdf, [("a", "CCO", {}), ("b", "CCO", {"pocket_id": "None"})])
        pop = make_population(read_sdf(sdf), "t", pocket_field="pocket_id")
        ids = [m.mol_id for m in pop.molecules]
        self.assertEqual(len(pop.molecules), 2)
        self.assertEqual(len(set(ids)), 2)
        self.assertEqual({m.pocket for m in pop.molecules}, {None, "None"})

    def test_missing_pocket_field_is_rejected_not_ignored(self):
        sdf = self.p("m.sdf")
        write_sdf(sdf, [("a", "CCO", {"pocket": "A"})])
        with self.assertRaises(ValueError):
            make_population(read_sdf(sdf), "t", pocket_field="pocket_id")      # typo

    def test_F14_title_cannot_escape_the_records_folder(self):
        sdf = self.p("m.sdf")
        write_sdf(sdf, [("../outside", "CCO", {}), (self.p("abs"), "CCN", {})])
        os.makedirs(self.p("r"))
        json.dump(direct("CCO", "outside"), open(self.p("outside.json"), "w"))
        json.dump(direct("CCN", "abs"), open(self.p("abs.json"), "w"))
        out = load_outcomes(make_population(read_sdf(sdf), "t"), [self.p("r")])
        self.assertEqual({o.outcome for o in out.values()}, {"pending"})

    def test_title_lookup_can_be_switched_off(self):
        sdf = self.p("m.sdf")
        write_sdf(sdf, [("ethanol", "CCO", {})])
        os.makedirs(self.p("r"))
        json.dump(direct("CCO", "ethanol"), open(self.p("r", "ethanol.json"), "w"))
        pop = make_population(read_sdf(sdf), "t")
        self.assertEqual(list(load_outcomes(pop, [self.p("r")]).values())[0].outcome, "direct_purchase")
        self.assertEqual(list(load_outcomes(pop, [self.p("r")], use_titles=False).values())[0].outcome, "pending")


class RecordSchema(Tmp):
    def test_F10_malformed_records_are_errors_not_chemistry(self):
        cases = {
            "solvable_string": {"smiles": "CCO", "solvable": "false", "n_steps": 0, "trees": []},
            "reaction_without_reactants": {"smiles": "CCO", "solvable": True, "n_steps": 1, "trees": [mol("CCO", [rxn()])]},
            "empty_object": {},
            "solvable_without_trees": {"smiles": "CCO", "solvable": True, "n_steps": 0, "trees": []},
            "leaf_without_stock_flag": {"smiles": "CCO", "solvable": True, "n_steps": 0, "trees": [{"type": "mol", "smiles": "CCO"}]},
        }
        for name, rec in cases.items():
            o = classify_planner_json(rec)
            self.assertEqual(o.outcome, "execution_error", name)
            self.assertTrue(o.error, name)

    def test_F15_wrong_json_shapes_do_not_raise(self):
        for i, text in enumerate(['[]', 'null', '"x"', '{"solvable": true, "trees": 1}', '{"solvable": false, "trees": [null]}',
                                  '{"solvable": false, "trees": [], "n_steps": []}', '{"solvable": fal', '']):
            f = self.write(f"c{i}.json", text)
            o = load_outcome(f)
            self.assertEqual(o.outcome, "execution_error", text)

    def test_wellformed_unsolved_is_still_unresolved(self):
        self.assertEqual(classify_planner_json({"smiles": "CCO", "solvable": False, "trees": []}).outcome, "unresolved")


class Resources(Tmp):
    def test_F13_concurrent_registration_keeps_every_entry(self):
        root = self.p("store")
        files = []
        for i in range(8):
            files.append(self.write(f"c{i}.smi", "CCO\nCCN\n"))
        ctx = multiprocessing.get_context("fork")
        procs = [ctx.Process(target=_register, args=(root, f"res{i}", files[i])) for i in range(8)]
        for q in procs:
            q.start()
        for q in procs:
            q.join()
        self.assertEqual([q.exitcode for q in procs], [0] * 8)
        self.assertEqual(sorted(Store(root).all()), [f"res{i}" for i in range(8)])

    def test_F13_forced_interleaving_loses_nothing(self):
        """Both writers read the registry, then one is slowed inside the save: without a transaction lock the second save
        overwrites the first (the review's forced schedule)."""
        import threading, time
        root = self.p("store")
        f1, f2 = self.write("a.smi", "CCO\n"), self.write("b.smi", "CCN\n")
        real_save = Store._save
        def slow_save(self_, reg):
            time.sleep(0.3)
            real_save(self_, reg)
        errors = []
        def reg(rid, f):
            try:
                Store(root).register(rid, f, fmt="smiles", hash_it=False)
            except Exception as e:             # noqa: BLE001
                errors.append(e)
        with mock.patch.object(Store, "_save", slow_save):
            ts = [threading.Thread(target=reg, args=(r, f)) for r, f in (("one", f1), ("two", f2))]
            for t in ts:
                t.start()
            for t in ts:
                t.join()
        self.assertEqual(errors, [])
        self.assertEqual(sorted(Store(root).all()), ["one", "two"])

    def test_F07_unusable_catalogues_are_invalid(self):
        os.makedirs(self.p("dir"))
        self.write("dir/README.txt", "This directory contains the catalogue.\nSee the website.\n")
        st = Store(self.p("s"))
        self.assertEqual(st.register("enamine_real", self.p("dir"), fmt="smiles_or_sdf_directory")["validation"]["status"], "invalid")
        fake_pq = self.write("fake.parquet", "PAR1" + "x" * 7 + "PAR1")
        fake_h5 = self.write("fake.hdf5", "\x89HDF\r\n\x1a\n" + "0" * 7)
        self.assertNotEqual(st.register("a", fake_pq, fmt="parquet")["validation"]["status"], "ok")
        self.assertNotEqual(st.register("b", fake_h5, fmt="aizynthfinder_hdf5")["validation"]["status"], "ok")
        # and the policy refuses to treat an invalid registered catalogue as usable
        cfg = load_config(None, {"resources": ["enamine_real"], "pillars": ["P3"]})
        r = resolve(cfg, load_manifest(), st, mode="replay")
        self.assertEqual(r.states["enamine_real"].state, "invalid")
        self.assertEqual(r.exit_code, 2)

    def test_F08_catalogue_reader_keeps_a_last_record_without_terminator(self):
        sdf = self.p("bb.sdf")
        write_sdf(sdf, [("a", "CCO", {}), ("b", "CCN", {})])
        txt = open(sdf).read()
        self.assertTrue(txt.rstrip().endswith("$$$$"))
        with open(sdf, "w") as fh:
            fh.write(txt.rstrip()[:-4].rstrip("\n") + "\n")       # drop the final terminator
        self.assertEqual(len(list(read_sdf(sdf))), 2)
        self.assertEqual(load_building_blocks(sdf), {"CCO", "CCN"})
        self.assertEqual(Store(self.p("s")).register("enamine_bb", sdf, fmt="sdf")["validation"]["n_records"], 2)

    def test_F06_changed_registered_file_is_not_trusted(self):
        f = self.write("bb.smi", "CCO\nCCN\n")
        st = Store(self.p("s"))
        st.register("enamine_bb", f, fmt="smiles")
        cfg = load_config(None, {"resources": ["enamine_bb"], "pillars": ["P3"]})
        self.assertEqual(resolve(cfg, load_manifest(), st, mode="replay").states["enamine_bb"].state, "ready")
        with open(f, "w") as fh:
            fh.write("CCCC\nCCCO\n")                            # same size, different content
        os.utime(f, (1, 1))
        r = resolve(cfg, load_manifest(), st, mode="replay")
        self.assertEqual(r.states["enamine_bb"].state, "invalid")
        self.assertIn("changed since import", r.states["enamine_bb"].reason)
        with open(f, "w") as fh:
            fh.write("this is not a catalogue\n" * 10)
        self.assertEqual(resolve(cfg, load_manifest(), st, mode="replay").states["enamine_bb"].state, "invalid")

    def test_F06_plan_fingerprint_sees_content_of_same_named_stocks(self):
        from synspeq.cli import _check_plan_fingerprint
        os.makedirs(self.p("m"))
        self.write("m/config.yml", "expansion: {}\n")
        self.write("m/uspto_model.onnx", "model-bytes")
        class St:
            def __init__(s, path): s.path, s.state = path, "ready"
        class Res:
            states = {"aizynth_models": St(self.p("m"))}
        class Store_:
            root = self.p("store")
        out = self.p("out"); os.makedirs(out)
        a = self.write("a/stock.hdf5", "AAAA"); b = self.write("b/stock.hdf5", "BBBB")
        self.assertEqual(_check_plan_fingerprint(out, Res, Store_, {"s": a}), 0)
        self.assertEqual(_check_plan_fingerprint(out, Res, Store_, {"s": a}), 0)
        self.assertEqual(_check_plan_fingerprint(out, Res, Store_, {"s": b}), 2)          # same basename, other content
        self.write("m/uspto_model.onnx", "other-model-bytes")                              # same size not needed: content differs
        self.assertEqual(_check_plan_fingerprint(out, Res, Store_, {"s": a}), 2)

    def test_F16_partial_planner_stock_blocks_planning(self):
        f = self.write("zinc.smi", "CCO\nCCN\n")
        st = Store(self.p("s"))
        st.register("aizynth_zinc", f, fmt="smiles", coverage={"fraction": 0.5})
        md = self.p("models"); os.makedirs(md)
        self.write("models/config.yml", "x: 1\n"); self.write("models/m.onnx", "x")
        st.register("aizynth_models", md, check={"dir_requires": ["config.yml"], "dir_requires_ext": ".onnx"})
        cfg = load_config(None, {"pillars": ["P1"]})
        r = resolve(cfg, load_manifest(), st, mode="plan")
        self.assertEqual(r.exit_code, 2)
        self.assertEqual(r.metrics["P1"].status, "blocked")

    def test_enamine_bb_is_a_purchasing_catalogue_unless_listed_as_planner_stock(self):
        f = self.write("bb.smi", "CCO\nCCN\n")
        st = Store(self.p("s"))
        st.register("enamine_bb", f, fmt="smiles")
        man = load_manifest()
        r = resolve(load_config(None, {"resources": ["enamine_bb"], "pillars": ["P1", "P3"]}), man, st, mode="plan")
        self.assertEqual(r.planner_stocks, [])                       # the default stock is still required; it is not installed here
        self.assertEqual(r.states["aizynth_zinc"].requested_by, "default")
        r = resolve(load_config(None, {"planner_stock": ["enamine_bb"], "pillars": ["P1"]}), man, st, mode="plan")
        self.assertEqual(r.planner_stocks, ["enamine_bb"])
        self.assertEqual(r.states["aizynth_zinc"].state, "not_requested")

    def test_F03_smiles_catalogue_becomes_an_inchikey_stock(self):
        from synspeq.resources.stocks import compiled_planner_stock
        f = self.write("bb.smi", "CCO\nCCN\nnot_a_smiles\n")
        out = compiled_planner_stock("enamine_bb", f, self.p("cache"))
        keys = {Chem.MolToInchiKey(Chem.MolFromSmiles(s)) for s in ("CCO", "CCN")}
        if out.endswith(".txt"):
            got = {l.strip() for l in open(out) if l.strip()}
        else:
            import pandas as pd
            got = set(pd.read_hdf(out, key="table")["inchi_key"])
        self.assertEqual(got, keys)
        try:
            from aizynthfinder.context.stock.queries import InMemoryInchiKeyQuery
            from aizynthfinder.chem import Molecule
        except Exception:
            return                                                    # AiZynthFinder is an optional dependency
        q = InMemoryInchiKeyQuery(out)
        self.assertIn(Molecule(smiles="CCO"), q)
        self.assertIn(Molecule(smiles="CCN"), q)
        self.assertNotIn(Molecule(smiles="CCC"), q)

    def test_F18_config_data_dir_selects_the_store(self):
        bb = self.write("bb.smi", "CCO\n")
        Store(self.p("good")).register("enamine_bb", bb, fmt="smiles")
        cfg = self.write("c.yaml", f"data_dir: {self.p('good')}\nresources: [enamine_bb]\npillars: [P1, P3]\n")
        sdf = self.p("m.sdf"); write_sdf(sdf, [("a", "CCO", {})])
        os.makedirs(self.p("r")); json.dump(direct("CCO", "a"), open(self.p("r", "a.json"), "w"))
        old = os.environ.get("SYNSPEQ_DATA_DIR")
        os.environ["SYNSPEQ_DATA_DIR"] = self.p("empty")
        try:
            import io, contextlib
            with contextlib.redirect_stdout(io.StringIO()):
                rc = main(["evaluate", sdf, "--routes-dir", self.p("r"), "--config", cfg, "--out", self.p("o")])
        finally:
            if old is None:
                os.environ.pop("SYNSPEQ_DATA_DIR")
            else:
                os.environ["SYNSPEQ_DATA_DIR"] = old
        self.assertEqual(rc, 0)

    def test_config_rejects_inactive_and_mistyped_options(self):
        for bad in ({"auto_download": True}, {"on_mising": "skip"}, {"resources": "enamine_bb"}, {"paths": {"a": 1}}):
            with self.assertRaises(ValueError):
                load_config(None, bad)


def _register(root, rid, path):
    Store(root).register(rid, path, fmt="smiles", hash_it=False)


class RealScan(Tmp):
    def _catalogue(self):
        self.write("REAL/part0.smi", "CCO\nCCC\n")
        self.write("REAL/part1.smi", "CCN\nCCCC\n")
        return self.p("REAL")

    def test_F02_partial_scan_keeps_unqueried_shards_unknown(self):
        real = self._catalogue()
        q = self.write("q.smi", "CCN\nOCCO\n")                     # CCN is in part1 only; OCCO is in neither
        scan_shard(q, real, 0, 2, self.p("h0.txt"))
        s = merge_shards(q, [self.p("h0.txt")], self.p("lk"))
        self.assertFalse(s["complete"])
        lk = RealLookup.load(self.p("lk"))
        self.assertEqual(lk.state("CCN"), "unknown")                # not a measured non-hit
        self.assertEqual(lk.state("OCCO"), "unknown")
        scan_shard(q, real, 1, 2, self.p("h1.txt"))
        s = merge_shards(q, [self.p("h0.txt"), self.p("h1.txt")], self.p("lk2"))
        self.assertTrue(s["complete"])
        lk = RealLookup.load(self.p("lk2"))
        self.assertEqual((lk.state("CCN"), lk.state("OCCO")), ("hit", "non_hit"))

    def test_hit_without_header_or_for_other_queries_is_refused(self):
        q = self.write("q.smi", "CCN\n")
        self.write("old.txt", "CCN\tpart1.smi\n")
        with self.assertRaises(ValueError):
            merge_shards(q, [self.p("old.txt")], self.p("lk"))
        real = self._catalogue()
        scan_shard(q, real, 0, 1, self.p("h.txt"))
        other = self.write("q2.smi", "CCC\n")
        with self.assertRaises(ValueError):
            merge_shards(other, [self.p("h.txt")], self.p("lk"))

    def test_changed_catalogue_is_not_mixed_into_one_lookup(self):
        real = self._catalogue()
        q = self.write("q.smi", "CCN\n")
        scan_shard(q, real, 0, 1, self.p("h.txt"))
        merge_shards(q, [self.p("h.txt")], self.p("lk"))
        self.write("REAL/part2.smi", "CCCCC\n")
        scan_shard(q, real, 0, 1, self.p("h2.txt"))
        with self.assertRaises(ValueError):
            merge_shards(q, [self.p("h2.txt")], self.p("lk"))


class Metrics(Tmp):
    def test_F11_unmeasured_on_demand_share_is_null(self):
        sdf = self.p("m.sdf"); write_sdf(sdf, [("a", "CCCCO", {})])
        t = mol("CCCCO", [rxn(mol("CCCO", stock=False, children=[rxn(mol("CCO", stock=True), mol("CO", stock=True))]), mol("CC", stock=True))])
        pop = make_population(read_sdf(sdf), "t")
        o = {pop.molecules[0].mol_id: classify_planner_json({"smiles": "CCCCO", "solvable": True, "n_steps": 2, "trees": [t]})}
        r = compute_p3(pop, o, Catalogue({"CCCO"}, None))
        self.assertIsNone(r["needs_on_demand"])
        self.assertEqual(r["needs_on_demand_status"], "not_assessed")
        self.assertEqual(r["median_residual_steps"], 1)
        r = compute_p3(pop, o, Catalogue({"CCCO"}, RealLookup()))
        self.assertEqual(r["needs_on_demand_status"], "bound")       # nothing was looked up in REAL yet

    def test_F12_prior_and_query_use_the_same_scaffold_policy(self):
        smiles = ["c1ccccc1"] * 6 + ["CCO"] * 4
        prior = scaffold_prior(smiles)
        self.assertNotIn("", prior)
        from synspeq.coverage import generic_scaffold
        self.assertAlmostEqual(cosine_to_prior([generic_scaffold(s) for s in smiles], prior), 1.0)

    def test_F20_undefined_distributions_are_null(self):
        self.assertEqual(n50_csr([]), (None, 0, None, 0))
        self.assertEqual(n50_csr(["", "", None]), (None, 0, None, 0))
        self.assertEqual(n50_csr(["A"])[0], 1)                       # a single cyclic molecule is defined
        self.assertIsNone(subsample_n50([""] * 10, n=5))
        sdf = self.p("m.sdf"); write_sdf(sdf, [("a", "CCO", {}), ("b", "CCN", {})])
        pop = make_population(read_sdf(sdf), "t")
        outs = {m.mol_id: classify_planner_json(direct(m.identity)) for m in pop.molecules}
        b = compute_p4(pop, outs)["all_distinct"]
        self.assertIsNone(b["n50"])
        self.assertIsNone(b["csr_auc"])
        self.assertEqual(b["n_acyclic"], 2)


class Caches(Tmp):
    def test_F17_interrupted_append_is_recovered_interior_corruption_is_not(self):
        from synspeq.support.engine import FeatureCache
        good = json.dumps({"rxn": "A>>B", "_db": "k", "prep_status": "ok"}) + "\n"
        f = self.write("f.jsonl", good + '{"rxn": ')
        c = FeatureCache(f, "k")
        self.assertIn("A>>B", c.rows)
        self.assertTrue(c.recovered_truncated_tail)
        self.assertEqual(open(f).read(), good)                      # truncated back to the last complete row
        c.get_or_compute("C>>D", lambda r: {"rxn": r, "prep_status": "ok"})
        self.assertEqual(set(FeatureCache(f, "k").rows), {"A>>B", "C>>D"})
        bad = self.write("g.jsonl", good + "{broken}\n" + good)
        with self.assertRaises(ValueError):
            FeatureCache(bad, "k")
        # rows made under another key are ignored, not reused
        self.assertEqual(FeatureCache(f, "other").rows, {})

    def test_F05_features_are_not_reused_across_databases(self):
        import pandas as pd
        out, cache = self.p("o"), self.p("c")
        db_a, db_b = self.p("a", "db.parquet"), self.p("b", "db.parquet")
        for path, rx in ((db_a, "CC>>CCO"), (db_b, "CC>>CCN")):
            os.makedirs(os.path.dirname(path))
            pd.DataFrame({"REACTION": [rx], "REF": [""], "TAG": ["t"], "TAG2": ["t2"], "rxn_str_patt_fp": ["0" * 334]}).to_parquet(path)
        calls = []
        def fake_features(rxn, prep, searcher):
            calls.append(searcher.index.reaction[0])
            return {"rxn": rxn, "prep_status": "ok", "v": 1 if searcher.index.reaction[0] == "CC>>CCO" else 0}
        from synspeq import pipeline
        with mock.patch("synspeq.support.engine.Preprocessor", lambda: object()), mock.patch("synspeq.support.engine.reaction_features", fake_features):
            ra = pipeline._p2_features(["X>>Y"], db_a, out, cache)
            rb = pipeline._p2_features(["X>>Y"], db_b, out, cache)
            ra2 = pipeline._p2_features(["X>>Y"], db_a, out, cache)
        self.assertEqual((ra["X>>Y"]["v"], rb["X>>Y"]["v"], ra2["X>>Y"]["v"]), (1, 0, 1))
        self.assertEqual(calls, ["CC>>CCO", "CC>>CCN"])             # A's row is reused for A, never for B

    def test_F04_identity_cache_belongs_to_one_database(self):
        import pandas as pd
        from synspeq.support.index import PrecedentIndex, key_hash, rxn_key
        cache = self.p("ids.parquet")
        db = {}
        for n, rx in (("a", "CC>>CCO"), ("b", "CC>>CCN")):
            db[n] = self.p(n + ".parquet")
            pd.DataFrame({"REACTION": [rx], "REF": [""], "TAG": ["t"], "TAG2": ["t2"], "rxn_str_patt_fp": ["1" * 334]}).to_parquet(db[n])
        ia = PrecedentIndex.from_parquet(db["a"], identity_cache=cache, cache_meta={"db_sha256": "A"})
        ib = PrecedentIndex.from_parquet(db["b"], identity_cache=cache, cache_meta={"db_sha256": "B"})
        self.assertEqual(ia.key_stereo[0], key_hash(rxn_key("CC>>CCO")))
        self.assertEqual(ib.key_stereo[0], key_hash(rxn_key("CC>>CCN")))      # not A's key
        ib2 = PrecedentIndex.from_parquet(db["b"], identity_cache=cache, cache_meta={"db_sha256": "B"})
        self.assertEqual(ib2.key_stereo[0], ib.key_stereo[0])                 # a matching cache is reused
        # without provenance a cache is never reused
        ic = PrecedentIndex.from_parquet(db["a"], identity_cache=cache)
        self.assertEqual(ic.key_stereo[0], key_hash(rxn_key("CC>>CCO")))


class Interfaces(Tmp):
    def test_F19_invalid_counts_are_rejected(self):
        sdf = self.p("m.sdf"); write_sdf(sdf, [(str(i), "C" * (i + 1), {}) for i in range(4)])
        for kw in ({"last": 0}, {"first": -1}, {"last": -1}):
            with self.assertRaises(ValueError):
                make_population(read_sdf(sdf), "t", **kw)
        from synspeq.planners.aizynth import plan_population
        for kw in ({"workers": 0}, {"workers": -1}, {"n_chunks": 0}, {"chunk": 2, "n_chunks": 2}, {"chunk": -1}):
            with self.assertRaises(ValueError):
                plan_population([("a", "CCO")], "cfg.yml", self.p("o"), **kw)
        with self.assertRaises(SystemExit):
            build_parser().parse_args(["validate", sdf, "--first", "0"])
        with self.assertRaises(SystemExit):
            build_parser().parse_args(["plan", sdf, "--out", "o", "--workers", "0"])

    def test_F22_bad_gzip_is_a_message_not_a_traceback(self):
        f = self.write("bad.sdf.gz", "this is not gzip")
        import contextlib, io
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rc = main(["validate", f])
        self.assertEqual(rc, 2)
        self.assertIn("synspeq: error", err.getvalue())

    def test_F21_every_documented_command_parses(self):
        """Commands printed in README.md and docs/DATA_RESOURCES_GUIDE.md must be accepted by the real argument parser."""
        parser = build_parser()
        n = 0
        for doc in ("README.md", os.path.join("docs", "DATA_RESOURCES_GUIDE.md")):
            text = open(os.path.join(ROOT, doc)).read()
            for block in re.findall(r"```(?:bash|sh)?\n(.*?)```", text, flags=re.S):
                joined = re.sub(r"\\\n\s*", " ", block)
                for line in joined.splitlines():
                    line = re.sub(r"\s#.*$", "", line).strip()
                    if not line.startswith("synspeq "):
                        continue
                    argv = [a for a in shlex.split(line)[1:] if a != "#"]
                    argv = [re.sub(r"^<.*>$", "x", a) for a in argv]
                    try:
                        parser.parse_args(argv)
                    except SystemExit:
                        self.fail(f"{doc}: documented command is not accepted: {line}")
                    n += 1
        self.assertGreater(n, 10)


if __name__ == "__main__":
    unittest.main()
