import contextlib, io, json, os, tempfile, unittest
from synspeq.cli import main


class CLI(unittest.TestCase):
    def run_cli(self, *argv):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = main(list(argv))
        return rc, buf.getvalue()

    def test_validate_and_p1_replay(self):
        d = tempfile.mkdtemp()
        with open(os.path.join(d, "m.csv"), "w") as fh:
            fh.write("mol_id,smiles,pocket_id\na,CCO,P\nb,OCC,P\nc,CCN,P\n")
        os.makedirs(os.path.join(d, "r"))
        with open(os.path.join(d, "r", "a.json"), "w") as fh:
            json.dump({"smiles": "CCO", "solvable": True, "n_steps": 0, "trees": [{"type": "mol", "smiles": "CCO", "in_stock": True}]}, fh)
        rc, out = self.run_cli("validate", os.path.join(d, "m.csv"), "--pocket-field", "pocket_id")
        self.assertEqual(rc, 0)
        self.assertEqual(json.loads(out)["n_distinct_molecules"], 2)
        rc, out = self.run_cli("p1-replay", os.path.join(d, "m.csv"), "--pocket-field", "pocket_id", "--routes-dir", os.path.join(d, "r"))
        self.assertEqual(rc, 0)
        self.assertIn("direct 1", out)
        self.assertIn("pending 1", out)             # record c has no planner record yet: pending, not unresolved
        self.assertIn("complete=False", out)

    def test_evaluate_p1_p3_p4_with_missing_real_not_assessed_parts(self):
        d = tempfile.mkdtemp()
        with open(os.path.join(d, "m.csv"), "w") as fh:
            fh.write("mol_id,smiles,pocket_id\na,CCCCO,P\nb,c1ccccc1CCO,P\nc,CCN,P\n")
        os.makedirs(os.path.join(d, "r"))
        t = {"type": "mol", "smiles": "CCCCO", "in_stock": False, "children": [{"type": "reaction", "children": [
            {"type": "mol", "smiles": "CCCO", "in_stock": False, "children": [{"type": "reaction", "children": [
                {"type": "mol", "smiles": "CCO", "in_stock": True}, {"type": "mol", "smiles": "CO", "in_stock": True}]}]},
            {"type": "mol", "smiles": "CC", "in_stock": True}]}]}
        with open(os.path.join(d, "r", "a.json"), "w") as fh:
            json.dump({"smiles": "CCCCO", "solvable": True, "n_steps": 2, "trees": [t]}, fh)
        with open(os.path.join(d, "r", "b.json"), "w") as fh:
            json.dump({"smiles": "c1ccccc1CCO", "solvable": False, "trees": []}, fh)
        with open(os.path.join(d, "r", "c.json"), "w") as fh:
            json.dump({"smiles": "CCN", "solvable": True, "n_steps": 0, "trees": [{"type": "mol", "smiles": "CCN", "in_stock": True}]}, fh)
        bb = os.path.join(d, "bb.smi")
        with open(bb, "w") as fh:
            fh.write("CCCO\n")
        sd = os.path.join(d, "store")
        self.assertEqual(self.run_cli("--data-dir", sd, "data", "import", "enamine_bb", bb)[0], 0)
        out = os.path.join(d, "out")
        rc, text = self.run_cli("--data-dir", sd, "evaluate", os.path.join(d, "m.csv"), "--pocket-field", "pocket_id", "--routes-dir", os.path.join(d, "r"),
                                "--pillars", "P1,P3,P4", "--resources", "enamine_bb", "--out", out)
        self.assertEqual(rc, 0, text)
        rep = json.load(open(os.path.join(out, "report.json")))
        self.assertEqual(rep["P1"]["counts"]["routed"], 1)
        self.assertEqual(rep["P3"]["median_residual_steps"], 1)         # frontier at the stocked intermediate CCCO
        self.assertEqual(rep["status"]["P3"]["sub"]["real_needed"]["status"], "not_assessed")
        self.assertEqual(rep["status"]["P2"]["status"], "not_assessed")
        self.assertIn("P4", rep)
        self.assertTrue(os.path.exists(os.path.join(out, "molecules.csv")))
        html = open(os.path.join(out, "report.html")).read()
        self.assertIn("Pillar 1: route recovery", html)
        self.assertIn("Pillar 3: deliverability", html)
        self.assertIn("P2: not_assessed", html)                  # the banner names what was not assessed
        self.assertNotIn("Pillar 2: chemical support", html)

    def test_fetch_verifies_checksum_and_registers(self):
        import hashlib
        from types import SimpleNamespace
        from synspeq.cli import _fetch
        from synspeq.resources import Store
        d = tempfile.mkdtemp()
        src = os.path.join(d, "src.smi")
        with open(src, "w") as fh:
            fh.write("CCO\nCCN\n")
        md5 = hashlib.md5(open(src, "rb").read()).hexdigest()
        r = SimpleNamespace(id="molport", format="smiles", fetch={"url": "file://" + src, "filename": "got.smi", "md5": md5})
        st = Store(os.path.join(d, "s"))
        self.assertEqual(_fetch(r, os.path.join(d, "dl"), st), 0)
        self.assertEqual(st.get("molport")["validation"]["status"], "ok")
        r.fetch["md5"] = "0" * 32
        self.assertEqual(_fetch(r, os.path.join(d, "dl2"), st), 2)           # mismatch: file removed, nothing registered
        self.assertFalse(os.path.exists(os.path.join(d, "dl2", "got.smi")))

    def test_scaffold_prior_and_real_scan_roundtrip(self):
        d = tempfile.mkdtemp()
        os.makedirs(os.path.join(d, "REAL"))
        for i in range(3):
            with open(os.path.join(d, "REAL", f"p{i}.smi"), "w") as fh:
                fh.write("c1ccccc1CC\nC1CCNCC1\nCCO\n")
        out = os.path.join(d, "prior.tsv")
        rc, _ = self.run_cli("scaffold-prior", "--real-dir", os.path.join(d, "REAL"), "--out", out, "--seed", "1", "--n-parts", "3", "--rows-per-part", "3")
        self.assertEqual(rc, 0)
        self.assertTrue(open(out).read().strip())
        q = os.path.join(d, "q.smi")
        with open(q, "w") as fh:
            fh.write("C1CCNCC1\nOCCO\n")
        h = os.path.join(d, "h0.txt")
        self.run_cli("real-scan", "--queries", q, "--real-dir", os.path.join(d, "REAL"), "--out", h)
        lk = os.path.join(d, "lk")
        self.run_cli("real-merge", "--queries", q, "--hits", h, "--out", lk)
        from synspeq.purchasing import RealLookup
        r = RealLookup.load(lk)
        self.assertEqual((r.state("C1CCNCC1"), r.state("OCCO"), r.state("CCCCCC")), ("hit", "non_hit", "unknown"))

    def test_aizynth_download_registers_models_and_stock(self):
        import stat
        d = tempfile.mkdtemp()
        shim = os.path.join(d, "bin"); os.makedirs(shim)
        script = os.path.join(shim, "download_public_data")
        with open(script, "w") as fh:                    # stand-in for AiZynthFinder's downloader
            fh.write("#!/bin/sh\nmkdir -p \"$1\"\nprintf 'expansion: {}\\n' > \"$1/config.yml\"\n"
                     "printf '\\211HDF\\r\\n\\032\\n0000' > \"$1/zinc_stock.hdf5\"\n")
        os.chmod(script, os.stat(script).st_mode | stat.S_IEXEC)
        old = os.environ.get("PATH", "")
        os.environ["PATH"] = shim + os.pathsep + old
        try:
            rc, out = self.run_cli("--data-dir", os.path.join(d, "s"), "data", "download", "aizynth_zinc", "--dest", os.path.join(d, "dl"), "--accept-terms")
        finally:
            os.environ["PATH"] = old
        self.assertEqual(rc, 0, out)
        from synspeq.resources import Store
        st = Store(os.path.join(d, "s"))
        self.assertEqual(st.get("aizynth_models")["validation"]["status"], "ok")
        self.assertEqual(st.get("aizynth_zinc")["validation"]["status"], "ok")

    def test_plan_records_are_not_mixed_across_stocks(self):
        from types import SimpleNamespace
        from synspeq.cli import _check_plan_fingerprint
        d = tempfile.mkdtemp()
        st = lambda p: SimpleNamespace(state="ready", path=p)
        res_a = SimpleNamespace(states={"aizynth_zinc": st("/x/zinc_stock.hdf5"), "aizynth_models": st("/x")})
        res_b = SimpleNamespace(states={"aizynth_zinc": st("/y/other_stock.hdf5"), "aizynth_models": st("/x")})
        ent = {"aizynth_zinc": {"sha256": "a" * 64}}
        self.assertEqual(_check_plan_fingerprint(d, res_a, ent), 0)
        self.assertEqual(_check_plan_fingerprint(d, res_a, ent), 0)          # same stock: reuse allowed
        self.assertEqual(_check_plan_fingerprint(d, res_b, {"aizynth_zinc": {"sha256": "b" * 64}}), 2)

    def test_data_list_and_info(self):
        d = tempfile.mkdtemp()
        rc, out = self.run_cli("--data-dir", d, "data", "list")
        self.assertEqual(rc, 0)
        self.assertIn("enamine_real", out)
        rc, out = self.run_cli("--data-dir", d, "data", "info", "enamine_real")
        self.assertIn("Enamine Terms of use", out)

    def test_check_exit_code_when_core_missing(self):
        d = tempfile.mkdtemp()
        rc, out = self.run_cli("--data-dir", d, "data", "check")
        self.assertEqual(rc, 2)
        self.assertIn("synspeq data download aizynth_zinc", out)


if __name__ == "__main__":
    unittest.main()
