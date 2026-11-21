"""Behaviour a new user meets first: names in the output, where the resource store lives, what the report opens with."""
import contextlib, csv, io, json, os, re, shutil, tempfile, unittest
from unittest import mock
from synspeq.cli import main
from synspeq.environment import environment_info, rdkit_note
from synspeq.resources import Store, load_manifest
from synspeq.resources.store import default_dir

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEMO = os.path.join(ROOT, "examples", "demo")


def run_cli(*argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        rc = main(list(argv))
    return rc, out.getvalue(), err.getvalue()


class Tmp(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)

    def p(self, *a):
        return os.path.join(self.d, *a)


class Names(Tmp):
    """The output uses the paper's names, never the internal identifiers P1-P4."""
    def evaluate_demo(self):
        store = self.p("store")
        self.assertEqual(run_cli("--data-dir", store, "data", "import", "enamine_bb", os.path.join(DEMO, "building_blocks.smi"))[0], 0)
        rc, out, err = run_cli("--data-dir", store, "evaluate", os.path.join(DEMO, "molecules.sdf"), "--routes-dir", os.path.join(DEMO, "routes"),
                               "--pocket-field", "pocket_id", "--resources", "enamine_bb", "--measures", "findability,residual-steps,scaffold-diversity",
                               "--out", self.p("out"))
        self.assertEqual(rc, 0, err)
        return out

    def test_report_files_use_measure_names(self):
        self.evaluate_demo()
        rep = json.load(open(self.p("out", "report.json")))
        for key in ("findability", "residual_steps", "scaffold_diversity"):
            self.assertIn(key, rep)
            self.assertIn(key, rep["status"])
        self.assertFalse({"P1", "P2", "P3", "P4"} & (set(rep) | set(rep["status"])))
        with open(self.p("out", "molecules.csv")) as fh:
            cols = next(csv.reader(fh))
        self.assertFalse([c for c in cols if re.match(r"p[1-4]_", c)], cols)
        self.assertTrue(any(c.startswith("residual_steps_") for c in cols))
        html = open(self.p("out", "report.html"), encoding="utf-8").read()
        self.assertFalse(re.search(r"\bP[1-4]\b", html))

    def test_report_opens_with_a_summary_of_the_four_measures(self):
        self.evaluate_demo()
        html = open(self.p("out", "report.html"), encoding="utf-8").read()
        self.assertIn("<h2>Summary</h2>", html)
        self.assertLess(html.index("<h2>Summary</h2>"), html.index("<h2>Input</h2>"))
        for name in ("Findability", "Precedent support", "Residual steps", "Scaffold diversity"):
            self.assertIn(name, html.split("<h2>Input</h2>")[0])

    def test_data_check_names_the_measures_and_shows_the_store(self):
        store = self.p("store")
        rc, out, err = run_cli("--data-dir", store, "data", "check")
        text = out + err
        self.assertEqual(rc, 2)
        self.assertIn("resource store:", out)
        self.assertFalse(re.search(r"\bP[1-4]\b", text), text)
        self.assertIn("findability:", out)

    def test_a_missing_resource_is_reported_once(self):
        store = self.p("store")
        rc, out, err = run_cli("--data-dir", store, "run", os.path.join(DEMO, "molecules.sdf"), "--out", self.p("out"))
        self.assertEqual(rc, 2)
        lines = [l for l in err.splitlines() if l.startswith("[error]") and "aizynth_models" in l]
        self.assertEqual(len(lines), 1, err)
        self.assertIn("resource store", err)


class StoreLocation(Tmp):
    def test_environment_variable_wins_then_a_local_folder_then_the_user_folder(self):
        cwd = os.getcwd()
        os.chdir(self.d)
        self.addCleanup(os.chdir, cwd)
        env = {k: v for k, v in os.environ.items() if k != "SYNSPEQ_DATA_DIR"}
        env["XDG_DATA_HOME"] = self.p("xdg")
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertEqual(default_dir(), self.p("xdg", "synspeq"))        # not tied to the current directory
            os.mkdir(self.p(".synspeq"))
            self.assertEqual(os.path.realpath(default_dir()), os.path.realpath(self.p(".synspeq")))
        with mock.patch.dict(os.environ, {**env, "SYNSPEQ_DATA_DIR": self.p("elsewhere")}, clear=True):
            self.assertEqual(default_dir(), self.p("elsewhere"))

    def test_extra_files_next_to_the_models_do_not_invalidate_them(self):
        models = self.p("models")
        os.makedirs(models)
        for name in ("config.yml", "uspto_model.onnx"):
            with open(os.path.join(models, name), "w") as fh:
                fh.write("x")
        man = load_manifest()
        store = Store(self.p("store"))
        entry = store.register("aizynth_models", models, fmt=man.resources["aizynth_models"].format, check=man.resources["aizynth_models"].check)
        self.assertEqual(entry["validation"]["status"], "ok")
        self.assertIsNone(Store.changed_since_import(store.get("aizynth_models")))
        with open(os.path.join(models, "zinc_stock.hdf5"), "w") as fh:        # the stock is downloaded into the same folder
            fh.write("stock")
        with open(os.path.join(models, "download.log"), "w") as fh:
            fh.write("log")
        self.assertIsNone(Store.changed_since_import(store.get("aizynth_models")))
        with open(os.path.join(models, "uspto_model.onnx"), "w") as fh:        # a model file really changed
            fh.write("changed model")
        self.assertIn("changed", Store.changed_since_import(store.get("aizynth_models")))


    def test_a_store_made_by_the_earlier_version_stays_valid(self):
        models = self.p("models")
        os.makedirs(models)
        for name, text in (("config.yml", "x"), ("uspto_model.onnx", "x"), ("zinc_stock.hdf5", "stock")):
            with open(os.path.join(models, name), "w") as fh:
                fh.write(text)
        man = load_manifest()
        store = Store(self.p("store"))
        store.register("aizynth_models", models, fmt=None, check=man.resources["aizynth_models"].check)
        entry = store.get("aizynth_models")
        every = [f for f in sorted(os.listdir(models))]                  # the earlier version described every file of the folder
        from synspeq.hashing import listing_id
        entry.update(listing_sha256=listing_id(models, every), n_files=len(every),
                     size_bytes=sum(os.path.getsize(os.path.join(models, f)) for f in every),
                     mtime_ns=max(os.stat(os.path.join(models, f)).st_mtime_ns for f in every))
        self.assertIsNone(Store.changed_since_import(entry))
        with open(os.path.join(models, "uspto_model.onnx"), "w") as fh:
            fh.write("changed model")
        self.assertIn("changed", Store.changed_since_import(entry))


class EnvironmentRecord(Tmp):
    def test_versions_are_recorded_in_the_report_and_the_manifest(self):
        store = self.p("store")
        run_cli("--data-dir", store, "data", "import", "enamine_bb", os.path.join(DEMO, "building_blocks.smi"))
        run_cli("--data-dir", store, "evaluate", os.path.join(DEMO, "molecules.sdf"), "--routes-dir", os.path.join(DEMO, "routes"),
                "--pocket-field", "pocket_id", "--measures", "findability,scaffold-diversity", "--out", self.p("out"))
        rep = json.load(open(self.p("out", "report.json")))
        man = json.load(open(self.p("out", "run_manifest.json")))
        for d in (rep["environment"], man["environment"]):
            self.assertIn("rdkit", d)
            self.assertIn("synspeq", d)
        self.assertEqual(environment_info()["rdkit"], rep["environment"]["rdkit"])

    def test_note_only_when_rdkit_is_not_the_validated_release(self):
        with mock.patch("rdkit.__version__", "2023.09.6"):
            self.assertIsNone(rdkit_note())
        with mock.patch("rdkit.__version__", "2025.03.1"):
            self.assertIn("2025.03.1", rdkit_note())

    def test_no_stray_lock_file_in_the_results_folder(self):
        from synspeq.support.engine import FeatureCache
        cache = FeatureCache(self.p("out", "precedent_support_features.jsonl"), "key")
        cache.get_or_compute("A>>B", lambda r: {"x": 1})
        self.assertEqual([f for f in os.listdir(self.p("out")) if f.endswith(".lock") and not f.startswith(".")], [])


if __name__ == "__main__":
    unittest.main()
