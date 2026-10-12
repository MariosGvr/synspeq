"""The repository must not contain machine-specific absolute paths (home directories, cluster mounts)."""
import os, re, subprocess, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PATTERN = re.compile(r"(/net/[a-z]+/|/home/[A-Za-z0-9_.-]+/|/Users/[A-Za-z0-9_.-]+/|C:\\\\Users|/scratch/|/mnt/[a-z])")
SKIP = ("LICENSE", ".svg", ".png", ".gz", ".pkl")


class Hygiene(unittest.TestCase):
    def test_no_absolute_paths_in_tracked_text_files(self):
        try:
            files = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
        except Exception:
            self.skipTest("not a git checkout")
        bad = []
        for f in files:
            if f.endswith(SKIP):
                continue
            try:
                txt = open(os.path.join(ROOT, f), encoding="utf-8").read()
            except (UnicodeDecodeError, OSError):
                continue
            for i, line in enumerate(txt.splitlines(), 1):
                if PATTERN.search(line) and "test_hygiene" not in f:
                    bad.append(f"{f}:{i}: {line.strip()[:100]}")
        self.assertEqual(bad, [], "absolute paths found:\n" + "\n".join(bad))


if __name__ == "__main__":
    unittest.main()
