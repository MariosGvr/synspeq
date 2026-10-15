"""Provenance, ambiguity and concurrent cache-recovery cases."""
import json
import os
import threading
from unittest import mock
from tests.unit.test_regressions import Tmp, direct
from tests.unit.test_core import write_sdf
from synspeq.io import read_sdf
from synspeq.populations import make_population
from synspeq.pipeline import load_outcomes
from synspeq.purchasing.realscan import scan_shard, merge_shards, check_lookup_matches_catalogue, HEADER
from synspeq.support.engine import FeatureCache


class ProvenanceBoundaries(Tmp):
    def lookup(self):
        self.write('A/part.smi', 'CCO\n')
        self.write('B/part.smi', 'CCN\n')
        q = self.write('q.smi', 'CCO\n')
        scan_shard(q, self.p('A'), 0, 1, self.p('h.txt'))
        merge_shards(q, [self.p('h.txt')], self.p('lk'))
        return q

    def test_content_is_always_verified(self):
        self.lookup()
        check_lookup_matches_catalogue(self.p('lk'), self.p('A'))
        with self.assertRaisesRegex(ValueError, 'differs in content'):
            check_lookup_matches_catalogue(self.p('lk'), self.p('B'))

    def test_verified_files_are_not_hashed_again_but_changes_are_caught(self):
        self.lookup()
        import synspeq.hashing as hashing
        real = hashing.file_sha256
        calls = []
        def counting(path, *a, **k):
            calls.append(path)
            return real(path, *a, **k)
        with mock.patch('synspeq.hashing.file_sha256', counting):
            check_lookup_matches_catalogue(self.p('lk'), self.p('A'))
            first = len(calls)
            check_lookup_matches_catalogue(self.p('lk'), self.p('A'))
            self.assertGreater(first, 0)
            self.assertEqual(len(calls), first, 'an unchanged file was hashed again')
            # same size, new content, new modification time: must be hashed and refused
            f = self.p('A', 'part.smi')
            with open(f, 'w') as fh:
                fh.write('CCN\n')
            os.utime(f, ns=(1, 1))
            with self.assertRaisesRegex(ValueError, 'differs in content'):
                check_lookup_matches_catalogue(self.p('lk'), self.p('A'))

    def test_the_verification_option_no_longer_exists(self):
        import inspect
        from synspeq.pipeline import evaluate
        from synspeq.cli import build_parser
        self.assertNotIn('verify_catalogue', inspect.signature(evaluate).parameters)
        with self.assertRaises(SystemExit):
            build_parser().parse_args(['evaluate', 'x.sdf', '--routes-dir', 'r', '--out', 'out', '--verify-catalogue'])

    def test_cli_rejects_lookup_from_same_size_different_catalogue(self):
        import contextlib
        import io
        from synspeq.cli import main
        from synspeq.resources import Store
        from tests.unit.test_regressions import mol, rxn
        self.lookup()
        write_sdf(self.p('m.sdf'), [('target', 'CCO', {})])
        record = {'smiles':'CCO', 'solvable':True, 'n_steps':1,
                  'trees':[mol('CCO', [rxn(mol('CC', stock=True), mol('O', stock=True))])]}
        self.write('routes/target.json', json.dumps(record))
        store = Store(self.p('store'))
        store.register('enamine_real', self.p('B'), fmt='smiles_or_sdf_directory')
        args = ['--data-dir', self.p('store'), 'evaluate', self.p('m.sdf'),
                '--routes-dir', self.p('routes'), '--resources', 'enamine_real',
                '--pillars', 'P1,P3', '--out', self.p('out')]
        err = io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            rc = main(args + ['--real-lookup', self.p('lk')])
        self.assertEqual(rc, 2)
        self.assertIn('differs in content', err.getvalue())
        self.assertFalse(os.path.exists(self.p('out/report.json')))

    def test_lookup_without_content_manifest_is_rejected(self):
        self.lookup()
        p = self.p('lk', 'meta.json')
        with open(p) as fh:
            meta = json.load(fh)
        meta.pop('files')
        self.write('lk/meta.json', json.dumps(meta))
        with self.assertRaisesRegex(ValueError, 'content provenance'):
            check_lookup_matches_catalogue(self.p('lk'), self.p('A'))

    def test_complete_scan_cannot_omit_file_provenance(self):
        q = self.lookup()
        with open(self.p('h.txt')) as fh:
            lines = fh.readlines()
        meta = json.loads(lines[0][len(HEADER):])
        meta['files'] = {}
        self.write('h_bad.txt', HEADER + json.dumps(meta) + '\n' + ''.join(lines[1:]))
        with self.assertRaises(ValueError):
            merge_shards(q, [self.p('h_bad.txt')], self.p('other'))

    def test_policy_cannot_change_inside_existing_lookup(self):
        q = self.lookup()
        scan_shard(q, self.p('A'), 0, 1, self.p('h2.txt'), recanonicalise=True)
        with self.assertRaisesRegex(ValueError, 'matching policy changed'):
            merge_shards(q, [self.p('h2.txt')], self.p('lk'))

    def test_partial_lookup_validates_only_scanned_content(self):
        self.write('REAL/one.smi', 'CCO\n')
        self.write('REAL/two.smi', 'CCN\n')
        q = self.write('q.smi', 'CCO\nCCN\n')
        scan_shard(q, self.p('REAL'), 0, 2, self.p('h.txt'))
        merge_shards(q, [self.p('h.txt')], self.p('lk'))
        check_lookup_matches_catalogue(self.p('lk'), self.p('REAL'))

    def test_legacy_title_cannot_bind_multiple_pockets_of_same_molecule(self):
        write_sdf(self.p('m.sdf'), [('same', 'CCO', {'pocket':'A'}), ('same', 'CCO', {'pocket':'B'})])
        pop = make_population(read_sdf(self.p('m.sdf')), 'x', pocket_field='pocket')
        self.write('r/same.json', json.dumps(direct('CCO')))
        results = load_outcomes(pop, [self.p('r')])
        self.assertEqual([r.outcome for r in results.values()], ['execution_error'] * 2)
        self.assertTrue(all('ambiguous' in r.error for r in results.values()))

    def test_conflicting_exact_records_do_not_depend_on_directory_order(self):
        write_sdf(self.p('m.sdf'), [('a', 'CCO', {})])
        pop = make_population(read_sdf(self.p('m.sdf')), 'x')
        mid = pop.molecules[0].mol_id
        self.write('a/' + mid + '.json', json.dumps(direct('CCO')))
        self.write('b/' + mid + '.json', json.dumps({'smiles':'CCO', 'solvable':False,'trees':[]}))
        for dirs in ([self.p('a'),self.p('b')], [self.p('b'),self.p('a')]):
            self.assertEqual(load_outcomes(pop, dirs)[mid].outcome, 'execution_error')
        self.write('b/' + mid + '.json', json.dumps(direct('CCO')))
        self.assertEqual(load_outcomes(pop, [self.p('a'),self.p('b')])[mid].outcome, 'direct_purchase')

    def test_cache_loader_takes_the_append_lock_before_reading(self):
        from synspeq.locking import file_lock
        p = self.write('cache.jsonl', '{"rxn":')
        started, acquired = threading.Event(), threading.Event()
        errors = []
        def observed_lock(*args, **kwargs):
            started.set()
            return file_lock(*args, **kwargs)
        def reader():
            try:
                FeatureCache(p, 'k')
                acquired.set()
            except Exception as e:
                errors.append(e)
        with file_lock(p + '.lock'):
            with mock.patch('synspeq.locking.file_lock', observed_lock):
                t = threading.Thread(target=reader)
                t.start()
                self.assertTrue(started.wait(3), 'loader did not request the append lock')
                self.assertFalse(acquired.is_set())
                # A writer completes the record while holding the append lock.
                self.write('cache.jsonl', json.dumps({'rxn':'A>>B','_db':'k'}) + '\n')
        t.join(3)
        self.assertFalse(t.is_alive())
        self.assertEqual(errors, [])
        self.assertIn('A>>B', FeatureCache(p, 'k').rows)

    def test_reading_a_missing_cache_creates_no_lock_file(self):
        p = self.p('fresh', 'cache.jsonl')
        FeatureCache(p, 'k')
        self.assertFalse(os.path.exists(p + '.lock'))
        self.assertFalse(os.path.exists(os.path.dirname(p)) and os.listdir(os.path.dirname(p)))

    def test_cache_in_a_read_only_folder_is_read_not_repaired(self):
        import stat
        row = json.dumps({'rxn': 'A>>B', '_db': 'k'})
        p = self.write('ro/cache.jsonl', row + '\n{"rxn":')              # one complete row and an interrupted append
        before = open(p).read()
        folder = os.path.dirname(p)
        os.chmod(folder, stat.S_IRUSR | stat.S_IXUSR)
        try:
            c = FeatureCache(p, 'k')
        finally:
            os.chmod(folder, stat.S_IRWXU)
        self.assertEqual(list(c.rows), ['A>>B'])
        self.assertTrue(c.recovered_truncated_tail)
        self.assertEqual(open(p).read(), before, 'a read-only cache must not be modified')
