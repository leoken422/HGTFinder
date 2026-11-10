"""Exercise the CLI using saved hits; never execute DIAMOND or network calls."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from hgtfinder import io

ROOT = Path(__file__).resolve().parents[1]


class SavedSummary(unittest.TestCase):
    def test_default_cli_veto_preserves_scores_and_zero_hit_queries(self):
        ex = ROOT / 'examples'
        with tempfile.TemporaryDirectory(dir=ROOT / 'tests') as temp:
            root = Path(temp)
            db = root / 'synthetic_manifest_only_db'
            db.mkdir()
            # Summarize only validates this file's digest; it never opens an index.
            (db / 'reference.dmnd').write_bytes(b'SYNTHETIC PLACEHOLDER - NOT A DIAMOND DATABASE')
            shutil.copyfile(ex / 'reference_metadata.tsv', db / 'reference_metadata.tsv')
            io.dump(db / 'db_manifest.json', dict(
                taxonomy_source='GTDB', protein_scope='all-cds',
                gtdb_release='fixture-not-real-GTDB', selection='synthetic-fixture',
                n_genomes=14, n_proteins=15,
                database=io.digest(db / 'reference.dmnd'),
                metadata=io.digest(db / 'reference_metadata.tsv')))
            args = [sys.executable, '-m', 'hgtfinder', 'summarize',
                    '--query', str(ex / 'query.faa'), '--query-taxonomy', str(ex / 'query_taxonomy.tsv'),
                    '--db', str(db), '--hits', str(ex / 'saved_diamond.tsv')]
            raw_before = io.digest(ex / 'saved_diamond.tsv')
            subprocess.run(args + ['--out', str(root / 'standard')], check=True, capture_output=True, text=True)
            config = root / 'baseline.toml'
            config.write_text((ROOT / 'config/v0.1.toml').read_text().replace(
                'exclude_best_hit_same_genus = true', 'exclude_best_hit_same_genus = false'))
            subprocess.run(args[:3] + ['--config', str(config)] + args[3:] +
                           ['--out', str(root / 'baseline')], check=True, capture_output=True, text=True)
            standard = {r['query_id']: r for r in io.table(root / 'standard/genes.tsv', ['query_id'])}
            baseline = {r['query_id']: r for r in io.table(root / 'baseline/genes.tsv', ['query_id'])}
            self.assertEqual(set(standard), {'q1', 'no_hit'})
            self.assertEqual(standard['q1']['status'], 'non_candidate')
            self.assertEqual(json.loads(standard['q1']['reason_codes']), ['best_hit_same_genus'])
            self.assertEqual(standard['q1']['base_distance_status'], 'candidate')
            self.assertEqual(baseline['q1']['status'], 'candidate')
            for key in ('taxonomy_distance_score', 'independent', 'all_hits', 'n_selected'):
                self.assertEqual(standard['q1'][key], baseline['q1'][key])
            for file in ('hits_top100.tsv', 'hits_audit.tsv', 'taxonomy_counts.tsv'):
                self.assertEqual((root / 'standard' / file).read_bytes(), (root / 'baseline' / file).read_bytes())
            self.assertEqual(standard['no_hit']['status'], 'pending')
            self.assertEqual(standard['no_hit']['taxonomy_distance_score'], 'NA')
            self.assertEqual(io.digest(ex / 'saved_diamond.tsv'), raw_before)
