import copy
from pathlib import Path
import tempfile
import unittest
from hgtfinder import io
from hgtfinder.scoring import score_hits

ROOT = Path(__file__).resolve().parents[1]


class ReleaseDefaults(unittest.TestCase):
    def test_configs_enable_veto_and_bundled_config_matches(self):
        for name in ('v0.1.toml', 'full_gtdb_r232.toml'):
            self.assertIs(io.config(ROOT / 'config' / name)['scoring']['exclude_best_hit_same_genus'], True)
        self.assertEqual((ROOT / 'config/v0.1.toml').read_bytes(),
                         (ROOT / 'src/hgtfinder/v0.1.toml').read_bytes())

    def test_non_boolean_veto_rejected(self):
        with tempfile.TemporaryDirectory(dir=ROOT / 'tests') as temp:
            p = Path(temp) / 'invalid.toml'
            p.write_text((ROOT / 'config/v0.1.toml').read_text().replace(
                'exclude_best_hit_same_genus = true', 'exclude_best_hit_same_genus = "false"'))
            with self.assertRaisesRegex(ValueError, 'must be true or false'):
                io.config(p)

    def test_best_means_bitscore_not_identity_or_any_close_hit(self):
        ex = ROOT / 'examples'
        cfg = io.config()
        meta = io.metadata(ex / 'reference_metadata.tsv')
        query = io.table(ex / 'query_taxonomy.tsv', ['sample_id'])[0]
        hits = io.hits(ex / 'debug_hits.tsv', ['q1'], meta, True)['q1']
        hits = copy.deepcopy(hits)
        for h in hits:
            if h['sseqid'] == 's6':
                h.update(bitscore=1000, pident=60)
        result = score_hits(hits, query, meta, cfg)['evidence']
        self.assertEqual(result['status'], 'candidate')
        self.assertFalse(result['best_hit_same_genus'])
        self.assertGreater(result['proximity_counts']['close'], 0)

    def test_quality_failed_best_does_not_veto(self):
        ex = ROOT / 'examples'
        cfg = io.config()
        meta = io.metadata(ex / 'reference_metadata.tsv')
        query = io.table(ex / 'query_taxonomy.tsv', ['sample_id'])[0]
        hits = io.hits(ex / 'debug_hits.tsv', ['q1'], meta, True)['q1']
        hits = copy.deepcopy(hits)
        for h in hits:
            if h['sseqid'] == 's0':
                h.update(bitscore=2000, qcov=10)
            if h['sseqid'] == 's6':
                h['bitscore'] = 1000
        result = score_hits(hits, query, meta, cfg)['evidence']
        self.assertEqual(result['best_hit'], 's6')
        self.assertEqual(result['status'], 'candidate')
