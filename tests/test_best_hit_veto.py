import copy
from pathlib import Path
import unittest
from hgtfinder import io
from hgtfinder.scoring import apply_best_hit_genus_veto, score_hits

ROOT = Path(__file__).resolve().parents[1]


class BestHitVetoTests(unittest.TestCase):
    def setUp(self):
        ex = ROOT / 'examples'
        self.cfg = io.config()
        self.cfg['scoring']['exclude_best_hit_same_genus'] = False
        self.meta = io.metadata(ex / 'reference_metadata.tsv')
        self.query = io.table(ex / 'query_taxonomy.tsv', ['sample_id'])[0]
        self.hits = io.hits(ex / 'debug_hits.tsv', ['q1'], self.meta, True)['q1']

    def test_veto_never_replaces_pending_or_unknown_with_negative(self):
        self.assertEqual(apply_best_hit_genus_veto('candidate', ['distance_threshold_met'], True),
                         ('non_candidate', ['best_hit_same_genus']))
        for same in (False, None):
            self.assertEqual(apply_best_hit_genus_veto('candidate', ['distance_threshold_met'], same)[0], 'candidate')
        self.assertEqual(apply_best_hit_genus_veto('pending', ['insufficient_scored_hits'], True),
                         ('pending', ['insufficient_scored_hits']))
        self.assertEqual(apply_best_hit_genus_veto('non_candidate', ['distance_below_threshold'], True)[1],
                         ['distance_below_threshold'])

    def test_production_rule_preserves_hits_scores_and_denominators(self):
        baseline = score_hits(self.hits, self.query, self.meta, self.cfg)
        cfg = copy.deepcopy(self.cfg)
        cfg['scoring']['exclude_best_hit_same_genus'] = True
        result = score_hits(self.hits, self.query, self.meta, cfg)
        self.assertEqual(baseline['evidence']['status'], 'candidate')
        self.assertEqual(result['evidence']['status'], 'non_candidate')
        self.assertTrue(result['evidence']['best_hit_same_genus'])
        self.assertTrue(result['evidence']['genus_veto_applied'])
        for key in ('independent', 'all_hits', 'taxonomy_distance_score', 'n_selected'):
            self.assertEqual(result['evidence'][key], baseline['evidence'][key])
        for key in ('selected', 'audit', 'taxonomy_counts'):
            self.assertEqual(result[key], baseline[key])
        self.assertNotIn('genus_veto_applied', baseline['evidence'])

    def test_distal_best_keeps_candidate_and_no_hits_stay_pending(self):
        cfg = copy.deepcopy(self.cfg)
        cfg['scoring']['exclude_best_hit_same_genus'] = True
        hits = copy.deepcopy(self.hits)
        for h in hits:
            if h['sseqid'] == 's6':
                h['bitscore'] = 1000
        e = score_hits(hits, self.query, self.meta, cfg)['evidence']
        self.assertEqual(e['best_hit'], 's6')
        self.assertEqual(e['status'], 'candidate')
        self.assertFalse(e['genus_veto_applied'])
        e = score_hits([], self.query, self.meta, cfg)['evidence']
        self.assertEqual(e['status'], 'pending')
        self.assertIsNone(e['best_hit_same_genus'])


if __name__ == '__main__':
    unittest.main()
