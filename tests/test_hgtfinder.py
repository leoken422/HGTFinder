import copy
import json
from pathlib import Path
import tempfile
import unittest
from hgtfinder import io
from hgtfinder.scoring import distance, lineage, score_hits, stats, validate_lineages, decide

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT/'examples'

class ScoringTests(unittest.TestCase):
    def setUp(self):
        self.cfg = io.config(ROOT/'config/v0.1.toml')
        self.meta = io.metadata(EX/'reference_metadata.tsv')
        self.q = io.table(EX/'query_taxonomy.tsv',['sample_id'])[0]
        self.raw = io.hits(EX/'debug_hits.tsv',['q1','no_hit'],self.meta,True)

    def test_distance_all_ranks_missing_and_conflicts(self):
        for i in range(8):
            self.assertEqual(distance(self.q['gtdb_lineage'],self.meta[f's{i}']['gtdb_lineage'],self.cfg['scoring']['points'])[0],i)
        self.assertEqual(distance('bad','bad',list(range(8)),True)[0],0)
        self.assertIsNone(distance(self.q['gtdb_lineage'],self.meta['unknown']['gtdb_lineage'],list(range(8)))[0])
        # Divergence at known domain fixes distance even when all lower ranks missing.
        self.assertEqual(distance(self.q['gtdb_lineage'],'d__Archaea;p__;c__;o__;f__;g__;s__',list(range(8)))[0],7)
        with self.assertRaises(ValueError):
            validate_lineages([self.q['gtdb_lineage'],self.q['gtdb_lineage'].replace('p__P','p__Conflict')])
        with self.assertRaises(ValueError):
            lineage('g__G')

    def test_hand_calculation_and_expected_debug(self):
        r=score_hits(self.raw['q1'],self.q,self.meta,self.cfg)
        e=r['evidence']
        self.assertEqual(e['independent']['numerator'],34)
        self.assertEqual(e['independent']['denominator'],11)
        self.assertEqual(e['taxonomy_distance_score'],34/11)
        self.assertEqual(e['all_hits']['numerator'],36)
        self.assertEqual(e['all_hits']['denominator'],12)
        self.assertEqual(e['n_unknown'],1)
        self.assertEqual(e['exclusions'],{'self_genome':1,'duplicate_subject':1,'quality_filter':1})
        self.assertEqual(e['status'],'non_candidate')
        self.assertEqual(e['base_distance_status'],'candidate')
        self.assertTrue(e['genus_veto_applied'])
        self.assertEqual(json.loads(json.dumps(r)),json.loads((EX/'debug_expected.json').read_text())['q1'])
        for mode,n in [('protein',13),('genome_best',12)]:
            for rank in ('domain','phylum','class','order','family','genus','species'):
                self.assertEqual(sum(v['count'] for v in r['taxonomy_counts'] if v['mode']==mode and v['rank']==rank),n)
        self.assertLessEqual(stats([0,1,2])['mean'],stats([0,2,2])['mean'])

    def test_no_hits_reference_absence_and_boundary(self):
        r=score_hits([],self.q,self.meta,self.cfg)['evidence']
        self.assertEqual(r['status'],'pending'); self.assertIsNone(r['taxonomy_distance_score'])
        self.assertEqual(r['n_selected'],0)
        cfg=self.cfg['scoring']
        self.assertEqual(decide(2,10,5,True,1,cfg)[0],'candidate')
        self.assertEqual(decide(1.999,10,5,True,1,cfg)[0],'non_candidate')
        self.assertEqual(decide(7,9,5,True,1,cfg)[0],'pending')
        self.assertEqual(decide(7,10,4,True,1,cfg)[0],'pending')
        self.assertIn('close_reference_absent',decide(7,10,5,True,0,cfg)[1])
        self.assertIn('query_close_unknown',decide(7,10,5,False,1,cfg)[1])

    def test_top100_deterministic_protein_not_genome(self):
        meta={}; hits=[]
        for i in range(105):
            sid=f'p{i:03}'
            meta[sid]=dict(subject_id=sid,genome_id='g2',gtdb_lineage=self.meta['s2']['gtdb_lineage'])
            h=dict(self.raw['q1'][0],sseqid=sid,bitscore=100)
            hits.append(h)
        r=score_hits(hits[::-1],self.q,meta,self.cfg)
        self.assertEqual([h['sseqid'] for h in r['selected']],sorted(meta)[:100])
        self.assertEqual(r['evidence']['n_unique_genomes'],1)
        self.assertIn('boundary_bitscore_tie',r['evidence']['flags'])
        cfg=copy.deepcopy(self.cfg); cfg['search']['pool_size']=100
        self.assertIn('search_pool_limit_reached',score_hits(hits,self.q,meta,cfg)['evidence']['flags'])

    def test_input_contracts(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'tests') as d:
            p=Path(d)/'bad.faa'; p.write_text('>a\nAA\n>a\nAA\n')
            with self.assertRaises(ValueError): io.fasta(p)
            p.write_text('subject_id\tgenome_id\tgtdb_lineage\nx\tg\t'+self.q['gtdb_lineage']+'\nx\tg\t'+self.q['gtdb_lineage']+'\n')
            with self.assertRaises(ValueError): io.metadata(p)
            p.write_text('x\tmissing\t100\t10\t10\t10\t1\t10\t1\t10\t0\t20\n')
            with self.assertRaises(ValueError): io.hits(p,['x'],self.meta)

if __name__ == '__main__': unittest.main()
