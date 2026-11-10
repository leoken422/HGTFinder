"""Small real DIAMOND process test against explicitly synthetic proteins."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from hgtfinder import io

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
EX = ROOT/'examples'
DIAMOND = Path(os.environ.get('HGTFINDER_TEST_DIAMOND', '/nonexistent/diamond'))

@unittest.skipUnless(DIAMOND.exists(), 'DIAMOND integration requires explicit HGTFINDER_TEST_DIAMOND')
class PipelineTest(unittest.TestCase):
    def cli(self,*args):
        subprocess.run([sys.executable,'-m','hgtfinder',*map(str,args)],check=True,capture_output=True,text=True)

    def test_search_resummarize_and_debug_match(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'tests') as temp:
            out=Path(temp); db=out/'db'
            self.cli('prepare-db','--reference',EX/'reference.faa','--metadata',EX/'reference_metadata.tsv',
                '--release','fixture-not-real-GTDB','--selection','synthetic-fixture','--provenance','synthetic, not real GTDB',
                '--protein-scope','all-cds','--diamond',DIAMOND,'--out',db)
            common=['--query',EX/'query.faa','--query-taxonomy',EX/'query_taxonomy.tsv','--db',db]
            self.cli('run',*common,'--diamond',DIAMOND,'--out',out/'run')
            raw=out/'run/search/diamond.tsv'; before=io.digest(raw)
            self.cli('summarize',*common,'--hits',raw,'--search-manifest',out/'run/search/search_manifest.json','--out',out/'again')
            for name in ['genes.tsv','hits_top100.tsv','hits_audit.tsv','taxonomy_counts.tsv','summary.json']:
                self.assertEqual((out/'run/summary'/name).read_bytes(),(out/'again'/name).read_bytes())
            self.assertEqual(io.digest(raw),before)
            rows=io.table(out/'again/genes.tsv',['query_id'])
            self.assertEqual({r['query_id'] for r in rows},{'q1','no_hit'})
            nohit=next(r for r in rows if r['query_id']=='no_hit')
            self.assertEqual(nohit['taxonomy_distance_score'],'NA')
            self.assertEqual(nohit['status'],'pending')
            self.cli('summarize',*common,'--hits',EX/'saved_diamond.tsv','--out',out/'saved')
            q=io.table(out/'saved/genes.tsv',['query_id'])[0]
            expected=json.loads((EX/'debug_expected.json').read_text())['q1']['evidence']
            for k,v in expected.items():
                actual=q[k]
                if isinstance(v,(dict,list)):
                    actual=json.loads(actual)
                    if k=='flags': actual.remove('search_provenance_unverified')
                    self.assertEqual(actual,json.loads(json.dumps(v)))
                else:
                    self.assertEqual(actual,'NA' if v is None else str(v))
            m=json.loads((out/'again/run_manifest.json').read_text())
            self.assertIn('diamond_version',m['search_manifest'])
            self.assertEqual(m['config'],io.config())
            # Reject stale DB metadata rather than silently changing classifications.
            with (db/'reference_metadata.tsv').open('a') as f: f.write('\n')
            with self.assertRaises(subprocess.CalledProcessError):
                self.cli('summarize',*common,'--hits',raw,'--out',out/'bad')
