import os
import gzip
import hashlib
import io as bytesio
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest
from hgtfinder import io
from hgtfinder.compact import CompactMetadata,make_context
from hgtfinder.full_db import convert_proteome
from hgtfinder.scoring import score_hits

ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT.parents[1]
DIAMOND=Path(os.environ.get('HGTFINDER_TEST_DIAMOND', '/nonexistent/diamond'))

class CompactTests(unittest.TestCase):
    def test_unique_ordinals_and_same_scoring(self):
        data=b'>same_name source1\nACDEFG\n>same_name source2\nACDEFG\n'
        fa,mapping,n,letters=convert_proteome(data,'g1')
        self.assertEqual(n,2);self.assertEqual(letters,12)
        self.assertIn(b'>g1_p2\n',fa);self.assertEqual(mapping.count(b'same_name'),2)
        with self.assertRaises(ValueError):convert_proteome(b'>bad\n123\n','g1')
        with tempfile.TemporaryDirectory(dir=ROOT/'tests') as tmp:
            p=Path(tmp)/'genomes.tsv'
            tax='d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__Q'
            far=tax.replace('g__G','g__Other').replace('s__Q','s__Other')
            io.write_tsv(p,[dict(genome_key='g1',genome_id='self',gtdb_lineage=tax,protein_count=2),dict(genome_key='g2',genome_id='other',gtdb_lineage=far,protein_count=1)])
            meta=CompactMetadata(p)
            for bad in ('g1_p0','g1_p3','g3_p1','g1_p01','not_a_subject'):
                self.assertNotIn(bad,meta)
            hits=[dict(qseqid='q',sseqid=sid,pident=90,length=6,qlen=6,slen=6,qstart=1,qend=6,sstart=1,send=6,evalue=1e-9,bitscore=100-i,qcov=100,scov=100) for i,sid in enumerate(('g1_p1','g1_p2','g2_p1'))]
            query=dict(genome_id='self',gtdb_lineage=tax);cfg=io.config()
            compact=score_hits(hits,query,meta,cfg,make_context(query,meta,cfg))
            expanded=score_hits(hits,query,{h['sseqid']:meta[h['sseqid']] for h in hits},cfg)
            self.assertEqual(compact,expanded)
            self.assertEqual(compact['evidence']['exclusions']['self_genome'],2)

    @unittest.skipUnless(DIAMOND.exists(),'DIAMOND integration requires explicit HGTFINDER_TEST_DIAMOND')
    def test_stream_archive_to_diamond_and_cli(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'tests') as tmp:
            d=Path(tmp);archive=d/'fixture.tar.gz';meta=d/'metadata.tsv.gz';taxonomy=d/'taxonomy.tsv.gz'
            sequence=io.fasta(ROOT/'examples/query.faa')['q1']
            tax='d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__Q'
            fields=['accession','gtdb_representative','gtdb_taxonomy','protein_count','coding_bases','ncbi_organism_name']
            rows=[dict(accession=g,gtdb_representative='t',gtdb_taxonomy=tax,protein_count=1,coding_bases=600,ncbi_organism_name=g) for g in ('query_genome','reference_genome')]
            import csv
            with gzip.open(meta,'wt') as f:
                w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
            with gzip.open(taxonomy,'wt') as f:
                f.write(''.join(r['accession']+'\t'+tax+'\n' for r in rows))
            with tarfile.open(archive,'w:gz') as t:
                for row in rows:
                    payload=gzip.compress(('>original_protein\n'+sequence+'\n').encode())
                    member=tarfile.TarInfo('protein_faa_reps/'+row['accession']+'_protein.faa.gz');member.size=len(payload)
                    t.addfile(member,bytesio.BytesIO(payload))
            md5=d/'MD5SUM.txt';md5.write_text(''.join(hashlib.md5(p.read_bytes()).hexdigest()+'  ./'+p.name+'\n' for p in (meta,taxonomy))+'0'*32+'  ./genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz\n')
            cmd=[sys.executable,'-m','hgtfinder.full_db','--metadata',str(meta),'--taxonomy',str(taxonomy),'--md5sum',str(md5),'--out',str(d/'db'),'--diamond',str(DIAMOND),'--local-prefix',str(archive),'--sample-genomes','2','--shard-letters','100']
            subprocess.run(cmd,check=True,capture_output=True,text=True)
            db=json.loads((d/'db/db_manifest.json').read_text());self.assertEqual(db['n_proteins'],2);self.assertFalse(db['archive_integrity_verified']);self.assertEqual(len(db['shards']),2)
            query=d/'query.faa';query.write_text('>q1\n'+sequence+'\n>zero\n'+'W'*200+'\n')
            qt=d/'query_taxonomy.tsv';io.write_tsv(qt,[dict(sample_id='synthetic',genome_id='query_genome',gtdb_release='R232',gtdb_lineage=tax)])
            subprocess.run([sys.executable,'-m','hgtfinder','run','--query',str(query),'--query-taxonomy',str(qt),'--db',str(d/'db'),'--diamond',str(DIAMOND),'--out',str(d/'run')],check=True,capture_output=True,text=True)
            search=json.loads((d/'run/search/search_manifest.json').read_text())
            self.assertTrue(all('--dbsize' in e['command'] for e in search['shards']))
            self.assertEqual(search['total_reference_letters'],400)
            result=io.table(d/'run/summary/genes.tsv',['query_id'])
            self.assertEqual(len(result),2);self.assertEqual(result[0]['n_selected'],'1');self.assertEqual(result[1]['taxonomy_distance_score'],'NA')
            self.assertEqual(json.loads(result[0]['exclusions'])['self_genome'],1)
            # Same complete reference in one file and two shards must score alike.
            mono=cmd[:-2].copy();mono[mono.index('--out')+1]=str(d/'mono_db')
            subprocess.run(mono,check=True,capture_output=True,text=True)
            subprocess.run([sys.executable,'-m','hgtfinder','run','--query',str(query),'--query-taxonomy',str(qt),'--db',str(d/'mono_db'),'--diamond',str(DIAMOND),'--out',str(d/'mono_run')],check=True,capture_output=True,text=True)
            self.assertEqual((d/'run/summary/genes.tsv').read_bytes(),(d/'mono_run/summary/genes.tsv').read_bytes())
