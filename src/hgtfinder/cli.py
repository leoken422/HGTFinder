import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sqlite3
from contextlib import ExitStack
import subprocess
import sys
from . import __version__
from . import io
from .scoring import score_hits, validate_lineages, reference_context
from .compact import load_metadata, make_context

LIMITATIONS = '未検証の探索設定。HGT確率・統計的有意性ではない。DB分類群偏り・パラログは未補正。top100に近縁ヒットがなくても近縁群での遺伝子欠如を意味しない。非候補は垂直遺伝の証明ではない。'


def now():
    return datetime.now(timezone.utc).isoformat()


def manifest(kind, cfg, inputs):
    source = Path(__file__).parent
    return dict(operation=kind, started=now(), hgtfinder_version=__version__,
        code={p.name:io.digest(p) for p in sorted(source.glob('*.py'))},
        config=cfg, inputs={k:io.digest(v) for k,v in inputs.items() if v},
        argv=sys.argv, python=sys.version)


def command(args, log):
    with open(log,'a') as handle:
        handle.write(json.dumps(args)+'\n')
        handle.flush()
        subprocess.run(args, stdout=handle, stderr=subprocess.STDOUT, check=True)


def version(binary):
    return subprocess.check_output([binary,'version'], text=True).strip()


def new_dir(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    if any(path.iterdir()):
        raise ValueError('Output directory must be empty: '+str(path))
    return path


def prepare(a, cfg):
    if a.protein_scope != 'all-cds':
        raise ValueError('Marker-only references are not supported')
    meta, seq = io.metadata(a.metadata), io.fasta(a.reference)
    if set(meta) != set(seq):
        raise ValueError('Reference FASTA and metadata IDs differ')
    out = new_dir(a.out)
    m = manifest('prepare-db',cfg,dict(reference=a.reference, metadata=a.metadata))
    m.update(gtdb_release=a.release, selection=a.selection, taxonomy_source='GTDB', protein_scope=a.protein_scope,
             provenance=a.provenance, n_proteins=len(meta), n_genomes=len({r['genome_id'] for r in meta.values()}), diamond_version=version(a.diamond))
    shutil.copyfile(a.metadata,out/'reference_metadata.tsv')
    cmd = [a.diamond,'makedb','--in',str(Path(a.reference).resolve()),'--db',str(out/'reference'),'--threads',str(cfg['search']['threads'])]
    command(cmd,out/'diamond.log')
    m.update(command=cmd, ended=now(), database=io.digest(out/'reference.dmnd'), metadata=io.digest(out/'reference_metadata.tsv'))
    io.dump(out/'db_manifest.json',m)


def db_info(path):
    root = Path(path)
    m = json.loads((root/'db_manifest.json').read_text())
    if m['taxonomy_source'] != 'GTDB' or m['protein_scope'] != 'all-cds':
        raise ValueError('Not a GTDB all-CDS database')
    if m.get('metadata_layout')=='genome_ordinal_v1':
        records=m['shards'] if 'shards' in m else [dict(m['database'],file='reference.dmnd')]
        for record in records:
            stat=(root/record['file']).stat()
            if stat.st_size!=record['bytes'] or stat.st_mtime_ns!=record['mtime_ns']:
                raise ValueError('Database shard changed; re-verify full checksum before use')
        if io.digest(root/m['metadata_filename'])['sha256']!=m['metadata']['sha256']:
            raise ValueError('Compact metadata integrity check failed')
    else:
        for key, filename in [('database','reference.dmnd'),('metadata','reference_metadata.tsv')]:
            if io.digest(root/filename)['sha256'] != m[key]['sha256']:
                raise ValueError('Database integrity check failed: '+filename)
    return root,m


def search(a,cfg):
    root, db = db_info(a.db)
    if db.get('shards'):
        from .sharded import search_shards
        return search_shards(a,cfg,root,db)
    seq = io.fasta(a.query)
    q = io.query_taxonomy(a.query_taxonomy,db['gtdb_release'])
    validate_lineages([q['gtdb_lineage']] + [m['gtdb_lineage'] for m in load_metadata(root,db).values()])
    out = new_dir(a.out)
    m = manifest('search',cfg,dict(query=a.query,query_taxonomy=a.query_taxonomy,db_manifest=root/'db_manifest.json'))
    s = cfg['search']
    cmd = [a.diamond,'blastp','--query',str(Path(a.query).resolve()),'--db',str(root/'reference.dmnd'),
        '--out',str(out/'diamond.tsv'),'--outfmt','6',*io.FIELDS,'--max-target-seqs',str(s['pool_size']),
        '--max-hsps','1','--evalue',str(s['evalue']),'--threads',str(s['threads']),'--'+s['sensitivity'], '--tmpdir',str(out/'tmp')]
    for key,flag in [('block_size','--block-size'),('index_chunks','--index-chunks')]:
        if key in s:cmd += [flag,str(s[key])]
    (out/'tmp').mkdir()
    m.update(diamond_version=version(a.diamond), command=cmd, db=db, n_queries=len(seq))
    command(cmd,out/'diamond.log')
    m.update(ended=now(),raw_output=io.digest(out/'diamond.tsv'))
    io.dump(out/'search_manifest.json',m)


def summarize(a,cfg):
    root,db = db_info(a.db)
    seq, meta = io.fasta(a.query), load_metadata(root,db)
    q = io.query_taxonomy(a.query_taxonomy,db['gtdb_release'])
    validate_lineages([q['gtdb_lineage']] + [m['gtdb_lineage'] for m in meta.values()])
    search_m = None
    if a.search_manifest:
        search_m = json.loads(Path(a.search_manifest).read_text())
        if search_m['raw_output']['sha256'] != io.digest(a.hits)['sha256'] or search_m['inputs']['query']['sha256'] != io.digest(a.query)['sha256'] or search_m['inputs']['db_manifest']['sha256'] != io.digest(root/'db_manifest.json')['sha256'] or search_m['inputs']['query_taxonomy']['sha256'] != io.digest(a.query_taxonomy)['sha256']:
            raise ValueError('Search manifest input mismatch')
        old = search_m['config']['search']
        for k in ('pool_size','evalue','sensitivity','block_size','index_chunks'):
            if cfg['search'].get(k) != old.get(k):
                raise ValueError('Search-affecting setting changed; rerun search: '+k)
    annotations = {}
    if a.genes:
        annotations = {r['query_id']:r for r in io.table(a.genes,['query_id'],'query_id')}
        if not set(annotations) <= set(seq):
            raise ValueError('Annotation contains unknown query IDs')
    out = new_dir(a.out)
    m = manifest('summarize',cfg,dict(query=a.query,query_taxonomy=a.query_taxonomy,hits=a.hits,db_manifest=root/'db_manifest.json', genes=a.genes,search_manifest=a.search_manifest))
    context = make_context(q,meta,cfg)
    # Stage by query on disk so millions of pool hits never accumulate in memory.
    staging = out/'hits_staging.sqlite'
    conn = sqlite3.connect(staging)
    conn.execute('PRAGMA temp_store=MEMORY')
    try:
        conn.execute('CREATE TABLE hits (query_id TEXT, payload TEXT)')
        with conn:
            for h in io.iter_hits(a.hits,seq,meta):
                if h['qlen'] != len(seq[h['qseqid']]):
                    raise ValueError('Search query length mismatch: '+h['qseqid'])
                conn.execute('INSERT INTO hits VALUES (?,?)',(h['qseqid'],json.dumps(h)))
        conn.execute('CREATE INDEX query_idx ON hits(query_id)')
        hit_fields = [*io.FIELDS,'qcov','scov',*dict.fromkeys(k for r in meta.values() for k in r),
            'query_lineage','distance_points','distance_reason','common_rank','proximity','selected','exclusion_reason']
        hit_fields = list(dict.fromkeys(hit_fields))
        statuses = Counter()
        with ExitStack() as stack:
            def sink(filename, fields):
                writer = io.TSVSink(out/filename,fields)
                stack.callback(writer.close)
                return writer
            tops = sink('hits_top100.tsv',hit_fields)
            audits = sink('hits_audit.tsv',hit_fields)
            counts = sink('taxonomy_counts.tsv',['query_id','mode','rank','taxon','count','fraction','denominator','known_denominator'])
            genes = None
            for query_id in seq:
                raw = [json.loads(row[0]) for row in conn.execute('SELECT payload FROM hits WHERE query_id=?',(query_id,))]
                result = score_hits(raw,q,meta,cfg,context)
                ev = result['evidence']
                if not search_m:
                    ev['flags'].append('search_provenance_unverified')
                record = dict(annotations.get(query_id,{}),query_id=query_id,query_length=len(seq[query_id]),**ev)
                if genes is None:
                    fields = list(dict.fromkeys([k for row in annotations.values() for k in row]+list(record)))
                    genes = sink('genes.tsv',fields)
                genes.write(record)
                statuses[ev['status']] += 1
                for row in result['taxonomy_counts']:
                    counts.write(dict(query_id=query_id,**row))
                for row in result['selected']:
                    tops.write(row)
                for row in result['audit']:
                    audits.write(row)
    finally:
        conn.close()
        staging.unlink(missing_ok=True)
    summary = dict(n_queries=len(seq),statuses=dict(statuses),db=dict(gtdb_release=db['gtdb_release'],selection=db['selection'],n_genomes=db['n_genomes'],n_proteins=db['n_proteins']),limitations=LIMITATIONS)
    io.dump(out/'summary.json',summary)
    (out/'summary.md').write_text('# HGTFinder v0.1\n\n'+json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    m.update(ended=now(),db=db,search_manifest=search_m,n_queries=len(seq),outputs={p.name:io.digest(p) for p in out.iterdir() if p.is_file()})
    io.dump(out/'run_manifest.json',m)
    (out/'run.log').write_text(f'{m["started"]} — {m["ended"]}: {len(seq)} queries summarized\n')


def debug(a,cfg):
    meta = io.metadata(a.metadata)
    rows_q = io.table(a.query_taxonomy,['sample_id','gtdb_release','gtdb_lineage'])
    if len(rows_q) != 1:
        raise ValueError('Require one query taxonomy row')
    q = rows_q[0]
    # Headered saved hits, same parser and same scoring function as production.
    rows = io.table(a.hits,io.FIELDS)
    raw = io.hits(a.hits,dict.fromkeys(r['qseqid'] for r in rows),meta,header=True)
    return {key:score_hits(value,q,meta,cfg) for key,value in raw.items()}


def main(argv=None):
    parser = argparse.ArgumentParser(description=LIMITATIONS)
    parser.add_argument('--config',default=None)
    sub = parser.add_subparsers(dest='operation',required=True)
    p = sub.add_parser('prepare-db')
    for name in ('reference','metadata','release','selection','provenance','protein-scope','out'):
        p.add_argument('--'+name,required=True)
    p.add_argument('--diamond',default='diamond')
    for name in ('search','summarize','run'):
        p = sub.add_parser(name)
        for arg in ('query','query-taxonomy','db','out'):
            p.add_argument('--'+arg,required=True)
        if name in ('search','run'):
            p.add_argument('--diamond',default='diamond')
        if name in ('summarize','run'):
            p.add_argument('--genes')
        if name == 'summarize':
            p.add_argument('--hits',required=True)
            p.add_argument('--search-manifest')
    p = sub.add_parser('debug-score')
    for name in ('hits','metadata','query-taxonomy'):
        p.add_argument('--'+name,required=True)
    a = parser.parse_args(argv)
    try:
        cfg = io.config(a.config)
        if a.operation == 'run':
            base = new_dir(a.out)
            a.out = str(base/'search')
            search(a,cfg)
            a.hits, a.search_manifest = str(base/'search/diamond.tsv'),str(base/'search/search_manifest.json')
            a.out = str(base/'summary')
            summarize(a,cfg)
        elif a.operation == 'debug-score':
            print(json.dumps(debug(a,cfg),ensure_ascii=False,indent=2,allow_nan=False))
        else:
            {'prepare-db':prepare,'search':search,'summarize':summarize}[a.operation](a,cfg)
    except (ValueError,KeyError,OSError,subprocess.CalledProcessError) as exc:
        parser.exit(2,f'error: {exc}\n')

if __name__ == '__main__':
    main()
