"""Search the complete reference in bounded-memory shards, then merge globally.

Every shard uses the TOTAL database size for comparable E-values. A top-k result
from the union must be in its shard's top-k. Thus union candidates are sorted
before self/quality filters, reproducing the global protein candidate-pool rule.
"""
from collections import OrderedDict
import json
from pathlib import Path
import shutil
import tempfile
from . import io
from .compact import load_metadata
from .scoring import validate_lineages


def merge_candidates(paths,query_sequences,metadata,pool_size,out,tmp_root):
    counts={q:0 for q in query_sequences};index={q:i for i,q in enumerate(query_sequences)}
    with tempfile.TemporaryDirectory(prefix='merge_buckets_',dir=tmp_root) as tmp:
        root=Path(tmp);handles=OrderedDict()
        try:
            for path in paths:
                with open(path,'rb') as source:
                    for line in source:
                        if not line.strip():continue
                        q=line.split(b'\t',1)[0].decode()
                        if q not in index:raise ValueError('Unmapped shard query: '+q)
                        counts[q]+=1
                        if q not in handles:
                            if len(handles)>=32:handles.popitem(last=False)[1].close()
                            handles[q]=(root/f'{index[q]:06d}.tsv').open('ab')
                        handles.move_to_end(q);handles[q].write(line)
        finally:
            for f in handles.values():f.close()
        with open(out,'w') as target:
            for q in query_sequences:
                p=root/f'{index[q]:06d}.tsv'
                if not p.exists():continue
                hits=list(io.iter_hits(p,{q},metadata))
                hits.sort(key=lambda h:(-h['bitscore'],h['evalue'],h['sseqid']))
                for h in hits[:pool_size]:
                    target.write('\t'.join(str(h[k]) for k in io.FIELDS)+'\n')
                p.unlink()
    return counts


def search_shards(a,cfg,root,db):
    from .cli import manifest,now,command,version
    seq=io.fasta(a.query);meta=load_metadata(root,db);q=io.query_taxonomy(a.query_taxonomy,db['gtdb_release'])
    validate_lineages([q['gtdb_lineage']]+[m['gtdb_lineage'] for m in meta.values()])
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True);(out/'tmp').mkdir(exist_ok=True)
    inputs=dict(query=a.query,query_taxonomy=a.query_taxonomy,db_manifest=root/'db_manifest.json')
    plan=out/'shard_search_state.json'
    if plan.exists():
        m=json.loads(plan.read_text())
        if m['config']!=cfg or any(m['inputs'][k]['sha256']!=io.digest(p)['sha256'] for k,p in inputs.items()):
            raise ValueError('Cannot resume search with changed inputs/configuration')
    else:
        if any(p.name!='tmp' for p in out.iterdir()):raise ValueError('Search output contains unrelated files')
        m=manifest('search',cfg,inputs)
        m.update(diamond_version=version(a.diamond),db=db,n_queries=len(seq),total_reference_letters=db['letters'],shards=[])
        for i,record in enumerate(db['shards']):m['shards'].append(dict(index=i,db=record['file'],state='planned'))
        io.dump(plan,m)
    s=cfg['search'];paths=[]
    for i,record in enumerate(db['shards']):
        entry=m['shards'][i];raw=out/f'shard_{i:05d}.tsv';paths.append(raw)
        if entry['state']=='complete':
            if io.digest(raw)['sha256']!=entry['raw_output']['sha256']:raise ValueError('Completed shard output changed')
            continue
        if raw.exists():raw.rename(out/f'shard_{i:05d}.incomplete.tsv')
        cmd=[str(Path(a.diamond).resolve()),'blastp','--query',str(Path(a.query).resolve()),'--db',str(root/record['file']),'--dbsize',str(db['letters']),
            '--out',str(raw),'--outfmt','6',*io.FIELDS,'--max-target-seqs',str(s['pool_size']),'--max-hsps','1','--evalue',str(s['evalue']),'--threads',str(s['threads']),'--'+s['sensitivity'],'--tmpdir',str(out/'tmp')]
        for key,flag in [('block_size','--block-size'),('index_chunks','--index-chunks')]:
            if key in s:cmd.extend([flag,str(s[key])])
        entry.update(state='searching',started=now(),command=cmd);io.dump(plan,m)
        command(cmd,out/f'shard_{i:05d}.log')
        entry.update(state='complete',ended=now(),raw_output=io.digest(raw));io.dump(plan,m)
    m['stage']='merging_global_candidate_pool';io.dump(plan,m)
    counts=merge_candidates(paths,seq,meta,s['pool_size'],out/'diamond.tsv',out/'tmp')
    m.update(ended=now(),stage='complete',raw_output=io.digest(out/'diamond.tsv'),per_query_union_candidate_counts=counts,
        candidate_policy='Each shard searched with total DB letters for E-values; union sorted by bit score descending, E-value ascending, subject ID ascending; keep global pool before self/quality/top100 filtering',command=[e['command'] for e in m['shards']])
    io.dump(out/'search_manifest.json',m);io.dump(plan,m)
