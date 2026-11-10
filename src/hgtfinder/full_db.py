"""Stream the official GTDB all-species-representative protein archive to DIAMOND.

No full decompressed FASTA or full compressed archive is kept. HTTP Range chunks
are retried independently while the gzip stream stays alive. Final archive MD5,
all genome IDs and per-genome protein counts must pass before publishing the DB.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime,timezone
import gzip
import hashlib
import io as pyio
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import threading
import time
from . import io
from .cli import manifest,now
from .scoring import validate_lineages

GIB=1024**3
AA=b'ABCDEFGHIJKLMNOPQRSTUVWXYZ*-'


def atomic(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp');io.dump(tmp,value);os.replace(tmp,path)


def representatives(metadata_paths,taxonomy_paths,md5sum):
    expected={}
    for line in Path(md5sum).read_text().splitlines():
        digest,name=line.split(maxsplit=1);expected[Path(name.strip().lstrip('*')).name]=digest
    paths=[*map(Path,metadata_paths),*map(Path,taxonomy_paths)]
    verified=[]
    for p in paths:
        h=hashlib.md5()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
        if h.hexdigest()!=expected.get(p.name):raise ValueError('GTDB MD5 mismatch: '+p.name)
        verified.append(dict(io.digest(p),md5=h.hexdigest()))
    rows={}
    for p in metadata_paths:
        with gzip.open(p,'rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r['gtdb_representative']=='t':
                    rows[r['accession']]=dict(genome_id=r['accession'],gtdb_lineage=r['gtdb_taxonomy'],protein_count=int(r['protein_count']),coding_bases=int(r['coding_bases']),ncbi_organism_name=r['ncbi_organism_name'])
    matched=set()
    for p in taxonomy_paths:
        with gzip.open(p,'rt') as f:
            for line in f:
                genome,tax=line.rstrip('\n').split('\t')
                if genome in rows:
                    if tax!=rows[genome]['gtdb_lineage']:raise ValueError('GTDB lineage mismatch: '+genome)
                    matched.add(genome)
    if matched!=set(rows):raise ValueError('Representative missing from official taxonomy')
    validate_lineages(r['gtdb_lineage'] for r in rows.values())
    for i,genome in enumerate(sorted(rows),1):rows[genome]['genome_key']=f'g{i:x}'
    return rows,verified,expected


def convert_proteome(data,genome_key):
    """Pure conversion: preserve every protein, assign globally unique ordinal ID."""
    if not data.startswith(b'>'):raise ValueError('Not a protein FASTA')
    fasta=[];ids=[];letters=0;count=0
    for record in data.split(b'>')[1:]:
        header,separator,body=record.partition(b'\n')
        if not separator or not header.split():raise ValueError('Invalid protein header')
        original=header.split()[0]
        seq=b''.join(body.split()).upper()
        if not seq or seq.translate(None,AA):raise ValueError('Invalid protein sequence: '+repr(original))
        count+=1;letters+=len(seq);sid=f'{genome_key}_p{count:x}'.encode()
        fasta += [b'>',sid,b'\n',seq,b'\n']
        ids += [sid,b'\t',original,b'\t',str(len(seq)).encode(),b'\n']
    return b''.join(fasta),b''.join(ids),count,letters


class RangeSource:
    def __init__(self,urls,total,cache,chunk_bytes,workers,stop):
        self.urls,self.total,self.cache,self.chunk_bytes,self.workers,self.stop=urls,total,Path(cache),chunk_bytes,workers,stop
        self.cache.mkdir(parents=True,exist_ok=True)
        self.pool=ThreadPoolExecutor(max_workers=workers)
        self.futures={};self.current=None;self.index=0;self.position=0
        self.md5=hashlib.md5();self.sha256=hashlib.sha256()
        self.n_chunks=math.ceil(total/chunk_bytes)
        for i in range(min(workers,self.n_chunks)):self.futures[i]=self.pool.submit(self.download,i)

    def download(self,index):
        start=index*self.chunk_bytes;end=min(self.total,start+self.chunk_bytes)-1;size=end-start+1
        path=self.cache/f'chunk_{index:06d}.part';headers=self.cache/f'chunk_{index:06d}.headers';log=self.cache/f'chunk_{index:06d}.log'
        for attempt in range(4):
            if self.stop.is_set():raise RuntimeError('Download stopped')
            url=self.urls[(index+attempt)%len(self.urls)]
            cmd=['curl','--fail','--location','--silent','--show-error','--connect-timeout','20','--max-time','600','--range',f'{start}-{end}','--max-filesize',str(size),'--dump-header',str(headers),'--output',str(path),url]
            with log.open('a') as error:
                error.write(f'{now()} range={start}-{end} attempt={attempt+1} url={url}\n');error.flush()
                proc=subprocess.Popen(cmd,stderr=error)
                while proc.poll() is None:
                    if self.stop.wait(1):proc.terminate();proc.wait();raise RuntimeError('Download stopped')
                code=proc.returncode
            if code==0:
                ranges=re.findall(r'content-range:\s*bytes\s+(\d+)-(\d+)/(\d+)',headers.read_text(),re.I)
                if ranges and tuple(map(int,ranges[-1]))==(start,end,self.total) and path.stat().st_size==size:
                    return path
            if attempt<3:self.stop.wait(2*(attempt+1))
        raise RuntimeError(f'Failed HTTP range {start}-{end}; see {log}')

    def read(self,size=-1):
        if size<0:raise ValueError('Unbounded archive read is not supported')
        output=[];remaining=size
        while remaining and self.position<self.total:
            if self.stop.is_set():raise RuntimeError('Streaming stopped')
            if self.current is None:self.current=self.futures[self.index].result().open('rb')
            block=self.current.read(remaining)
            if block:
                self.md5.update(block);self.sha256.update(block);self.position+=len(block);output.append(block);remaining-=len(block)
            else:
                path=Path(self.current.name);self.current.close();self.current=None;path.unlink()
                del self.futures[self.index]
                new=self.index+self.workers
                if new<self.n_chunks:self.futures[new]=self.pool.submit(self.download,new)
                self.index+=1
        return b''.join(output)

    def close(self):
        self.stop.set()
        if self.current:self.current.close()
        self.pool.shutdown(wait=True,cancel_futures=True)


def build(a,cfg):
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'db_manifest.json').exists() or (out/'reference.building.dmnd').exists() or list(out.glob('reference_*.dmnd')):
        raise ValueError('DB output already exists')
    started=time.time();state=dict(stage='catalog_validation',started=now(),genomes_processed=0,proteins_processed=0,letters_processed=0,completed_shards=0)
    atomic(out/'progress.json',state)
    rows,verified,md5s=representatives(a.metadata,a.taxonomy,a.md5sum)
    expected_md5=md5s['gtdb_proteins_aa_reps_r232.tar.gz']
    state.update(expected_genomes=len(rows),expected_proteins=sum(r['protein_count'] for r in rows.values()),estimated_coding_letters=sum(r['coding_bases'] for r in rows.values())//3)
    target=out/'reference.building.dmnd';tmp=out/'tmp';tmp.mkdir(exist_ok=True)
    stop=threading.Event();monitor_stop=threading.Event();source=None;proc=None;monitor=None;log=None;ids=None
    hashes=ThreadPoolExecutor(max_workers=1);hash_jobs=[];completed=[];seen=set();observed=[]
    shard_letters=shard_proteins=0;shard_index=0
    command=[str(Path(a.diamond).resolve()),'makedb','--db',str(out/'reference.building'),'--threads',str(cfg['search']['threads']),'--tmpdir',str(tmp)]
    m=manifest('stream-prepare-db',cfg,{**{f'metadata_{i}':p for i,p in enumerate(a.metadata)},**{f'taxonomy_{i}':p for i,p in enumerate(a.taxonomy)},'md5sum':a.md5sum})
    m.update(gtdb_release='R232',taxonomy_source='GTDB',protein_scope='all-cds',metadata_layout='genome_ordinal_v1',metadata_filename='genomes.tsv',selection='all_GTDB_species_representatives' if not a.sample_genomes else f'official_archive_prefix_probe_{a.sample_genomes}_genomes',archive_urls=a.url,archive_bytes=a.archive_bytes,official_archive_md5=expected_md5,command=command,official_metadata=verified,diamond_version=subprocess.check_output([a.diamond,'version'],text=True).strip(),shard_letter_limit=a.shard_letters)
    def digest_shard(path,record):
        record.update(io.digest(path));stat=path.stat();record.update(bytes=stat.st_size,mtime_ns=stat.st_mtime_ns)
        return record
    def open_shard():
        nonlocal proc,log,ids,shard_letters,shard_proteins
        shard_letters=shard_proteins=0
        log=(out/f'diamond_makedb_{shard_index:05d}.log').open('w')
        proc=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=log,stderr=subprocess.STDOUT,cwd=out)
        ids=gzip.open(out/f'protein_ids_{shard_index:05d}.tsv.gz.building','wb',compresslevel=1)
        ids.write(b'subject_id\tsource_protein_id\tsequence_length\n')
    def close_shard(final=False):
        nonlocal proc,log,ids
        proc.stdin.close()
        if proc.wait()!=0:raise RuntimeError('DIAMOND makedb failed; inspect shard log')
        log.close();log=None;ids.close();ids=None
        final_path=out/('reference.dmnd' if final and shard_index==0 else f'reference_{shard_index:05d}.dmnd')
        target.rename(final_path)
        mapping=out/f'protein_ids_{shard_index:05d}.tsv.gz'
        Path(str(mapping)+'.building').rename(mapping)
        text=(out/f'diamond_makedb_{shard_index:05d}.log').read_text()
        match=re.search(r'Database letters\s+([0-9,]+)',text)
        letters=int(match.group(1).replace(',','')) if match else shard_letters
        record=dict(shard_index=shard_index,file=final_path.name,n_proteins=shard_proteins,letters=letters,source_letters=shard_letters,protein_id_map=mapping.name)
        completed.append(record);hash_jobs.append(hashes.submit(digest_shard,final_path,record.copy()))
        atomic(out/'completed_shards.json',completed);state['completed_shards']=len(completed)
    def progress():
        while not monitor_stop.wait(10):
            state.update(updated=now(),elapsed_seconds=round(time.time()-started),free_bytes=shutil.disk_usage(out).free,database_bytes=sum((out/r['file']).stat().st_size for r in completed)+(target.stat().st_size if target.exists() else 0))
            if source is not None and hasattr(source,'position'):state.update(archive_bytes_read=source.position,archive_total_bytes=a.archive_bytes)
            atomic(out/'progress.json',state)
            if shutil.disk_usage(out).free<a.reserve_gib*GIB:
                state['failure']='Storage reserve reached';stop.set()
                if proc is not None and proc.poll() is None:proc.terminate()
    try:
        source=open(a.local_prefix,'rb') if a.local_prefix else RangeSource(a.url,a.archive_bytes,out/'range_cache',a.chunk_mib*1024**2,a.download_workers,stop)
        open_shard();monitor=threading.Thread(target=progress,daemon=True);monitor.start();state['stage']='download_and_build'
        with (out/'observed_genomes.tsv').open('w') as ledger:
            writer=csv.DictWriter(ledger,['genome_key','genome_id','protein_count','letters','archive_member','shard_index'],delimiter='\t',lineterminator='\n');writer.writeheader()
            with tarfile.open(fileobj=source,mode='r|gz') as archive:
                for member in archive:
                    if stop.is_set():raise RuntimeError(state.get('failure','Stopped'))
                    if not member.isfile() or not member.name.endswith('_protein.faa.gz'):continue
                    genome=Path(member.name).name.removesuffix('_protein.faa.gz')
                    if genome not in rows or genome in seen:raise ValueError('Unexpected/duplicate archive genome: '+genome)
                    row=rows[genome]
                    with gzip.GzipFile(fileobj=archive.extractfile(member)) as proteins:data=proteins.read()
                    converted,mapping,count,letters=convert_proteome(data,row['genome_key'])
                    if count!=row['protein_count']:raise ValueError(f'Protein count mismatch for {genome}: {count} vs {row["protein_count"]}')
                    if shard_proteins and shard_letters+letters>a.shard_letters:
                        close_shard();shard_index+=1;open_shard()
                    proc.stdin.write(converted);ids.write(mapping)
                    shard_letters+=letters;shard_proteins+=count
                    seen.add(genome);observed.append(dict(row,shard_index=shard_index))
                    writer.writerow(dict(genome_key=row['genome_key'],genome_id=genome,protein_count=count,letters=letters,archive_member=member.name,shard_index=shard_index))
                    state.update(genomes_processed=len(seen),proteins_processed=state['proteins_processed']+count,letters_processed=state['letters_processed']+letters,shard_index=shard_index)
                    if len(seen)%50==0:ledger.flush()
                    if a.sample_genomes and len(seen)>=a.sample_genomes:break
        if not a.sample_genomes:
            state['stage']='archive_integrity_verification'
            while source.read(1024*1024):pass
            if source.position!=a.archive_bytes or source.md5.hexdigest()!=expected_md5:raise ValueError('Complete archive MD5 verification failed')
            if seen!=set(rows):raise ValueError(f'Missing representative genomes: {len(set(rows)-seen)}')
            state['archive_bytes_read']=source.position
            m.update(archive_md5=source.md5.hexdigest(),archive_sha256=source.sha256.hexdigest(),archive_integrity_verified=True)
        else:m.update(archive_integrity_verified=False,probe_prefix=io.digest(a.local_prefix))
        state['stage']='finalizing_diamond';close_shard(final=True)
        source.close();source=None
        io.write_tsv(out/'genomes.tsv',observed,['genome_id','gtdb_lineage','protein_count','coding_bases','ncbi_organism_name','genome_key','shard_index'])
        state['stage']='database_sha256';records=[future.result() for future in hash_jobs]
        maps=[io.digest(out/r['protein_id_map']) for r in records]
        m.update(ended=now(),n_genomes=len(seen),n_proteins=state['proteins_processed'],letters=sum(r['letters'] for r in records),source_letters=state['letters_processed'],shards=records,database_total_bytes=sum(r['bytes'] for r in records),metadata=io.digest(out/'genomes.tsv'),protein_id_maps=maps,integrity_check_on_use='Full SHA256 per shard at build; size+mtime on routine use, metadata SHA256 each use',provenance='Official GTDB R232 Prodigal all-species-representative amino-acid archive; streamed without retaining duplicate archive/FASTA. All genome memberships and protein counts verified. Genome/ordinal IDs map to original IDs in per-shard protein_ids_*.tsv.gz. All shards jointly form the complete reference.')
        io.dump(out/'db_manifest.json',m);state.update(stage='complete',ended=now(),database_bytes=m['database_total_bytes'])
    except BaseException as exc:
        state.update(stage='failed',error=str(exc),ended=now());atomic(out/'failure.json',state)
        if proc is not None and proc.poll() is None:proc.terminate();proc.wait()
        raise
    finally:
        stop.set();monitor_stop.set()
        if source is not None:source.close()
        if ids is not None:ids.close()
        if monitor is not None:monitor.join(timeout=11)
        if log is not None:log.close()
        hashes.shutdown(wait=True,cancel_futures=True);atomic(out/'progress.json',state)
    return m

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metadata',nargs='+',required=True);parser.add_argument('--taxonomy',nargs='+',required=True)
    parser.add_argument('--md5sum',required=True);parser.add_argument('--out',required=True);parser.add_argument('--diamond',required=True)
    parser.add_argument('--config');parser.add_argument('--url',action='append',default=[]);parser.add_argument('--archive-bytes',type=int,default=131946537881)
    parser.add_argument('--chunk-mib',type=int,default=64);parser.add_argument('--download-workers',type=int,default=4);parser.add_argument('--reserve-gib',type=float,default=8)
    parser.add_argument('--shard-letters',type=int,default=8_000_000_000)
    parser.add_argument('--local-prefix');parser.add_argument('--sample-genomes',type=int)
    a=parser.parse_args(argv)
    if bool(a.local_prefix)!=bool(a.sample_genomes):parser.error('local-prefix and sample-genomes are required together for a probe')
    if not a.local_prefix and not a.url:parser.error('Full archive URLs required')
    import signal
    signal.signal(signal.SIGTERM,lambda signum,frame:(_ for _ in ()).throw(KeyboardInterrupt('Terminated')))
    build(a,io.config(a.config))

if __name__=='__main__':main()
