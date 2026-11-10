"""Validated input and reproducibility helpers."""
import csv
import hashlib
import json
import math
from pathlib import Path
import tomllib
from .scoring import RANKS, validate_lineages

FIELDS = ['qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore']
DEFAULT_CONFIG = Path(__file__).resolve().parents[2] / 'config' / 'v0.1.toml'
if not DEFAULT_CONFIG.exists():
    DEFAULT_CONFIG = Path(__file__).with_name('v0.1.toml')


def config(path=None):
    with open(path or DEFAULT_CONFIG, 'rb') as handle:
        cfg = tomllib.load(handle)
    s, c = cfg['search'], cfg['scoring']
    if not isinstance(c.get('exclude_best_hit_same_genus'), bool):
        raise ValueError('exclude_best_hit_same_genus must be true or false')
    if not (isinstance(s['top_k'], int) and 0 < s['top_k'] <= s['pool_size']):
        raise ValueError('Require 0 < top_k <= pool_size')
    if s['threads'] < 1 or s['evalue'] <= 0 or any(not 0 <= s[k] <= 100 for k in ('query_coverage','subject_coverage')):
        raise ValueError('Invalid search thresholds')
    if s['sensitivity'] not in ('fast','mid-sensitive','sensitive','more-sensitive','very-sensitive','ultra-sensitive'):
        raise ValueError('Invalid sensitivity')
    if c['close_rank'] not in RANKS or len(c['points']) != 8 or c['points'][0] != 0 or c['points'] != sorted(c['points']):
        raise ValueError('Invalid distance configuration')
    if c['minimum_hits'] < 1 or c['minimum_genomes'] < 1:
        raise ValueError('Invalid minimum support')
    if any(not math.isfinite(v) for v in c['points'] + [c['candidate_mean'], s['evalue'], s['query_coverage'], s['subject_coverage']]):
        raise ValueError('Nonfinite configuration')
    return cfg


def table(path, required, key=None):
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or not set(required) <= set(reader.fieldnames) or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f'{path}: required columns {required}')
        rows = list(reader)
    seen = set()
    for row in rows:
        if None in row or any(row[k] is None or not row[k].strip() for k in required):
            raise ValueError(f'{path}: malformed/empty required field')
        if key:
            if row[key] in seen:
                raise ValueError(f'{path}: duplicate {key}: {row[key]}')
            seen.add(row[key])
    return rows


def fasta(path):
    records, name, seq = {}, None, []
    def finish():
        if name is not None:
            sequence = ''.join(seq)
            if not sequence or any(c not in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*-' for c in sequence.upper()):
                raise ValueError(f'Invalid protein sequence: {name}')
            records[name] = sequence
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                finish()
                if not line[1:].split():
                    raise ValueError('Empty FASTA identifier')
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError('Duplicate FASTA ID: '+name)
                seq = []
            elif name is None:
                raise ValueError('FASTA sequence before header')
            else:
                seq.append(line)
    finish()
    if not records:
        raise ValueError('Empty FASTA')
    return records


def metadata(path):
    rows = table(path, ['subject_id','genome_id','gtdb_lineage'], 'subject_id')
    validate_lineages(r['gtdb_lineage'] for r in rows)
    genomes = {}
    for r in rows:
        if r['genome_id'] in genomes and genomes[r['genome_id']] != r['gtdb_lineage']:
            raise ValueError('Inconsistent genome lineage: '+r['genome_id'])
        genomes[r['genome_id']] = r['gtdb_lineage']
    return {r['subject_id']:r for r in rows}


def query_taxonomy(path, release):
    rows = table(path, ['sample_id','gtdb_release','gtdb_lineage'])
    if len(rows) != 1 or rows[0]['gtdb_release'] != release:
        raise ValueError('Require one query taxonomy row and matching GTDB release')
    validate_lineages([rows[0]['gtdb_lineage']])
    return rows[0]


def iter_hits(path, query_ids, meta, header=False):
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t', fieldnames=None if header else FIELDS)
        if header and not set(FIELDS) <= set(reader.fieldnames or []):
            raise ValueError('Missing hit columns')
        for row in reader:
            if None in row or any(row.get(k) is None for k in FIELDS):
                raise ValueError('Invalid DIAMOND row')
            if row['qseqid'] not in query_ids or row['sseqid'] not in meta:
                raise ValueError('Unmapped query/subject: '+row['qseqid']+'/'+row['sseqid'])
            h = {k: row[k] if k in ('qseqid','sseqid') else float(row[k]) for k in FIELDS}
            if any(not math.isfinite(v) for v in h.values() if isinstance(v,float)):
                raise ValueError('Nonfinite hit')
            if h['qlen'] <= 0 or h['slen'] <= 0 or h['length'] <= 0 or h['evalue'] < 0 or h['bitscore'] < 0 or not 0 <= h['pident'] <= 100:
                raise ValueError('Invalid hit metrics')
            for prefix in ('q','s'):
                if not 1 <= min(h[prefix+'start'],h[prefix+'end']) <= max(h[prefix+'start'],h[prefix+'end']) <= h[prefix+'len']:
                    raise ValueError('Invalid alignment coordinates')
            h['qcov'] = 100*(abs(h['qend']-h['qstart'])+1)/h['qlen']
            h['scov'] = 100*(abs(h['send']-h['sstart'])+1)/h['slen']
            yield h



def hits(path, query_ids, meta, header=False):
    result = {q: [] for q in query_ids}
    for h in iter_hits(path, result, meta, header):
        result[h['qseqid']].append(h)
    return result


class TSVSink:
    def __init__(self, path, fields):
        self.handle = open(path, 'w', newline='')
        self.writer = csv.DictWriter(self.handle, fields, delimiter='\t', lineterminator='\n')
        self.writer.writeheader()

    def write(self, row):
        self.writer.writerow({k: 'NA' if v is None else json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})

    def close(self):
        self.handle.close()


def digest(path):
    h = hashlib.sha256()
    with open(path,'rb') as handle:
        for chunk in iter(lambda:handle.read(1024*1024), b''):
            h.update(chunk)
    return dict(path=str(Path(path).resolve()), sha256=h.hexdigest())


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    with open(path,'w',newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({k:'NA' if v is None else json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})
