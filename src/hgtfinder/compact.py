"""Compact genome/ordinal IDs for very large GTDB protein references.

Only one taxonomy row per genome is held in memory. Protein IDs are validated
against each genome's observed protein count; original IDs remain in a gzip map.
"""
from . import io
from .scoring import reference_context


class CompactMetadata:
    def __init__(self,path):
        rows=io.table(path,['genome_key','genome_id','gtdb_lineage','protein_count'],'genome_key')
        self.genomes={r['genome_key']:dict(subject_id=r['genome_key'],genome_id=r['genome_id'],gtdb_lineage=r['gtdb_lineage']) for r in rows}
        self.counts={r['genome_key']:int(r['protein_count']) for r in rows}
        if len({r['genome_id'] for r in rows})!=len(rows) or any(n<1 for n in self.counts.values()):
            raise ValueError('Invalid compact genome table')

    def genome_key(self,sid):
        if len(sid)>40:raise KeyError(sid)
        try:
            key,number=sid.split('_p')
            ordinal=int(number,16)
            if key not in self.counts or not 1<=ordinal<=self.counts[key] or number!=f'{ordinal:x}':raise KeyError(sid)
            return key
        except (ValueError,TypeError):raise KeyError(sid) from None

    def __contains__(self,sid):
        try:self.genome_key(sid);return True
        except KeyError:return False

    def __getitem__(self,sid):
        return dict(self.genomes[self.genome_key(sid)],subject_id=sid)

    def values(self):
        # Deliberately genome rows, not hundreds of millions of protein records.
        return self.genomes.values()


class CompactClassifications:
    def __init__(self,metadata,by_genome):self.metadata,self.by_genome=metadata,by_genome
    def __getitem__(self,sid):return self.by_genome[self.metadata.genome_key(sid)]


def load_metadata(root,db):
    if db.get('metadata_layout')=='genome_ordinal_v1':return CompactMetadata(root/db['metadata_filename'])
    return io.metadata(root/'reference_metadata.tsv')


def make_context(query,metadata,config):
    if isinstance(metadata,CompactMetadata):
        classifications,count=reference_context(query,metadata.genomes,config)
        return CompactClassifications(metadata,classifications),count
    return reference_context(query,metadata,config)
