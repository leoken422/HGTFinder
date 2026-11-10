"""Pure taxonomy, hit selection and candidate calculations; no I/O."""
from collections import Counter
from statistics import median

RANKS = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
PREFIXES = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')


def lineage(value):
    parts = value.split(';')
    if len(parts) != 7:
        raise ValueError('GTDB lineage requires seven ordered ranks: ' + value)
    result = []
    for prefix, part in zip(PREFIXES, parts):
        part = part.strip()
        if not part.startswith(prefix) or len(part.split('__')) != 2:
            raise ValueError('Invalid GTDB rank: ' + part)
        result.append(part[len(prefix):] or None)
    return tuple(result)


def validate_lineages(values):
    # A named taxon must not map to conflicting known ancestors.
    parents = {}
    for value in values:
        tax = lineage(value)
        for i, name in enumerate(tax):
            if name is None:
                continue
            key = (i, name)
            known = parents.setdefault(key, [None] * i)
            for j in range(i):
                if known[j] is not None and tax[j] is not None and known[j] != tax[j]:
                    raise ValueError('Conflicting GTDB lineage for ' + name)
                if tax[j] is not None:
                    known[j] = tax[j]


def distance(query, reference, points, self_hit=False):
    if self_hit:
        return 0, 'self_genome', None
    q, r = lineage(query), lineage(reference)
    common = None
    for i, (a, b) in enumerate(zip(q, r)):
        if a is None or b is None:
            return None, 'taxonomy_incomplete', common
        if a != b:
            return points[7 - i], 'scored', common
        common = RANKS[i]
    return points[0], 'scored', common


def proximity(query, reference, rank):
    q, r = lineage(query), lineage(reference)
    for a, b in zip(q[:RANKS.index(rank)+1], r[:RANKS.index(rank)+1]):
        if a is None or b is None:
            return 'unknown'
        if a != b:
            return 'distal'
    return 'close'


def stats(values):
    values = [v for v in values if v is not None]
    return dict(numerator=sum(values), denominator=len(values),
                mean=sum(values)/len(values) if values else None,
                median=median(values) if values else None,
                maximum=max(values) if values else None,
                votes=dict(sorted(Counter(values).items())))


def decide(score, n_hits, n_genomes, query_close_known, close_genomes, cfg):
    reasons = []
    if not query_close_known:
        reasons.append('query_close_unknown')
    if close_genomes == 0:
        reasons.append('close_reference_absent')
    if n_hits < cfg['minimum_hits']:
        reasons.append('insufficient_scored_hits')
    if n_genomes < cfg['minimum_genomes']:
        reasons.append('insufficient_scored_genomes')
    if reasons:
        return 'pending', reasons
    if score >= cfg['candidate_mean']:
        return 'candidate', ['distance_threshold_met']
    return 'non_candidate', ['distance_below_threshold']


def reference_context(query, metadata, config):
    validate_lineages([query['gtdb_lineage']] + [v['gtdb_lineage'] for v in metadata.values()])
    sc = config['scoring']
    classifications = {}
    close = set()
    for sid, m in metadata.items():
        own = bool(query.get('genome_id')) and query['genome_id'] == m['genome_id']
        result = distance(query['gtdb_lineage'], m['gtdb_lineage'], sc['points'], own)
        prox = proximity(query['gtdb_lineage'], m['gtdb_lineage'], sc['close_rank'])
        classifications[sid] = (own, result, prox)
        if not own and prox == 'close':
            close.add(m['genome_id'])
    return classifications, len(close)


def hit_sort_key(h):
    return (-h['bitscore'], h['evalue'], h['sseqid'],
            -h['length'], -h['pident'], -h['qcov'], -h['scov'],
            h.get('qstart', 0), h.get('qend', 0),
            h.get('sstart', 0), h.get('send', 0), h['qlen'], h['slen'])


def apply_best_hit_genus_veto(status, reason_codes, best_hit_same_genus):
    """Exclude supported distance candidates whose selected best hit is same-genus.

    Pending evidence stays pending; unknown genus is never interpreted as distal.
    """
    if status == 'candidate' and best_hit_same_genus is True:
        return 'non_candidate', ['best_hit_same_genus']
    return status, list(reason_codes)


def score_hits(hits, query, metadata, config, context=None):
    """hits have parsed DIAMOND fields; metadata keyed by unique subject ID."""
    sc, search = config['scoring'], config['search']
    classifications, close_count = context if context is not None else reference_context(query, metadata, config)
    audit, seen, selected = [], set(), []
    # Include additional deterministic tie breakers for repeated subject/HSP rows.
    ordered = sorted(hits, key=hit_sort_key)
    for h in ordered:
        m = metadata[h['sseqid']]
        self_hit, (point, reason, common), prox = classifications[h['sseqid']]
        row = dict(h, **m, query_lineage=query['gtdb_lineage'], distance_points=point,
                   distance_reason=reason, common_rank=common,
                   proximity=prox)
        exclusion = None
        if self_hit:
            exclusion = 'self_genome'
        elif h['evalue'] > search['evalue'] or h['qcov'] < search['query_coverage'] or h['scov'] < search['subject_coverage']:
            exclusion = 'quality_filter'
        elif h['sseqid'] in seen:
            exclusion = 'duplicate_subject'
        else:
            seen.add(h['sseqid'])
            if len(selected) >= search['top_k']:
                exclusion = 'top_k_truncation'
        row['selected'] = exclusion is None
        row['exclusion_reason'] = exclusion
        audit.append(row)
        if row['selected']:
            selected.append(row)
    unique = {}
    for row in selected:
        unique.setdefault(row['genome_id'], row)
    independent = list(unique.values())
    main = stats([h['distance_points'] for h in independent])
    all_stats = stats([h['distance_points'] for h in selected])
    q = lineage(query['gtdb_lineage'])
    known = all(v is not None for v in q[:RANKS.index(sc['close_rank'])+1])
    status, reasons = decide(main['mean'], all_stats['denominator'], main['denominator'], known, close_count, sc)
    base_status, base_reasons = status, list(reasons)
    best_genus = proximity(query['gtdb_lineage'], selected[0]['gtdb_lineage'], 'genus') if selected else 'unknown'
    best_same_genus = None if best_genus == 'unknown' else best_genus == 'close'
    if sc.get('exclude_best_hit_same_genus', False):
        status, reasons = apply_best_hit_genus_veto(status, reasons, best_same_genus)
    flags = ['db_sampling_bias_unadjusted', 'paralogs_unresolved', 'no_close_presence_inference']
    n_subjects = len({h['sseqid'] for h in hits})
    if n_subjects >= search['pool_size']:
        flags.append('search_pool_limit_reached')
    truncated = [r for r in audit if r['exclusion_reason'] == 'top_k_truncation']
    if truncated:
        flags.append('top_k_truncated')
        if selected and truncated[0]['bitscore'] == selected[-1]['bitscore']:
            flags.append('boundary_bitscore_tie')
    if any(r['exclusion_reason'] == 'quality_filter' for r in audit):
        flags.append('quality_filtered')
    if not close_count:
        flags.append('close_reference_absent')
    if not query.get('genome_id'):
        flags.append('self_genome_id_unknown')
    if len(selected) == search['top_k']:
        flags.append('top_k_limit_reached')
        if len(selected) > 1 and selected[-1]['bitscore'] == selected[-2]['bitscore'] and 'boundary_bitscore_tie' not in flags:
            flags.append('boundary_bitscore_tie')
    prox = Counter(h['proximity'] for h in selected)
    comparable = prox['close'] + prox['distal']
    maxima = {k: max((h['bitscore'] for h in selected if h['proximity'] == k), default=None) for k in ('close','distal')}
    evidence = dict(status=status, reason_codes=reasons, flags=flags, raw_hits=len(hits),
        n_selected=len(selected), n_unique_genomes=len(unique),
        n_unique_species=len({tuple(lineage(h['gtdb_lineage'])) for h in selected if all(lineage(h['gtdb_lineage']))}),
        taxonomy_distance_score=main['mean'], independent=main, all_hits=all_stats,
        n_unknown=len(selected)-all_stats['denominator'], n_unknown_unique_genomes=len(independent)-main['denominator'],
        exclusions=dict(Counter(h['exclusion_reason'] for h in audit if h['exclusion_reason'])),
        db_close_genomes=close_count, proximity_counts=dict(prox),
        distal_fraction=prox['distal']/comparable if comparable else None,
        proximity_denominator=comparable, best_hit=selected[0]['sseqid'] if selected else None,
        best_bitscore=selected[0]['bitscore'] if selected else None,
        close_max_bitscore=maxima['close'], distal_max_bitscore=maxima['distal'],
        distal_minus_close=maxima['distal']-maxima['close'] if all(v is not None for v in maxima.values()) else None,
        top_group_fraction=max(prox.values())/len(selected) if selected else None)
    if sc.get('exclude_best_hit_same_genus', False):
        applied = base_status == 'candidate' and status == 'non_candidate'
        evidence.update(base_distance_status=base_status, base_distance_reason_codes=base_reasons,
                        best_hit_same_genus=best_same_genus, genus_veto_applied=applied)
        if applied:
            evidence['flags'].append('best_hit_same_genus_veto')
    counts = []
    for mode, rows in [('protein', selected), ('genome_best', independent)]:
        for i, rank in enumerate(RANKS):
            taxa = Counter(';'.join(PREFIXES[j]+lineage(h['gtdb_lineage'])[j] for j in range(i+1))
                           if all(lineage(h['gtdb_lineage'])[:i+1]) else 'unknown' for h in rows)
            for taxon, count in sorted(taxa.items()):
                counts.append(dict(mode=mode, rank=rank, taxon=taxon, count=count,
                    fraction=count/len(rows), denominator=len(rows), known_denominator=sum(v for k,v in taxa.items() if k != 'unknown')))
    return dict(evidence=evidence, selected=selected, audit=audit, taxonomy_counts=counts)
