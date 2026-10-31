#!/usr/bin/env python3
"""Check pair counts against per-query rows and search-task coverage."""
import argparse
import collections
import csv
import gzip
from pathlib import Path


def validate(root):
    root = Path(root)
    stem = root / 'summary/kmerseek_botryllus_swissprot'
    pairs = list(csv.DictReader(Path(str(stem) + '.pairs.tsv').read_text().splitlines(), delimiter='\t'))
    assert len(pairs) == 38
    counts, strong, names = collections.Counter(), collections.Counter(), set()
    with gzip.open(str(stem) + '.best_per_protein_per_pair.tsv.gz', 'rt') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            key = (r['alphabet'], r['ksize'])
            identity = (*key, r['query_name'])
            assert identity not in names
            names.add(identity)
            counts[key] += 1
            strong[key] += int(r['n_targets_strong']) > 0
    for p in pairs:
        key = (p['alphabet'], p['ksize'])
        assert int(p['n_botryllus_proteins']) == 45339
        assert counts[key] == int(p['n_with_match'])
        assert strong[key] == int(p['n_with_strong_match'])
        assert abs(strong[key]/45339 - float(p['share_with_strong_match'])) <= .000051
    with gzip.open(str(stem) + '.per_protein.tsv.gz', 'rt') as fh:
        per_protein = list(csv.DictReader(fh, delimiter='\t'))
    assert len(per_protein) == len({r['query_name'] for r in per_protein}) == 45339
    trace = list(csv.DictReader((root/'trace.txt').read_text().splitlines(), delimiter='\t'))
    search = [r for r in trace if r['name'].startswith('KMERSEEK_SEARCH')]
    assert len(search) == len({r['tag'] for r in search}) == 874
    assert all(r['status'] in ('COMPLETED', 'CACHED') and r['exit']=='0' for r in search)
    assert len({r['tag'].rsplit('.',1)[0] for r in search}) == 38
    assert all(v==23 for v in collections.Counter(r['tag'].rsplit('.',1)[0] for r in search).values())
    print('PASS: 874 unique successful searches, 38 pairs, 23 chunks per pair, 45,339 unique queries; all pair counts agree with per-query rows.')


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('root');args=ap.parse_args();validate(args.root)
