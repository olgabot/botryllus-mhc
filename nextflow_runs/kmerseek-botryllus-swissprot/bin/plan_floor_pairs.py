#!/usr/bin/env python3
"""Add floor-rounded pairs anchored to hp_lehninger2 k=24 and k=27 (31 October)."""
import csv
import math
from pathlib import Path
from plan_pairs import CLUSTERS, EXTENSION, swissprot_bits

ROOT = Path(__file__).resolve().parent.parent


def build_rows():
    bits = {a: swissprot_bits(c) for a, c in CLUSTERS.items()}
    ranges = {r['alphabet']: r for r in csv.DictReader(
        (ROOT / 'assets/kmerseek_ksize_range_per_alphabet.tsv').read_text().splitlines(), delimiter='\t')}
    rows = []
    for ref in (24, 27):
        target = ref * bits['hp_lehninger2']
        for a, classes in CLUSTERS.items():
            k = math.floor(target / bits[a] + 1e-9)
            assert int(ranges[a]['k_min']) <= k <= int(ranges[a]['k_max'])
            assert k * bits[a] <= target + 1e-9 < (k + 1) * bits[a]
            rows.append(dict(alphabet=a, n_classes=len(classes), bits_per_position=f'{bits[a]:.8f}',
                             protein20_k=math.floor(target / bits['protein20']), ksize=k,
                             bits_per_kmer=f'{k * bits[a]:.8f}', penalty=EXTENSION[a][0], xdrop=EXTENSION[a][1],
                             reference_alphabet='hp_lehninger2', reference_k=ref,
                             reference_bits=f'{target:.8f}', rounding='floor'))
    return rows


def write_table(path, rows):
    with path.open('w') as out:
        w = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def main():
    rows = build_rows()
    old = {(r['alphabet'], int(r['ksize'])) for r in csv.DictReader(
        (ROOT / 'assets/pairs.tsv').read_text().splitlines(), delimiter='\t')}
    reused = [r for r in rows if (r['alphabet'], r['ksize']) in old]
    assert len(rows) == 38 and len(reused) == 8
    write_table(ROOT / 'assets/pairs_floor_lehninger.tsv', rows)
    write_table(ROOT / 'assets/reuse_floor_lehninger.tsv', reused)
    print('38 pairs: 30 new, 8 reused; 690 new searches, 184 reused')


if __name__ == '__main__':
    main()
