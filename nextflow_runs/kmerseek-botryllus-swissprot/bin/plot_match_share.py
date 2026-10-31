#!/usr/bin/env python3
"""Compare complete E-value results at the two protein20 information thresholds."""
import argparse
import csv
from pathlib import Path
import pubfig as pf


def reference(r):
    return r.get('reference_k', r['protein20_k'])


def read_rows(path):
    rows = list(csv.DictReader(Path(path).read_text().splitlines(), delimiter='\t'))
    keys = {(r['alphabet'], reference(r)) for r in rows}
    if len(rows) != 38 or len(keys) != 38 or len({a for a, _ in keys}) != 19:
        raise ValueError('Expected 38 distinct pairs in 19 alphabets')
    if len({reference(r) for r in rows}) != 2:
        raise ValueError('Expected two reference sizes')
    for r in rows:
        if int(r['n_botryllus_proteins']) != 45339:
            raise ValueError('Unexpected threshold or query count')
        if r['fit'] not in ('fitted', 'unfitted'):
            raise ValueError('Missing fit status')
        if r['fit'] == 'fitted' and r['statistic'] != 'region_evalue':
            raise ValueError('Fitted pair must use E-values')
        if not 0 <= int(r['n_with_strong_match']) <= int(r['n_with_match']) <= 45339:
            raise ValueError('Invalid match count')
    return rows


def plot(rows, output):
    # What: alphabet and threshold are categories, matched share is a percentage.
    # Why: compare both thresholds for every alphabet on the same scale.
    # How: horizontally aligned points; colour and shape identify the threshold.
    pf.use_style()
    fig, ax = pf.figure(pf.TWO_COLUMN_MM, 115)
    order = list(dict.fromkeys(r['alphabet'] for r in rows))
    lookup = {(r['alphabet'], reference(r)): r for r in rows}
    thresholds = sorted({reference(r) for r in rows}, key=int)
    anchor = rows[0].get('reference_alphabet', 'protein20')
    rounding = ' (floor)' if rows[0].get('rounding') == 'floor' else ''
    for threshold, color, marker, offset in (
        (thresholds[0], pf.OKABE_ITO['blue'], 'o', -0.13),
        (thresholds[1], pf.OKABE_ITO['vermillion'], 's', 0.13),
    ):
        selected = [(i, lookup[(a, threshold)]) for i, a in enumerate(order)]
        fitted = [(i, r) for i, r in selected if r['fit'] == 'fitted']
        ax.scatter([100 * int(r['n_with_strong_match']) / 45339 for _, r in fitted],
                   [i + offset for i, _ in fitted], c=color, marker=marker, s=12,
                   label=f'Bits of {anchor} k = {threshold}{rounding}')
    for i, alphabet in enumerate(order):
        missing = [lookup[(alphabet, k)]['ksize'] for k in thresholds
                   if lookup[(alphabet, k)]['fit'] == 'unfitted']
        if missing:
            ax.text(102, i, 'No E-value fit: k = ' + ', '.join(missing), va='center')
    labels = [f'{a}  (k = {lookup[(a, thresholds[0])]["ksize"]}, {lookup[(a, thresholds[1])]["ksize"]})' for a in order]
    ax.set_yticks(range(19), labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel('Botryllus proteins with a Swiss-Prot match at E ≤ 0.001 (%)')
    ax.set_ylabel('Alphabet (k at the two bit thresholds)')
    ax.grid(axis='y', color='#dddddd', linewidth=0.3)
    ax.set_axisbelow(True)
    pf.shared_legend(fig, ncol=2)
    # Checklist: separate markers prevent overlap; legend above marks; full 0–100% scale;
    # missing fits are words outside the scale; denominator and k are explicit.
    pf.save(fig, output)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='results/summary/kmerseek_botryllus_swissprot.pairs.tsv')
    ap.add_argument('--output', default='figures/kmerseek_botryllus_swissprot_match_share')
    args = ap.parse_args()
    plot(read_rows(args.input), args.output)


if __name__ == '__main__':
    main()
