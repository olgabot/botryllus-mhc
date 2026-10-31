#!/usr/bin/env python3
"""Size floor-rounded searches from current and historical observed peaks."""
import argparse
import csv
import math
import re
from pathlib import Path
from measure_resources import amount
from plan_floor_pairs import write_table


def table(path):
    return list(csv.DictReader(Path(path).read_text().splitlines(), delimiter='\t'))


def measurements(path, historical=False):
    result = {}
    for row in table(path):
        if ('kmerseekSearch' if historical else 'KMERSEEK_SEARCH') not in row['name']:
            continue
        if row['status'] not in ('COMPLETED', 'CACHED'):
            raise ValueError(f'Diagnose unfinished/failed task first: {row["name"]}')
        pattern = r'\.([a-z0-9_]+)\.k(\d+)\.s1\.' if historical else r'\(([a-z0-9_]+)\.k(\d+)\.'
        match = re.search(pattern, row['name'])
        if not match:
            continue
        peak = amount(row['peak_rss'], {'GB': 1, 'MB': 1/1024, 'KB': 1/1024**2, 'B': 1/1024**3})
        hours = amount(row['realtime'], {'d': 24, 'h': 1, 'm': 1/60, 's': 1/3600, 'ms': 1/3600000})
        key = (match[1], int(match[2]))
        prev = result.get(key, (0, 0))
        result[key] = max(prev[0], peak), max(prev[1], hours)
    return result


def main():
    ap = argparse.ArgumentParser()
    for name in ('pairs', 'current', 'historical', 'output'):
        ap.add_argument('--' + name, required=True)
    args = ap.parse_args()
    current, historical = measurements(args.current), measurements(args.historical, True)
    rows = []
    for pair in table(args.pairs):
        a, k = pair['alphabet'], int(pair['ksize'])
        candidates = [n for name, n in historical if name == a and n <= k]
        if not candidates:
            candidates = [n for name, n in historical if name == a]
        hk = max(candidates) if any(n <= k for n in candidates) else min(candidates)
        ck = min((n for name, n in current if name == a), key=lambda n: (abs(n-k), n))
        # A smaller k may need more memory/time: double a larger-k measurement.
        hp, ht = historical[(a, hk)]
        cp, ct = current[(a, ck)]
        hp, ht = (hp * 2, ht * 2) if hk > k else (hp, ht)
        cp, ct = (cp * 2, ct * 2) if ck > k else (cp, ct)
        memory = max(32, math.ceil(1.5 * max(hp, cp) / 8) * 8)
        hours = max(2, math.ceil(1.5 * max(ht, ct)))
        if memory > 500:
            raise ValueError(f'{a} k{k} needs {memory} GB; exceeds ceiling')
        rows.append(dict(alphabet=a, ksize=k, memory_gb=memory, time_h=hours,
                         basis=f'dark-set k{hk}, current k{ck}; larger-k measurement doubled; 1.5x margin',
                         historical_peak_gb=historical[(a, hk)][0], current_peak_gb=current[(a, ck)][0]))
    write_table(Path(args.output), rows)
    for row in rows:
        print(row['alphabet'], row['ksize'], row['memory_gb'], row['time_h'])


if __name__ == '__main__':
    main()
