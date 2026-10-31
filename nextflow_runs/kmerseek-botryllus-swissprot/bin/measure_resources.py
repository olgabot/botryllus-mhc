#!/usr/bin/env python3
"""Raise baseline search requests using saved traces; never rewrite tracked assets."""
import argparse
import csv
import math
import re
from pathlib import Path


def amount(text, units):
    parts = re.findall(r'([\d.]+)\s*([A-Za-z]+)', text)
    if not parts:
        raise ValueError(f'Missing measurement: {text}')
    return sum(float(n) * units[u] for n, u in parts)


def measure(baseline, history):
    rows = list(csv.DictReader(Path(baseline).read_text().splitlines(), delimiter='\t'))
    by_pair = {(r['alphabet'], r['ksize']): r for r in rows}
    for trace in sorted(Path(history).glob('trace-*.tsv')):
        for task in csv.DictReader(trace.read_text().splitlines(), delimiter='\t'):
            if not task['name'].startswith('KMERSEEK_SEARCH'):
                continue
            alphabet, k, _ = task['tag'].split('.')
            row = by_pair[(alphabet, k[1:])]
            mem = amount(task['memory'], {'GB': 1, 'MB': 1/1024, 'TB': 1024})
            hours = amount(task['time'], {'h': 1, 'm': 1/60, 's': 1/3600, 'd': 24})
            need_mem, need_time = float(row['memory_gb']), float(row['time_h'])
            if task['status'] in ('COMPLETED', 'CACHED'):
                if task['peak_rss'] != '-':
                    peak = amount(task['peak_rss'], {'GB': 1, 'MB': 1/1024, 'KB': 1/1024**2, 'B': 1/1024**3})
                    need_mem = max(need_mem, 1.25 * peak)
                if task['realtime'] != '-':
                    duration = amount(task['realtime'], {'h': 1, 'm': 1/60, 's': 1/3600, 'ms': 1/3600000, 'd': 24})
                    need_time = max(need_time, 1.5 * duration)
            elif task['status'] == 'FAILED':
                # An unexplained failure must be diagnosed before a new submission.
                err = Path(task['workdir']) / '.command.err'
                killed = task['exit'] in ('137', '135') or (err.exists() and 'search 137,' in err.read_text())
                if killed:
                    need_mem = max(need_mem, 2 * mem)
                elif task['exit'] in ('140', '-', '2147483647'):
                    need_time = max(need_time, 2 * hours)
                else:
                    raise ValueError(f'Unexplained failure: {task["tag"]}, {err}')
            else:
                raise ValueError(f'Unfinished task in saved trace: {task["tag"]}')
            row['memory_gb'] = str(math.ceil(need_mem / 8) * 8)
            row['time_h'] = str(math.ceil(need_time))
            if int(row['memory_gb']) > 500:
                raise ValueError(f'Request exceeds 500 GB: {row}')
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--baseline', required=True)
    ap.add_argument('--history', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    rows = measure(args.baseline, args.history)
    with open(args.output, 'w') as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row['alphabet'], row['ksize'], row['memory_gb'] + ' GB', row['time_h'] + ' h')


if __name__ == '__main__':
    main()
