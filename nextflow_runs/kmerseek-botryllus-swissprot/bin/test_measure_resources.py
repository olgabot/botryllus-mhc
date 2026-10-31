"""Exercise the wrapped memory kill and refusal to relaunch unexplained failures."""
import csv
import tempfile
import unittest
from pathlib import Path
from measure_resources import amount, measure


class ResourceTests(unittest.TestCase):
    def test_units(self):
        self.assertAlmostEqual(amount('1h 30m', {'h': 1, 'm': 1/60}), 1.5)

    def test_memory_kill_and_unknown_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            baseline = root / 'baseline.tsv'
            baseline.write_text('alphabet\tksize\tmemory_gb\ttime_h\nx\t26\t32\t2\n')
            task = dict(name='KMERSEEK_SEARCH (x.k26.chunk_0000)', tag='x.k26.chunk_0000',
                        status='FAILED', exit='1', memory='32 GB', time='2h',
                        workdir=str(root))
            with (root / 'trace-test.tsv').open('w') as out:
                writer = csv.DictWriter(out, fieldnames=list(task), delimiter='\t')
                writer.writeheader()
                writer.writerow(task)
            (root / '.command.err').write_text('search 137, filter 0, zstd 0\n')
            self.assertEqual(measure(baseline, root)[0]['memory_gb'], '64')
            (root / '.command.err').write_text('wrong column\n')
            with self.assertRaisesRegex(ValueError, 'Unexplained failure'):
                measure(baseline, root)


if __name__ == '__main__':
    unittest.main()
