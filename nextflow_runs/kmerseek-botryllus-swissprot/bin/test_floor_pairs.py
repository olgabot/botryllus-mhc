import math
import unittest
from plan_floor_pairs import build_rows
from plan_pairs import CLUSTERS, swissprot_bits


class FloorPairsTests(unittest.TestCase):
    def test_floor_rule(self):
        rows = build_rows()
        self.assertEqual(len(rows), 38)
        self.assertEqual(len({(r['alphabet'], r['ksize']) for r in rows}), 38)
        for r in rows:
            target = r['reference_k'] * swissprot_bits(CLUSTERS['hp_lehninger2'])
            bits = swissprot_bits(CLUSTERS[r['alphabet']])
            self.assertLessEqual(r['ksize']*bits, target + 1e-9)
            self.assertGreater((r['ksize']+1)*bits, target)
        self.assertEqual([r['ksize'] for r in rows if r['alphabet']=='protein20'], [5, 6])
        self.assertEqual([r['ksize'] for r in rows if r['alphabet']=='hp_lehninger2'], [24, 27])


if __name__ == '__main__':
    unittest.main()
