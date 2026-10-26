"""Sparse tree probe, its literal source, and its complete forward XOR query.

This is a linear-count source with a linear-depth native baseline. The
constant-depth out-of-place query does not implement its inverse in place.
"""
from collections import Counter
import unittest

from tests.exact_matrices import CO
from tests.test_coherent_feedback_sources import _shared_ccz
from tests.test_common_probe_source import _native
from tests.test_dirty_cubic_phase_compiler import _gray_phase_batch
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _edges(p):
    return tuple(((subset & -subset).bit_length() - 1, subset,
                  subset ^ (subset & -subset)) for subset in range(1, 1 << p))


def _tree_step(x, c, p):
    result = 0
    for i, child, parent in _edges(p):
        result ^= (((x >> i) & (c >> child) & 1) << parent)
    return result


def _tree_forward(x, c, p):
    return c ^ _tree_step(x, c, p)


def _tree_inverse(x, c, p):
    result, term = c, c
    for _ in range(p):
        term = _tree_step(x, term, p)
        result ^= term
    assert _tree_step(x, term, p) == 0
    return result


def _matching_step(x, c, p, i):
    for label, child, parent in _edges(p):
        if label == i:
            c ^= (((x >> i) & (c >> child) & 1) << parent)
    return c


def _tree_source(p):
    word = []
    for i in reversed(range(p)):
        pairs = [(p + child, p + parent) for label, child, parent in _edges(p)
                 if label == i]
        wall = [('H', target) for _, target in pairs]
        word += wall + _shared_ccz(i, pairs) + list(reversed(wall))
    return word


def _tree_query(p):
    m = 1 << p
    triples = [(i, p + child, p + m + parent) for i, child, parent in _edges(p)]
    helpers = tuple(range(p + 2 * m, p + 3 * m - 1))
    wall = [('H', target) for target in sorted({triple[-1] for triple in triples})]
    word = wall + _gray_phase_batch(triples, helpers, [1] * (m - 1), 3)
    word += list(reversed(wall)) + [('CX', p + i, p + m + i) for i in range(m)]
    return word, p + 3 * m - 1


class MinimalTreeProbeTests(unittest.TestCase):
    def test_unique_paths_source_direction_and_reversed_matching_order(self):
        for p in range(1, 6):
            m = 1 << p
            carriers = range(1 << m) if p <= 3 else (0, (1 << m) - 1, *(1 << i for i in range(m)))
            for x in range(m):
                for c in carriers:
                    forward = _tree_forward(x, c, p)
                    actual = c
                    for i in reversed(range(p)):
                        actual = _matching_step(x, actual, p, i)
                    self.assertEqual(actual, forward)
                    inverse = _tree_inverse(x, c, p)
                    self.assertEqual(_tree_forward(x, inverse, p), c)
                    self.assertEqual(_tree_inverse(x, forward, p), c)
                    actual = c
                    for i in range(p):
                        actual = _matching_step(x, actual, p, i)
                    self.assertEqual(actual, inverse)
                    target = sum((c >> subset) & 1 for subset in range(m)
                                 if subset & x == subset) & 1
                    self.assertEqual(inverse & 1, target)
        # A forward/inverse swap already loses the first long path at p=2.
        self.assertEqual(_tree_forward(3, 1 << 3, 2) & 1, 0)
        self.assertEqual(_tree_inverse(3, 1 << 3, 2) & 1, 1)

    def test_literal_native_source_and_actual_inverse_have_common_ledger(self):
        for p in (1, 2, 3):
            m = 1 << p
            source, reverse = _tree_source(p), _inverse(_tree_source(p))
            carriers = range(1 << m) if p <= 2 else (0, 1, 128, 85, 170, 255)
            for x in range(m):
                for c in carriers:
                    basis = x | (c << p)
                    target = x | (_tree_forward(x, c, p) << p)
                    self.assertEqual(_native(source, {basis: CO}), {target: CO})
                    self.assertEqual(_native(reverse, {target: CO}), {basis: CO})
        for p in range(1, 9):
            m = 1 << p
            for native in (_tree_source(p), _inverse(_tree_source(p))):
                counts = Counter(gate[0] for gate in native)
                t = counts['T'] + counts['TDG']
                self.assertEqual(t, 6 * m - 5)
                self.assertLessEqual(len(native) - t, 12 * (m - 1) + 5)
                self.assertLessEqual(len(native), 18 * m - 12)
                self.assertLessEqual(len(_bands(native)), 4 * p)
                self.assertLess(max(wire for gate in native for wire in gate[1:]), p + m)

    def test_complete_forward_query_returns_arbitrary_helpers_literally(self):
        for p in (1, 2):
            m = 1 << p
            word, width = _tree_query(p)
            reverse = _inverse(word)
            if p == 1:
                inputs = range(1 << width)
            else:
                # Cover all addresses, carriers and private helper values,
                # with arbitrary old outputs depending on the other registers.
                inputs = (x | (c << p) | (((c ^ h ^ 10) & 15) << (p + m))
                          | (h << (p + 2 * m))
                          for x in range(m) for c in range(1 << m)
                          for h in range(1 << (m - 1)))
            for basis in inputs:
                x, c = basis & (m - 1), (basis >> p) & ((1 << m) - 1)
                target = basis ^ (_tree_forward(x, c, p) << (p + m))
                self.assertEqual(_native(word, {basis: CO}), {target: CO})
                self.assertEqual(_native(reverse, {target: CO}), {basis: CO})
        for p in range(1, 7):
            m = 1 << p
            word, width = _tree_query(p)
            for native in (word, _inverse(word)):
                counts = Counter(gate[0] for gate in native)
                self.assertEqual(counts['T'] + counts['TDG'], 8 * (m - 1))
                self.assertEqual(len(native), 18 * m - 16)
                self.assertLessEqual(len(_bands(native)), 8)
                self.assertEqual(width - (p + 2 * m), m - 1)


if __name__ == '__main__':
    unittest.main()
