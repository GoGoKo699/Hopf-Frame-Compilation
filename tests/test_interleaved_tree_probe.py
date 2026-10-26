"""Literal pulled-back tree probe and its nondegenerate target Pauli algebra."""
from collections import Counter
import unittest

from tests.exact_matrices import CO
from tests.test_common_probe_source import _native
from tests.test_fixed_readout_syndrome_bound import _rank
from tests.test_joint_lookup_pauli_source import _bands, _inverse
from tests.test_minimal_tree_probe import _tree_forward, _tree_query


def _lookup_vector(p, x):
    return sum(1 << subset for subset in range(1 << p)
               if subset & x == subset)


def _pulled_back_label(p, x):
    m = 1 << p
    v = _lookup_vector(p, x)
    transpose_v = sum(((_tree_forward(x, 1 << j, p) & v).bit_count() & 1) << j
                      for j in range(m))
    return v | (_tree_forward(x, v, p) << m), transpose_v | (v << m)


def _pair(left, right):
    return ((left[0] & right[1]).bit_count()
            + (left[1] & right[0]).bit_count()) & 1


def _kernel(p, x, y):
    value = 0
    for i in range(p):
        if (x ^ y) >> i & 1 and not ((x & y) >> (i + 1)):
            value ^= 1
    return value


def _gram_rows(p):
    labels = [_pulled_back_label(p, x) for x in range(1 << p)]
    return [sum(_pair(left, right) << y for y, right in enumerate(labels))
            for left in labels]


def _interleaved_word(p):
    word, width = _tree_query(p)
    m = 1 << p
    return (word + [('H', p + j) for j in range(m)]
            + [('CX', p + j, p + m + j) for j in range(m)]), width


def _canonical_pairs(labels):
    remaining, pairs = list(labels), []
    while remaining:
        a = remaining.pop()
        index = next(i for i, b in enumerate(remaining) if _pair(a, b))
        b = remaining.pop(index)
        updated = []
        for value in remaining:
            use_a, use_b = _pair(value, b), _pair(value, a)
            updated.append((value[0] ^ (a[0] if use_a else 0)
                            ^ (b[0] if use_b else 0),
                            value[1] ^ (a[1] if use_a else 0)
                            ^ (b[1] if use_b else 0)))
        remaining = updated
        pairs.append((a, b))
    return pairs


class InterleavedTreeProbeTests(unittest.TestCase):
    def test_complete_native_probe_literal_sign_and_arbitrary_helpers(self):
        for p in (1, 2):
            m = 1 << p
            word, width = _interleaved_word(p)
            reverse = _inverse(word)
            if p == 1:
                inputs = range(1 << width)
            else:
                # Every address and helper column, with varied nonzero old
                # carrier/output values; p=1 covers every full basis column.
                inputs = (x | (((3 * h + x) & 15) << p)
                          | (((11 * h + 2 * x + 5) & 15) << (p + m))
                          | (h << (p + 2 * m))
                          for x in range(m) for h in range(1 << (m - 1)))
            for basis in inputs:
                x = basis & (m - 1)
                a, b = _pulled_back_label(p, x)
                state = _native(word, {basis: CO})
                phase = [('S', p + m + subset) for subset in range(m)
                         if subset & x == subset for _ in range(2)]
                actual = _native(reverse, _native(phase, state))
                sign = -1 if ((basis >> p) & b).bit_count() & 1 else 1
                self.assertEqual(actual, {basis ^ (a << p): sign * CO})
                self.assertEqual(_native(reverse, state), {basis: CO})

    def test_native_wall_and_actual_inverse_have_complete_common_ledger(self):
        for p in range(1, 8):
            m = 1 << p
            word, width = _interleaved_word(p)
            for native in (word, _inverse(word)):
                counts = Counter(gate[0] for gate in native)
                t = counts['T'] + counts['TDG']
                self.assertEqual(t, 8 * (m - 1))
                self.assertEqual(len(native) - t, 12 * m - 8)
                self.assertEqual(len(native), 20 * m - 16)
                self.assertLessEqual(len(_bands(native)), 8)
                self.assertEqual(width, p + 3 * m - 1)

    def test_labels_are_independent_hermitian_and_keep_root_minus_sign(self):
        for p in range(1, 8):
            m = 1 << p
            labels = [_pulled_back_label(p, x) for x in range(m)]
            self.assertEqual(_rank(a | (b << (2 * m)) for a, b in labels), m)
            self.assertEqual(_rank(a & ((1 << m) - 1) for a, _ in labels), m)
            for a, b in labels:
                self.assertEqual(b & ((1 << m) - 1), 1)
                self.assertEqual((a & b).bit_count() % 2, 0)
            a, b = labels[0]
            self.assertEqual(a, 1 | (1 << m))
            self.assertEqual(a, b)
            self.assertEqual((-1) ** ((a & b).bit_count() // 2), -1)

    def test_direct_label_gram_matches_kernel_block_recursion_and_full_rank(self):
        previous = None
        for p in range(1, 8):
            m = 1 << p
            rows = _gram_rows(p)
            self.assertEqual(_rank(rows), m)
            for x, row in enumerate(rows):
                self.assertEqual(row, sum(_kernel(p, x, y) << y for y in range(m)))
                self.assertEqual((row >> x) & 1, 0)
            if previous is not None:
                half, all_ones = m // 2, (1 << (m // 2)) - 1
                self.assertEqual(rows[:half],
                                 [row | ((row ^ all_ones) << half) for row in previous])
                self.assertEqual(rows[half:], [row ^ all_ones for row in previous])
                self.assertEqual(_rank(row ^ all_ones for row in previous), half)
            previous = rows

    def test_symplectic_reduction_gives_exact_target_logical_pairs(self):
        for p in range(1, 7):
            m = 1 << p
            labels = [_pulled_back_label(p, x) for x in range(m)]
            pairs = _canonical_pairs(labels)
            self.assertEqual(len(pairs), m // 2)
            flat = [label for pair in pairs for label in pair]
            self.assertEqual(_rank(a | (b << (2 * m)) for a, b in flat), m)
            for i, left in enumerate(flat):
                for j, right in enumerate(flat):
                    self.assertEqual(_pair(left, right), int(i ^ j == 1))


if __name__ == '__main__':
    unittest.main()
