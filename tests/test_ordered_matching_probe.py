"""Exact ordered-matching probe algebra and returned native components.

The all-size source has a paid linear-depth baseline. These finite checks
do not emit a sublogarithmic lookup source or establish a depth lower bound.
"""
from collections import Counter
from itertools import combinations
import unittest

from tests.exact_matrices import CI, CO
from tests.test_coherent_feedback_sources import _shared_ccz
from tests.test_common_probe_source import ROOT_HALF, _native
from tests.test_dirty_cubic_phase_compiler import _gray_phase_batch
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _pairs(p):
    return tuple(combinations(range(p), 2))


def _matching(subset, p):
    labels = [i for i in range(p) if subset >> i & 1]
    if len(labels) % 2:
        labels = labels[1:]
    return tuple(zip(labels[::2], labels[1::2]))


def _occurrences(p, family=None):
    pairs = _pairs(p)
    allowed = set(pairs if family is None else family)
    return [(pairs.index(edge), subset, subset ^ (1 << edge[0]) ^ (1 << edge[1]))
            for subset in range(1 << p)
            for edge in _matching(subset, p) if edge in allowed]


def _rows(p, key):
    rows = [1 << subset for subset in range(1 << p)]
    for index, edge in enumerate(_pairs(p)):
        if key >> index & 1:
            for _, source, target in _occurrences(p, (edge,)):
                rows[target] ^= rows[source]
    return tuple(rows)


def _row_times(row, matrix):
    result = 0
    while row:
        bit = row & -row
        result ^= matrix[bit.bit_length() - 1]
        row ^= bit
    return result


def _product(left, right):
    return tuple(_row_times(row, right) for row in left)


def _apply(rows, carrier):
    return sum(((row & carrier).bit_count() & 1) << i for i, row in enumerate(rows))


def _clique(p, address):
    return sum((((address >> i) & (address >> j)) & 1) << k
               for k, (i, j) in enumerate(_pairs(p)))


def _lookup_row(p, address):
    return sum(1 << subset for subset in range(1 << p) if subset & address == subset)


def _linear_row(p, address):
    return 1 | sum(1 << (1 << i) for i in range(p) if address >> i & 1)


def _families(p):
    # Zero-based odd labels are the even integers in the proof.
    return tuple(tuple(edge for edge in _pairs(p)
                       if next(k for k in range(edge[0], edge[1] + 1) if k % 2) == pivot)
                 for pivot in range(1, p, 2))


def _cut_word(p, family, address_wires=0):
    q, m = len(_pairs(p)), 1 << p
    carrier, private = address_wires + q, address_wires + q + m
    occurrences = _occurrences(p, family)
    terms = [(address_wires + key, carrier + source, carrier + target)
             for key, source, target in occurrences]
    targets = sorted({term[2] for term in terms})
    wall = [('H', target) for target in targets]
    return (wall + _gray_phase_batch(terms, list(range(private, private + len(terms))),
                                     [1] * len(terms), 3) + list(reversed(wall)))


def _source_word(p):
    pairs, m = _pairs(p), 1 << p
    q = len(pairs)
    carrier, private = p + q, p + q + m
    root_wall = [('H', carrier)]
    linear = root_wall + _shared_ccz(carrier, [(i, carrier + (1 << i)) for i in range(p)])
    linear += root_wall
    gauge = [gate for family in _families(p) for gate in _cut_word(p, family, p)]
    key_wall = [('H', p + k) for k in range(q)]
    toggle = key_wall + _gray_phase_batch(
        [(i, j, p + k) for k, (i, j) in enumerate(pairs)],
        list(range(private, private + q)), [1] * q, 3) + list(reversed(key_wall))
    width = private + max(q, max(map(lambda family: len(_occurrences(p, family)), _families(p))))
    return linear + gauge + toggle + _inverse(gauge), width


class OrderedMatchingProbeTests(unittest.TestCase):
    def test_commuting_generator_algebra_and_every_p4_pair_key(self):
        for p in (2, 3, 4, 6):
            identity = tuple(1 << i for i in range(1 << p))
            generators = [_rows(p, 1 << k) for k in range(len(_pairs(p)))]
            for left in generators:
                self.assertEqual(_product(left, left), identity)
                for right in generators:
                    self.assertEqual(_product(left, right), _product(right, left))
            self.assertEqual(len(_occurrences(p)), (p - 1) * (1 << (p - 2)))
        matrices = [_rows(4, key) for key in range(64)]
        for s in range(64):
            for t in range(64):
                self.assertEqual(_product(matrices[s], matrices[t]), matrices[s ^ t])

    def test_complete_probe_rows_with_arbitrary_pair_keys(self):
        for p in (2, 3, 4, 6, 8):
            q, m = len(_pairs(p)), 1 << p
            keys = range(1 << q) if p <= 4 else (0, 1, (1 << q) - 1, (1 << q) // 3)
            addresses = range(1 << p) if p <= 4 else (0, 1, 3, (1 << p) - 1, (1 << p) // 3)
            for x in addresses:
                k = _clique(p, x)
                self.assertEqual(_row_times(_linear_row(p, x), _rows(p, k)), _lookup_row(p, x))
                for s in keys:
                    # The physical gauge echo uses every input key; its carrier
                    # action is independent of that key before the probe readout.
                    self.assertEqual(_product(_rows(p, s ^ k), _rows(p, s)), _rows(p, k))
                    # This compares every carrier column, not only a zero bank.
                    row = _row_times(_linear_row(p, x),
                                     _product(_rows(p, s ^ k), _rows(p, s)))
                    self.assertEqual(row, _lookup_row(p, x))
            self.assertEqual(len(_rows(p, 0)), m)

    def test_selected_matching_differs_from_pfaffian_and_source_order_matters(self):
        p, pairs = 4, _pairs(4)
        key = (1 << pairs.index((0, 2))) | (1 << pairs.index((1, 3)))
        self.assertEqual((_rows(p, key)[0] >> 15) & 1, 0)
        # Exterior multiplication has its (13)(24) matching term here.
        self.assertEqual(((key >> pairs.index((0, 2))) &
                          (key >> pairs.index((1, 3)))) & 1, 1)
        identity = tuple(1 << i for i in range(16))
        ordered = [tuple(row ^ identity[i] for i, row in enumerate(_rows(p, 1 << pairs.index(edge))))
                   for edge in ((0, 2), (1, 3))]
        self.assertEqual(_product(*ordered), (0,) * 16)
        exterior = [tuple((1 << (subset | edge)) if not subset & edge else 0
                          for subset in range(16)) for edge in (5, 10)]
        self.assertEqual(_product(*exterior)[0], 1 << 15)
        p, x = 3, 7
        linear = [1 << i for i in range(1 << p)]
        linear[0] = _linear_row(p, x)
        correct = _product(linear, _rows(p, _clique(p, x)))
        wrong = _product(_rows(p, _clique(p, x)), linear)
        self.assertEqual(correct[0], _lookup_row(p, x))
        self.assertEqual((correct[0] >> 7) & 1, 1)
        self.assertEqual((wrong[0] >> 7) & 1, 0)

    def test_interval_cut_batches_are_native_shears_with_returned_dirty_work(self):
        for p in (2, 3, 4, 6):
            total = 0
            for family in _families(p):
                occurrences = _occurrences(p, family)
                sources, targets = ({item[column] for item in occurrences} for column in (1, 2))
                self.assertFalse(sources & targets)
                word = _cut_word(p, family)
                total += len(occurrences)
                for native in (word, _inverse(word)):
                    counts = Counter(gate[0] for gate in native)
                    self.assertEqual(counts['T'] + counts['TDG'], 8 * len(occurrences))
                    self.assertEqual(len(native), 16 * len(occurrences) + 2 * len(targets))
                    self.assertLessEqual(len(_bands(native)), 8)
            self.assertEqual(total, (p - 1) * (1 << (p - 2)))
        p, q, m = 3, 3, 8
        family = _families(p)[0]
        word, private = _cut_word(p, family), q + m
        for key in range(8):
            for carrier, dirty in ((0, 0), (1, 15), (0x5a, 6), (0xff, 9)):
                basis = key | (carrier << q) | (dirty << private)
                expected = key | (_apply(_rows(p, key), carrier) << q) | (dirty << private)
                actual = _native(word, {basis: CO})
                self.assertEqual(actual, {expected: CO})
                self.assertEqual(_native(_inverse(word), actual), {basis: CO})
        # Reference coherence is retained: no measurement or classical branch.
        state = {3 | (0x37 << q): ROOT_HALF,
                 6 | (0xa2 << q) | (13 << private): CI * ROOT_HALF}
        self.assertEqual(_native(_inverse(word), _native(word, state)), state)

    def test_complete_p2_native_source_all_columns_and_actual_probe(self):
        p, q, m = 2, 1, 4
        word, width = _source_word(p)
        carrier_start = p + q
        reverse = _inverse(word)
        self.assertEqual(width, 8)
        for basis in range(1 << width):
            x, s = basis & 3, (basis >> p) & 1
            c = (basis >> carrier_start) & 15
            linear_c = c ^ (((c & _linear_row(p, x)).bit_count() & 1) ^ (c & 1))
            expected_c = _apply(_rows(p, _clique(p, x)), linear_c)
            expected = (basis & ~(((1 << m) - 1) << carrier_start))
            expected ^= _clique(p, x) << p
            expected |= expected_c << carrier_start
            actual = _native(word, {basis: CO})
            self.assertEqual(actual, {expected: CO})
            self.assertEqual(_native(reverse, actual), {basis: CO})
            inverse = _native(reverse, {basis: CO})
            phased = {label: amplitude * ((-1) ** ((label >> carrier_start) & 1))
                      for label, amplitude in inverse.items()}
            probe = _native(word, phased)
            sign = (-1) ** ((c & _lookup_row(p, x)).bit_count() & 1)
            self.assertEqual(probe, {basis: CO * sign})
            self.assertIn(s, (0, 1))  # Both arbitrary key sectors were included.
        for p in (2, 4, 6):
            word, width = _source_word(p)
            q, m, h = len(_pairs(p)), 1 << p, p // 2
            edges, parity = (p - 1) * (1 << (p - 2)), p % 2
            for native in (word, _inverse(word)):
                counts = Counter(gate[0] for gate in native)
                self.assertLessEqual(counts['T'] + counts['TDG'], 16 * edges + 8 * q + 6 * p + parity)
                self.assertLessEqual(len(_bands(native)), 16 * h + 12)
                self.assertLessEqual(len(native), 32 * edges + 4 * m * h + 18 * q + 16 * p + parity + 5)
            self.assertLessEqual(width - p - m, q + max(edges, q))


if __name__ == '__main__':
    unittest.main()
