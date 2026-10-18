"""Finite evidence for the fully returned diagonal translation echo theorem.

The general degree bound is analytic. These checks retain literal phases,
all physical inputs, actual inverses, dependent translation directions, and
all possible source truth tables in a bounded linear-algebra calculation.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
import unittest


def _ccz(a, b, c):
    schedule = [('T', [('T', a), ('T', b), ('T', c)])]
    for control, target in combinations((a, b, c), 2):
        compute = [('CX', control, target)]
        schedule += [('C', compute), ('T', [('TDG', target)]),
                     ('C', list(reversed(compute)))]
    compute = [('CX', a, c), ('CX', b, c)]
    return schedule + [('C', compute), ('T', [('T', c)]),
                       ('C', list(reversed(compute)))]


def _inverse(schedule):
    inverse_name = {'T': 'TDG', 'TDG': 'T', 'CX': 'CX'}
    return [(kind, [(inverse_name[g[0]], *g[1:]) for g in reversed(gates)])
            for kind, gates in reversed(schedule)]


def _source():
    return _ccz(3, 7, 8) + _ccz(4, 6, 8) + _ccz(5, 6, 7)


def _echo():
    source = _source()
    inverse = _inverse(source)
    a = [('C', [('CX', i, i + 3) for i in range(3)])]
    b = [('C', [('CX', i, i + 6) for i in range(3)])]
    # Chronological order for A (D B D†) A† (D B† D†).
    return inverse + b + source + a + inverse + b + source + a


def _column(schedule, column):
    """Exact permutation and phase in units of pi/4; no dense matrices."""
    phase = 0
    for _, gates in schedule:
        for gate in gates:
            if gate[0] == 'CX':
                if column >> gate[1] & 1:
                    column ^= 1 << gate[2]
            else:
                phase += (1 if gate[0] == 'T' else -1) * (column >> gate[1] & 1)
    return column, phase % 8


def _subsets(mask):
    sub = mask
    while sub:
        yield sub
        sub = (sub - 1) & mask


def _linear(images, x):
    value = 0
    for i, image in enumerate(images):
        if x >> i & 1:
            value ^= image
    return value


def _delta(table, shift):
    return [entry ^ table[z ^ shift] for z, entry in enumerate(table)]


def _anf(values):
    coefficients = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for x in range(len(values)):
            if x >> bit & 1:
                coefficients[x] ^= coefficients[x ^ (1 << bit)]
    return coefficients


def _additive_delta(table, shift):
    """Difference in R/Z, with exact rational phases and the literal sign."""
    return [(table[z ^ shift] - value) % 1 for z, value in enumerate(table)]


def _additive_anf(values):
    coefficients = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for x in range(len(values)):
            if x >> bit & 1:
                coefficients[x] = (coefficients[x]
                                   - coefficients[x ^ (1 << bit)]) % 1
    return coefficients


def _phase_echo_column(source, left, right, z):
    """Execute the actual D†, B, D, A, D†, B, D, A word."""
    phase = -source[z]
    z ^= right
    phase += source[z]
    z ^= left
    phase -= source[z]
    z ^= right
    phase += source[z]
    z ^= left
    return z, phase % 1


class ReturnedDiagonalEchoDegreeTests(unittest.TestCase):
    def test_literal_cubic_native_word_and_full_resource_ledger(self):
        source = _source()
        inverse = _inverse(source)
        schedule = _echo()
        layers = [gates for kind, gates in schedule if kind == 'T']
        self.assertEqual(len(layers), 60)
        self.assertEqual(sum(map(len, layers)), 84)
        self.assertEqual(sum(len(gates) for kind, gates in schedule
                             if kind == 'C'), 132)
        for layer in layers:
            targets = [gate[1] for gate in layer]
            self.assertEqual(len(targets), len(set(targets)))
            self.assertTrue(all(3 <= target < 9 for target in targets))
        for kind, gates in schedule:
            if kind == 'C':
                self.assertTrue(all(g[0] == 'CX' and 3 <= g[2] < 9
                                    for g in gates))
        for column in range(512):
            p = [(column >> (i + 3)) & 1 for i in range(3)]
            q = [(column >> (i + 6)) & 1 for i in range(3)]
            h = (p[0] * q[1] * q[2] ^ p[1] * q[0] * q[2]
                 ^ p[2] * q[0] * q[1])
            self.assertEqual(_column(source, column), (column, 4 * h))
            self.assertEqual(_column(inverse, column), (column, 4 * h))
            self.assertEqual(_column(source + inverse, column), (column, 0))
            self.assertEqual(_column(schedule, column),
                             (column, 4 * int(column & 7 == 7)))

    def test_derivative_coefficients_allow_dependent_overlapping_maps(self):
        # L has dependent directions and a zero column; M shares directions
        # with L. The source is arbitrary and is not assumed fully returned.
        left = (1, 2, 3, 0)
        right = (3, 2, 1, 4)
        source = [(0b10110101 >> z) & 1 for z in range(8)]
        fibres = [_delta(_delta(source, _linear(left, x)), _linear(right, x))
                  for x in range(16)]
        coefficients = [_anf([fibres[x][z] for x in range(16)])
                        for z in range(8)]
        for mask in range(1, 16):
            predicted = [0] * 8
            for a in _subsets(mask):
                for b in _subsets(mask):
                    if a | b != mask:
                        continue
                    term = source
                    for i in range(4):
                        if a >> i & 1:
                            term = _delta(term, left[i])
                        if b >> i & 1:
                            term = _delta(term, right[i])
                    predicted = [x ^ y for x, y in zip(predicted, term)]
            self.assertEqual(predicted, [coefficients[z][mask] for z in range(8)])

    def test_additive_derivative_sign_and_coefficients_with_overlaps(self):
        # Arbitrary Z/16 phases: the source need not return a scalar.
        # Zero, dependent and shared columns exercise the unrestricted maps.
        source = [Fraction(value, 16) for value in (0, 1, 7, 3, 11, 5, 2, 13)]
        left = (1, 2, 3, 0)
        right = (3, 2, 1, 4)
        fibres = []
        for x in range(16):
            l, r = _linear(left, x), _linear(right, x)
            difference = _additive_delta(_additive_delta(source, l), r)
            fibres.append(difference)
            for z in range(8):
                self.assertEqual(_phase_echo_column(source, l, r, z),
                                 (z, -difference[z] % 1))
        for shift in set(left + right):
            once = _additive_delta(source, shift)
            self.assertEqual(_additive_delta(once, shift),
                             [(-2 * value) % 1 for value in once])
        coefficients = [_additive_anf([fibres[x][z] for x in range(16)])
                        for z in range(8)]
        for mask in range(1, 16):
            predicted = [Fraction(0)] * 8
            for a in _subsets(mask):
                for b in _subsets(mask):
                    if a | b != mask:
                        continue
                    term = source
                    for i in range(4):
                        if a >> i & 1:
                            term = _additive_delta(term, left[i])
                        if b >> i & 1:
                            term = _additive_delta(term, right[i])
                    predicted = [(u + v) % 1 for u, v in zip(predicted, term)]
            self.assertEqual(predicted, [coefficients[z][mask] for z in range(8)])
        # Signed differences are essential for this non-Boolean source.
        self.assertNotEqual(_additive_delta(source, 1),
                            [(-value) % 1 for value in _additive_delta(source, 1)])

    def test_non_boolean_source_returns_exact_cubic_phase(self):
        # This is an identity witness, not a native compilation claim. The
        # arbitrary-denominator one-block terms cancel; a quarter-phase
        # parity straddles both blocks and contributes the linear x_1 term.
        source = []
        for z in range(64):
            p = [(z >> i) & 1 for i in range(3)]
            q = [(z >> (i + 3)) & 1 for i in range(3)]
            cubic = (p[0] * q[1] * q[2] ^ p[1] * q[0] * q[2]
                     ^ p[2] * q[0] * q[1])
            source.append((Fraction(cubic, 2) + Fraction(p[0] ^ q[0], 4)
                           + Fraction(p[1], 7) + Fraction(q[1] * q[2], 11)) % 1)
        self.assertTrue(any(2 * value % 1 for value in source))
        for x in range(8):
            expected = Fraction(int(x == 7) ^ (x & 1), 2)
            for z in range(64):
                self.assertEqual(_phase_echo_column(source, x, x << 3, z),
                                 (z, expected))
            difference = _additive_delta(_additive_delta(source, x), x << 3)
            self.assertEqual(difference, [expected] * 64)
        # Coincident directions, including zero, also force a Boolean scalar.
        quarter_parity = [Fraction(z & 1, 4) for z in range(4)]
        for shift in range(4):
            for z in range(4):
                self.assertEqual(_phase_echo_column(quarter_parity, shift, shift, z),
                                 (z, Fraction(shift & 1, 2)))

    def test_all_returned_truth_tables_have_no_quartic_coefficient(self):
        # Unknowns are all h(p,q) and q(x) entries, with no source-degree
        # ansatz. Gaussian elimination tests every exact return equation.
        for n, expected_rank in enumerate((2, 9, 46, 216), start=1):
            sources = 1 << (2 * n)
            addresses = 1 << n
            pivots = {}

            def reduce(row):
                while row:
                    pivot = row.bit_length() - 1
                    if pivot not in pivots:
                        break
                    row ^= pivots[pivot]
                return row

            for x in range(addresses):
                for z in range(sources):
                    row = ((1 << z) ^ (1 << (z ^ x))
                           ^ (1 << (z ^ (x << n)))
                           ^ (1 << (z ^ x ^ (x << n)))
                           ^ (1 << (sources + x)))
                    row = reduce(row)
                    if row:
                        pivots[row.bit_length() - 1] = row
            self.assertEqual(len(pivots), expected_rank)
            highest_possible_degree = 0
            for monomial in range(addresses):
                row = 1 << sources  # The zero-address entry belongs to every ANF sum.
                for x in _subsets(monomial):
                    row ^= 1 << (sources + x)
                if reduce(row):
                    highest_possible_degree = max(highest_possible_degree,
                                                  monomial.bit_count())
                if monomial.bit_count() >= 4:
                    self.assertEqual(reduce(row), 0)
            self.assertEqual(highest_possible_degree, min(n, 3))


if __name__ == '__main__':
    unittest.main()
