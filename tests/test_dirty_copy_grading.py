"""Exact copy grading, spectral obstructions and noncommuting examples.

The all-word criterion and all-size proofs live in DIRTY_COPY_GRADING.md.
These small full-column fixtures do not compile the field-coordinate
encoding or replace its missing native seed/copy implementation.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, Q2, Gate, add, adjoint, apply_word, eye, multiply,
    native_columns, scale, zeros,
)


def _inverse(word):
    """Actual Clifford inverse, reducing adjacent S powers literally."""
    expanded = []
    for gate in reversed(word):
        expanded.extend([gate] * (3 if gate[0] == 'S' else 1))
    reduced = []
    index = 0
    while index < len(expanded):
        gate = expanded[index]
        end = index + 1
        if gate[0] == 'S':
            while end < len(expanded) and expanded[end] == gate:
                end += 1
            reduced.extend([gate] * ((end - index) % 4))
        else:
            reduced.append(gate)
        index = end
    return reduced


def _act(word, basis):
    native = [Gate(name, tuple(wires)) for name, *wires in word]
    state = apply_word(native, {basis: CO})
    return {row: value for row, value in state.items() if value != CZ}


def _native_x(wire):
    return [('H', wire), ('S', wire), ('S', wire), ('H', wire)]


def _weighted_flip(x, a=1, b=2):
    # R_x = K (X-Y)/sqrt(2) K^-1.  Middle chronology is H,S^-1,H,S,H.
    k = [('CX', x ^ a, x), ('CX', x ^ b, x)]
    middle = ([('H', x)] + [('S', x)] * 3
              + [('H', x), ('S', x), ('H', x)])
    return k + middle + _inverse(k)


def _translate(p, address, basis):
    return sum(((basis >> wire) & 1) << (wire ^ address)
               for wire in range(1 << p))


class DirtyCopyGradingTests(unittest.TestCase):
    def test_noncommuting_pauli_pair_has_an_exact_clifford_copy(self):
        copy = [('CX', 0, 1), ('H', 0), ('CX', 0, 2), ('H', 0)]
        reverse = _inverse(copy)
        sources = (_native_x(0), [('S', 0), ('S', 0)])
        for index, source in enumerate(sources):
            for basis in range(8):
                self.assertEqual(_act(reverse + source + copy, basis),
                                 _act(source + _native_x(index + 1), basis))
                self.assertEqual(_act(copy + reverse, basis), {basis: CO})
        x, z = (native_columns(1, source) for source in sources)
        self.assertEqual(multiply(x, z), scale(multiply(z, x), -1))
        for word in (x, z, multiply(x, z)):
            self.assertEqual(sum((word[i][i] for i in range(2)), CZ), CZ)

    def test_weighted_flip_seed_phase_covariance_and_nonscalar_commutator(self):
        root = Q2(0, Fraction(1, 2))
        omega = (CO + CI) * root
        matrices = []
        for x in range(4):
            word = _weighted_flip(x)
            reverse = _inverse(word)
            self.assertEqual(len(word), 11)
            self.assertEqual(len(reverse), 11)
            self.assertTrue(all(gate[0] in ('H', 'S', 'CX') for gate in word))
            actual = native_columns(4, word)
            expected = zeros(16)
            for column in range(16):
                parity = sum((column >> wire) & 1
                             for wire in (x, x ^ 1, x ^ 2)) % 2
                phase = omega if parity else omega.conj()
                expected[column ^ (1 << x)][column] = phase
                self.assertEqual(_act(word, column),
                                 {column ^ (1 << x): phase})
            self.assertEqual(actual, expected)
            self.assertEqual(native_columns(4, reverse), adjoint(actual))
            self.assertEqual(adjoint(actual), actual)
            self.assertEqual(multiply(actual, actual), eye(16))
            matrices.append(actual)
        for address in range(4):
            for x in range(4):
                for column in range(16):
                    source = _translate(2, address, column)
                    transported = {_translate(2, address, row): value
                                   for row, value in _act(_weighted_flip(x), source).items()}
                    self.assertEqual(transported,
                                     _act(_weighted_flip(x ^ address), column))
        commutator = multiply(multiply(matrices[0], matrices[1]),
                              multiply(matrices[0], matrices[1]))
        z_pair = native_columns(4, [('S', 2), ('S', 2), ('S', 3), ('S', 3)])
        self.assertEqual(commutator, z_pair)
        self.assertNotEqual(commutator, eye(16))
        self.assertNotEqual(commutator, scale(eye(16), -1))

    def test_weighted_flip_copy_covers_every_dirty_carrier_and_output_column(self):
        copy = [('CX', wire, 4 + wire) for wire in range(4)]
        reverse = _inverse(copy)
        for x in range(4):
            source = _weighted_flip(x)
            for column in range(256):
                self.assertEqual(_act(reverse + source + copy, column),
                                 _act(source + _native_x(4 + x), column))
        for column in range(256):
            self.assertEqual(_act(copy + reverse, column), {column: CO})
        # The all-word proof uses the permutation signature. These
        # ordered representatives and an interleaved word check it.
        words = [[x for x in range(4) if mask >> x & 1] for mask in range(1, 16)]
        words.append([0, 1, 0, 2, 1])
        for indices in words:
            parity = 0
            native = []
            for x in indices:
                parity ^= 1 << x
                native += _weighted_flip(x)
            self.assertNotEqual(parity, 0)
            trace = CZ
            for column in range(16):
                state = _act(native, column)
                self.assertEqual(set(state), {column ^ parity})
                trace += state.get(column, CZ)
            self.assertEqual(trace, CZ)

    def test_small_word_trace_does_not_remove_the_spectral_error_floor(self):
        # Individual generators are balanced; their length-two word is
        # an unbalanced reflection of trace 4 in dimension 4096.
        half = 1 << 11
        j = [1] * (half // 2 + 1) + [-1] * (half // 2 - 1)
        r0 = [1] * half + [-1] * half
        r1 = j + [-value for value in j]
        word = [left * right for left, right in zip(r0, r1)]
        self.assertEqual(sum(r0), 0)
        self.assertEqual(sum(r1), 0)
        self.assertEqual(sum(word), 4)
        self.assertEqual(Fraction(sum(word), len(word)), Fraction(1, 1024))
        # Conjugation preserves the left multiplicity. Nonzero output
        # X parity balances the right one: distance 2, copy error >=1.
        left_positive = 4 * word.count(1)
        right_positive = 2 * len(word)
        self.assertNotEqual(left_positive, right_positive)
        self.assertEqual(left_positive - right_positive, 8)

    def test_unbalanced_shear_has_a_common_commuting_sector_and_sharp_pair_gap(self):
        # A one-address field shear with one arbitrary identity helper.
        query = native_columns(3, [('CX', 0, 1)])
        r = native_columns(3, [('H', 0)])
        s = multiply(query, multiply(r, query))
        self.assertEqual(sum((query[i][i] for i in range(8)), CZ), CO * 4)
        self.assertNotEqual(multiply(r, s), multiply(s, r))
        anticommutator = add(multiply(r, s), multiply(s, r))
        root = Q2(0, Fraction(1, 2))
        for u in range(2):
            for helper in range(2):
                column = [[CZ] for _ in range(8)]
                column[u | (helper << 2)][0] = CO * root
                column[u | 2 | (helper << 2)][0] = CO * root
                self.assertEqual(multiply(anticommutator, column), scale(column, 2))
        z = native_columns(1, [('S', 0), ('S', 0)])
        a = native_columns(1, [('H', 0)])
        b = multiply(z, multiply(a, z))
        self.assertEqual(add(multiply(a, b), multiply(b, a)), zeros(2))
        for target in (a, b):
            difference = add(target, scale(z, -1))
            self.assertEqual(multiply(difference, difference),
                             scale(eye(2), Q2(2, -1)))


if __name__ == '__main__':
    unittest.main()
