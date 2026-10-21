"""Exact guards for the normalized, thin scalar-kernel source.

The finite native case includes both coherent addresses, all bank and
private-helper columns, literal signal phases, and the actual inverse.
It uses an exact two-control zero XOR, not the asymptotic approximate
expander emitter; the latter's error ledger is analytic.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, adjoint, eye, multiply, native_columns,
)
from tests.test_bent_gram_source import ry_t
from tests.test_half_bilinear_gram_source import _private_toffoli
from tests.test_majorana_walsh_query import _inverse


def x_word(wire):
    return [('H', wire), ('S', wire), ('S', wire), ('H', wire)]


def small_source():
    # x,C0,C1,D0,D1,a,phase. The helper is arbitrary, not initialized.
    x, c0, c1, d0, d1, signal, helper = range(7)
    # Physical-to-dual B=C_E S. At p=1, E^{-T} maps(b0,b1) to
    # (b0,b0+b1), and t is the transformed physical C+D bank.
    basis = [('CX', c0, d0), ('CX', c1, d1),
             ('CX', c0, c1), ('CX', d0, d1)]
    query = (_private_toffoli(x, c1, d0, helper)
             + [('CX', c1, d1)])
    zero = (x_word(d0) + x_word(d1)
            + _private_toffoli(d0, d1, signal, helper)
            + x_word(d1) + x_word(d0))
    query_f = zero + query + zero + _inverse(query)
    middle = query_f + x_word(signal)
    core = (_inverse(middle) + ry_t(signal, inverse=True)
            + middle + ry_t(signal, cube=True))
    return basis + core + _inverse(basis)


def small_f(column):
    x, c0, c1, d0, d1 = ((column >> i) & 1 for i in range(5))
    z1 = c0 ^ c1
    t = (c0 ^ d0) | ((c0 ^ c1 ^ d0 ^ d1) << 1)
    shift = (x & z1) | (z1 << 1)
    return (t == 0) ^ ((t ^ shift) == 0)


class ScalarKernelFilterSourceTests(unittest.TestCase):
    def test_literal_source_all_columns_and_actual_inverse(self):
        signal, width = 5, 7
        word = small_source()
        for emitted, cliffords in ((word, 177), (_inverse(word), 179)):
            self.assertEqual(sum(g[0] in {'T', 'TDG'} for g in emitted), 66)
            self.assertEqual(sum(g[0] not in {'T', 'TDG'} for g in emitted), cliffords)
        actual = native_columns(width, word)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(multiply(adjoint(actual), actual), eye(1 << width))
        signal_even = native_columns(1, ry_t(0) + ry_t(0))
        root_half = Q2(0, Fraction(1, 2))
        omega_inverse = CQ2(root_half, -root_half)
        signal_odd = [[CZ, -CI], [CI, CZ]]
        for column in range(1 << width):
            f = small_f(column)
            for row in range(1 << width):
                if (row & ~(1 << signal)) != (column & ~(1 << signal)):
                    expected = CZ
                elif f:
                    expected = signal_even[(row >> signal) & 1][(column >> signal) & 1]
                else:
                    expected = omega_inverse * signal_odd[(row >> signal) & 1][(column >> signal) & 1]
                self.assertEqual(actual[row][column], expected)
        logical = [column for column in range(1 << width)
                   if not column >> signal & 1]
        accepted = [[actual[row][column] for column in logical] for row in logical]
        gram = multiply(adjoint(accepted), accepted)
        half = CQ2(Q2(Fraction(1, 2)))
        accepted_by_address = [0, 0]
        for i, row in enumerate(logical):
            if small_f(row):
                accepted_by_address[row & 1] += 1
                c = (row >> 1) & 3
                d = (row >> 3) & 3
                self.assertEqual(((c ^ d) >> (row & 1)) & 1, 0)
            for j, column in enumerate(logical):
                self.assertEqual(gram[i][j],
                                 half if i == j and small_f(column) else CZ)
        # Four bank columns per address, each with both helper values.
        self.assertEqual(accepted_by_address, [8, 8])
        self.assertFalse(small_f(0))  # Missing even-charge all-zero input.

    def test_corank_one_support_count_and_shared_zero_eigenvalue(self):
        for n in (1, 2, 3, 5):
            size = 1 << n
            # K e0=0; K ej=e0+ej for j>0. Its image is even parity.
            def kernel_image(z):
                tail = z & ~1
                return tail | (tail.bit_count() & 1)
            accepted = 0
            even = 0
            missing_even = 0
            rejected_odd = 0
            for z in range(size):
                shift = kernel_image(z)
                for t in range(size):
                    f = (t == 0) ^ ((t ^ shift) == 0)
                    charge = t.bit_count() & 1
                    if f:
                        self.assertEqual(charge, 0)
                        accepted += 1
                    if charge == 0:
                        even += 1
                        missing_even += not f
                    else:
                        self.assertFalse(f)
                        rejected_odd += 1
            self.assertEqual(accepted, 2 ** (n + 1) - 4)
            self.assertEqual(even, 2 ** (2 * n - 1))
            self.assertGreater(missing_even, 0)
            self.assertGreater(rejected_odd, 0)

    def test_fixed_diagonal_phase_halfspace_rigidity(self):
        # Three carrier bits: two specified bank wires and one arbitrary
        # extra carrier. Exhaust every real diagonal unitary (256 choices).
        n, bank = 3, 2
        size = 1 << n
        satisfying = 0
        for signs in range(1 << size):
            phase = [(-1 if signs >> t & 1 else 1) for t in range(size)]
            cancellation = all(
                phase[t ^ h] == phase[t]
                for x in range(bank)
                for t in range(size) if t >> x & 1
                for h in range(size) if not h >> x & 1)
            outside = [t for t in range(size) if t & ((1 << bank) - 1)]
            constant_outside = len({phase[t] for t in outside}) == 1
            self.assertEqual(cancellation, constant_outside)
            if cancellation:
                satisfying += 1
                common = phase[outside[0]]
                self.assertLessEqual(sum(value != common for value in phase),
                                     1 << (n - bank))
        self.assertEqual(satisfying, 8)


if __name__ == '__main__':
    unittest.main()
