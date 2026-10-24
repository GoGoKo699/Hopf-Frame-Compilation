"""Exact complete-space guards for the Bell-envelope acceptance transfer.

The envelope is built from its literal Clifford words.  No auxiliary
input is initialized or traced out, and the non-diagonal effect below
retains coherent auxiliary matrix entries.  These finite identities
support the analytic reduction; they do not synthesize its missing effect.
"""
from fractions import Fraction
import unittest

import numpy as np

from tests.exact_matrices import CQ2, native_columns
from tests.test_corank_fourier_bell_source import _kernel, _w_word


def _pair_numerator(rows):
    """Return q^2 W_pair in the C,D,Y,Z register ordering, for n=2."""
    n, q = len(rows), 1 << len(rows)
    assert n == 2
    native = native_columns(2*n, _w_word(rows, [0, 1], [2, 3]))
    one = np.empty((q*q, q*q), dtype=np.int64)
    for i, row in enumerate(native):
        for j, value in enumerate(row):
            if value == CQ2(Fraction(1, q)):
                one[i, j] = 1
            else:
                assert value == CQ2(Fraction(-1, q))
                one[i, j] = -1
    labels = np.arange(q**4)
    c, d, y, z = [(labels >> (n*j)) & (q-1) for j in range(4)]
    cy, dz = c | (y << n), d | (z << n)
    pair = one[cy[:, None], cy[None, :]] * one[dz[:, None], dz[None, :]]
    return pair, labels, y ^ z


class BellAcceptanceTransferTests(unittest.TestCase):
    def test_literal_paired_envelope_moves_the_complete_charge_to_auxiliaries(self):
        for x in (0, 1):
            rows, vector, n = _kernel(1, x)
            pair, labels, difference = _pair_numerator(rows)
            dimension = len(labels)
            np.testing.assert_array_equal(pair.T @ pair,
                                          dimension*np.eye(dimension, dtype=np.int64))
            # Q=X_C(v)X_D(v), including identity on arbitrary Y,Z.
            moved_rows = labels ^ vector ^ (vector << n)
            actual = pair.T @ pair[moved_rows]
            signs = np.array([1-2*((int(value) & vector).bit_count() & 1)
                              for value in difference], dtype=np.int64)
            np.testing.assert_array_equal(actual, dimension*np.diag(signs))
            # The right null vector e_0 would give the wrong nonzero-address
            # effect: the transfer uses the complete left kernel v_x.
            if x:
                wrong = np.array([1-2*(int(value) & 1) for value in difference])
                self.assertFalse(np.array_equal(signs, wrong))

    def test_nondiagonal_auxiliary_effect_retains_the_exact_whole_space_defect(self):
        rows, vector, n = _kernel(1, 1)
        pair, labels, difference = _pair_numerator(rows)
        dimension, q = len(labels), 1 << n
        # E=|psi><psi|/2 with psi=(|0,0>+|1,1>)/sqrt(2).
        # E has coherent off-diagonal entries across both auxiliary banks.
        numerator = np.zeros((q*q, q*q), dtype=np.int64)
        support = (0, 1 | (1 << n))
        numerator[np.ix_(support, support)] = 1  # E=numerator/4.
        effect = np.kron(numerator, np.eye(q*q, dtype=np.int64))
        self.assertEqual(effect[(1 | (1 << n)) << (2*n), 0], 1)
        moved_rows = labels ^ vector ^ (vector << n)
        charge = np.eye(dimension, dtype=np.int64)[moved_rows]
        signs = np.array([1-2*((int(value) & vector).bit_count() & 1)
                          for value in difference], dtype=np.int64)
        gram_numerator = pair @ effect @ pair.T
        defect_numerator = gram_numerator - dimension*(np.eye(dimension, dtype=np.int64) + charge)
        # Both sides are exact integer numerators.  The equality includes
        # every coherent auxiliary input; diagonal-block tests would miss E.
        transported = pair.T @ defect_numerator @ pair
        expected = dimension**2*(effect - np.eye(dimension, dtype=np.int64) - np.diag(signs))
        np.testing.assert_array_equal(transported, expected)

    def test_address_independent_effect_has_the_sharp_quarter_floor(self):
        for p in (1, 2):
            vectors = [_kernel(p, x)[1] for x in range(1 << p)]
            n = _kernel(p, 0)[2]
            labels = range(1 << n)
            # A single auxiliary difference y+z suffices to represent
            # all diagonal entries; each has the omitted multiplicity 2^n.
            target_numerators = [np.array([1 + (1-2*((vector & y).bit_count() & 1))
                                          for y in labels], dtype=np.int64)
                                 for vector in vectors]
            for index, left in enumerate(target_numerators):
                for right in target_numerators[index+1:]:
                    self.assertEqual(int(np.max(np.abs(left-right))), 2)
                # I/4 is a positive contraction and simultaneously has
                # distance exactly 1/4 from every target half-projector.
                self.assertEqual(int(np.max(np.abs(left-1))), 1)

    def test_hyperplane_transvections_fix_both_uniform_cosets(self):
        n, vector = 3, 0b101
        labels = np.arange(1 << n)
        parity = np.array([(int(y) & vector).bit_count() & 1 for y in labels])
        even, odd = (parity == 0).astype(np.int64), (parity == 1).astype(np.int64)
        count = 0
        for direction in range(1, 1 << n):
            if (direction & vector).bit_count() & 1:
                continue
            for covector in range(1, 1 << n):
                if (direction & covector).bit_count() & 1:
                    continue
                images = np.array([int(y) ^ (direction if (int(y) & covector).bit_count() & 1 else 0)
                                   for y in labels])
                self.assertEqual(len(set(images)), len(labels))
                np.testing.assert_array_equal(parity[images], parity)
                np.testing.assert_array_equal(even[images], even)
                np.testing.assert_array_equal(odd[images], odd)
                count += 1
        self.assertEqual(count, 9)
        # The two normalized coherent vectors are orthogonal, but their
        # desired accepted-Gram expectations differ by exactly one half.
        self.assertEqual(int(even @ odd), 0)
        self.assertEqual(Fraction(int(even @ even), 2*int(even.sum())), Fraction(1, 2))
        self.assertEqual(int(odd @ even), 0)


if __name__ == '__main__':
    unittest.main()
