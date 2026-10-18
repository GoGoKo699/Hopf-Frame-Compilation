"""Bounded checks of the one-clean three-layer degree argument.

The width-independent theorem is analytic. These checks target the
stabilizer-overlap gap, its empty-intersection exception, and the precise
codimension used in the degree proof; they do not extrapolate a finite
circuit search into a lower bound.
"""
import unittest

import numpy as np

from compiler_robust_hopf.native_coarse_fixture import commutator_word
from tests.test_operator_source_compiler import _word_matrix
from tests.test_read_only_clean_rigidity import _degree


def _pauli(width, label):
    x, z = label & ((1 << width) - 1), label >> width
    word = [('X', bit) for bit in range(width) if x & (1 << bit)]
    word += [('Z', bit) for bit in range(width) if z & (1 << bit)]
    return _word_matrix(width, word)


def _code_isometry(projector):
    values, vectors = np.linalg.eigh(projector)
    return vectors[:, values > 0.5]


class ReadOnlyThreeLayerCleanRigidityTests(unittest.TestCase):
    def test_pauli_gap_on_nonzero_commuting_code_intersections(self):
        width, z_label = 3, 1 << 5  # Z on the top, clean wire.
        identity = np.eye(1 << width, dtype=complex)
        z = _pauli(width, z_label)
        # Independent, signed, and identical positive stabilizers.
        for other_label, sign in ((1, 1), (1, -1), (z_label ^ 1, 1),
                                  (z_label, 1)):
            other = sign * _pauli(width, other_label)
            projector = (identity + z) @ (identity + other) / 4
            j = _code_isometry(projector)
            self.assertGreater(j.shape[1], 0)
            span = {0, z_label, other_label, z_label ^ other_label}
            for label in range(1 << (2 * width)):
                if label in span:
                    continue
                p = _pauli(width, label)
                self.assertAlmostEqual(abs(np.trace(j.conj().T @ p @ j)), 0, places=12)
                distance = min(np.linalg.norm(j - phase * p @ j, 2)
                               for phase in (1, -1, 1j, -1j))
                self.assertAlmostEqual(distance, np.sqrt(2), places=12)

    def test_opposite_stabilizers_must_be_excluded(self):
        identity = np.eye(8, dtype=complex)
        z = _pauli(3, 1 << 5)
        intersection = (identity + z) @ (identity - z) / 4
        np.testing.assert_array_equal(intersection, np.zeros_like(identity))
        # There is no Pauli-label constraint on an empty intersection.
        outside_span = _pauli(3, 1)
        np.testing.assert_array_equal((identity - outside_span) @ intersection,
                                      np.zeros_like(identity))

    def test_anticommuting_orbit_extends_with_no_extra_error(self):
        # A small native occupied-block perturbation; the assembled block
        # unitary tests the structural lemma, not a three-T-layer compiler.
        identity = np.eye(4, dtype=complex)
        n = _word_matrix(2, [('X', 0)])
        r = _word_matrix(2, [('Z', 0)])
        word = [(name, *wires) for name, wires in commutator_word(3, 0)]
        occupied = n @ _word_matrix(2, word)
        empty = r @ occupied @ r.conj().T
        zero = np.zeros((4, 4), dtype=complex)
        unitary = np.block([[occupied, zero], [zero, empty]])
        orbit_pauli = np.block([[zero, r.conj().T], [r, zero]])
        clean_z = np.block([[identity, zero], [zero, -identity]])
        full_pauli = np.block([[n, zero], [zero, -n]])
        np.testing.assert_allclose(unitary @ clean_z, clean_z @ unitary,
                                   atol=2e-12, rtol=0)
        np.testing.assert_allclose(unitary @ orbit_pauli, orbit_pauli @ unitary,
                                   atol=2e-12, rtol=0)
        occupied_error = np.linalg.norm(occupied - n, 2)
        self.assertGreater(occupied_error, 1e-8)
        self.assertLess(occupied_error, 1 / 64)
        self.assertAlmostEqual(np.linalg.norm(unitary - full_pauli, 2),
                               occupied_error, places=12)

    def test_codimension_two_derivative_boundary(self):
        # Exhaust all four-variable truth tables. G is spanned by the first
        # two coordinate directions; affine derivatives there force degree≤2.
        affine_tables = set()
        for linear in range(16):
            for constant in (0, 1):
                table = [constant ^ ((linear & x).bit_count() & 1)
                         for x in range(16)]
                affine_tables.add(sum(value << x for x, value in enumerate(table)))
        admitted = 0
        for table_code in range(1 << 16):
            shifted_one = ((table_code & 0x5555) << 1) | ((table_code & 0xAAAA) >> 1)
            shifted_two = ((table_code & 0x3333) << 2) | ((table_code & 0xCCCC) >> 2)
            if (table_code ^ shifted_one not in affine_tables
                    or table_code ^ shifted_two not in affine_tables):
                continue
            admitted += 1
            table = [(table_code >> x) & 1 for x in range(16)]
            self.assertLessEqual(_degree(table), 2)
        self.assertEqual(admitted, 1 << 11)  # All degree≤2 polynomials in four bits.
        # Codimension three would be insufficient: a cubic can live wholly
        # in the quotient while every G derivative vanishes.
        cubic = [int((x & 14) == 14) for x in range(16)]
        self.assertEqual(_degree(cubic), 3)
        self.assertTrue(all(cubic[x] == cubic[x ^ 1] for x in range(16)))

    def test_scalar_derivative_conversion_keeps_literal_address_phase(self):
        h = _word_matrix(3, [('H', 2)])
        for code in range(256):
            table = [(code >> x) & 1 for x in range(8)]
            scalar = np.diag([(-1) ** value for value in table])
            expected = np.zeros((8, 8), dtype=complex)
            for column in range(8):
                x = column & 3
                derivative = table[x] ^ table[x | 4]
                expected[column ^ (derivative << 2), column] = (-1) ** table[x]
            np.testing.assert_allclose(h @ scalar @ h, expected, atol=1e-14, rtol=0)


if __name__ == '__main__':
    unittest.main()
