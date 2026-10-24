"""Exact fixtures for a four-band broad source with a wrong Pauli image.

The retained coefficient excludes this emitted commutator, not all shifted
syndrome words, four-layer circuits, or approximate complete sources.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, CQ2, Q2, add, adjoint, multiply, native_columns, scale,
)
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _shifted_x_band(p, shift=1, signs=None):
    """Signed X-axis rotations, with scalar exp(i pi sum(signs)/8)."""
    size = 1 << p
    assert p >= 1 and 0 <= shift < size
    signs = (1,) * size if signs is None else signs
    assert len(signs) == size and set(signs) <= {-1, 1}
    bank = range(p, p + size)
    wall = [('CX', p + row, bit)
            for row in range(size) for bit in range(p)
            if ((row ^ shift) >> bit) & 1]
    h = [('H', wire) for wire in bank]
    t = [('T' if signs[row] == 1 else 'TDG', p + row)
         for row in range(size)]
    return (wall + h + t + h
            + list(reversed(wall)))


def _unshifted_z_band(p, signs=None):
    """Signed Z-axis rotations, with scalar exp(i pi sum(signs)/8)."""
    size = 1 << p
    signs = (1,) * size if signs is None else signs
    assert len(signs) == size and set(signs) <= {-1, 1}
    wall = [('CX', bit, p + row)
            for row in range(size) for bit in range(p)
            if (row >> bit) & 1]
    t = [('T' if signs[row] == 1 else 'TDG', p + row)
         for row in range(size)]
    return (wall + t
            + list(reversed(wall)))


def _commutator(p, a_signs=None, b_signs=None):
    a = _shifted_x_band(p, signs=a_signs)
    b = _unshifted_z_band(p, signs=b_signs)
    # Chronological B, A, B-dagger, A-dagger gives E=A†B†AB.
    return b + a + _inverse(b) + _inverse(a)


def _pauli_echo_coefficient(matrix, bit):
    """Normalized Tr(Z_bit E Z_bit E†), with all native entries exact."""
    return sum((
        (-1 if ((row ^ column) >> bit) & 1 else 1)
        * value * value.conj()
        for row, entries in enumerate(matrix)
        for column, value in enumerate(entries)), CZ) / len(matrix)


class ShiftedSyndromeCommutatorTests(unittest.TestCase):
    def test_complete_small_source_inverse_and_retained_pauli_coefficient(self):
        expected = {
            1: 3 * CO / 4,
            2: CQ2(Q2(Fraction(1, 2), Fraction(1, 8))),
        }
        for p in (1, 2):
            word = _commutator(p)
            actual = native_columns(p + (1 << p), word)
            self.assertEqual(
                native_columns(p + (1 << p), _inverse(word)),
                adjoint(actual))
            coefficient = _pauli_echo_coefficient(actual, p)
            self.assertEqual(coefficient, expected[p])
            self.assertNotEqual(coefficient, CO / (1 << p))

    def test_full_echo_orientation_and_mixed_t_orientations(self):
        # The correlation alone is unchanged by E -> E†. Compare the full
        # p=1 operator to A†(P+B†RB)A/sqrt(2), with R=-i P_0 P explicit.
        p = 1
        width = p + (1 << p)
        size = 1 << width
        a = native_columns(width, _shifted_x_band(p))
        b = native_columns(width, _unshifted_z_band(p))
        e = native_columns(width, _commutator(p))
        probe = [[((-CO if (column >> p) & 1 else CO)
                   if row == column else CZ)
                  for column in range(size)] for row in range(size)]
        r = [[((-CI if not (column >> p) & 1 else CI)
               if row == column ^ (1 << p) ^ 1 else CZ)
              for column in range(size)] for row in range(size)]
        r_prime = multiply(adjoint(b), multiply(r, b))
        rhs = multiply(adjoint(a), multiply(
            scale(add(probe, r_prime), CQ2(Q2(0, Fraction(1, 2)))), a))
        self.assertEqual(multiply(e, multiply(probe, adjoint(e))), rhs)
        self.assertNotEqual(multiply(adjoint(e), multiply(probe, e)), rhs)

        signed = _commutator(2, (-1, 1, -1, 1), (1, -1, -1, 1))
        self.assertEqual(
            _pauli_echo_coefficient(native_columns(6, signed), 2),
            CQ2(Q2(Fraction(1, 2), Fraction(1, 8))))

    def test_native_bands_ledger_and_anticommuting_axis_count(self):
        for p in (1, 2, 3, 6):
            size = 1 << p
            for word in (_commutator(p), _inverse(_commutator(p))):
                bands = _bands(word)
                self.assertEqual(list(map(len, bands)), [size] * 4)
                for band in bands:
                    self.assertEqual(len({gate[1] for gate in band}), size)
                counts = Counter(gate[0] for gate in word)
                self.assertEqual(counts['T'] + counts['TDG'], 4 * size)
                self.assertEqual(counts['CX'] + counts['H'],
                                 4 * p * size + 4 * size)
                self.assertEqual(set(counts), {'CX', 'H', 'T', 'TDG'})
            # v=e_0; the j=0 term and the v.j=1 labels are disjoint.
            anticommuting = [row for row in range(size)
                             if (row & 1) ^ int(row == 0)]
            self.assertEqual(len(anticommuting), size // 2 + 1)


if __name__ == '__main__':
    unittest.main()
