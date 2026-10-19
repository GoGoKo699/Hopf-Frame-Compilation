"""Exact native fixture for a nonconstant determinant defect and its kernel.

There is one initialized wire.  The other two wires are arbitrary data;
all full physical columns and the actual reversed words are retained.
The occupied target is identity, so this tests the correction mechanism,
not a difficult target or the general degree theorem.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _clifford_columns, _conjugate, _tau


class DeterminantDefectDegreeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.width = 3
        # K=fanout(c -> data) H_c H_1, so K X_1 K†=Z_c Z_1.
        cls.k_word = [('H', 1), ('H', 0), ('CX', 0, 1), ('CX', 0, 2)]
        cls.e_word = [('T', wire) for wire in range(3)]
        # D Z_1 D†=X_1 and D Z_c D†=Z_c X_1.
        d = [('CX', 1, 0), ('H', 1)]
        cls.e0_word = _inverse(d) + [('TDG', 0), ('T', 1)] + d
        cls.correction_word = [('CX', 0, 2)]
        cls.uc_word = cls.k_word + cls.e_word
        cls.uab_word = cls.correction_word + cls.e0_word + cls.uc_word
        cls.displayed_word = cls.uab_word + _inverse(cls.uc_word)
        cls.baseline_word = cls.correction_word + cls.e0_word
        cls.uc = native_columns(3, cls.uc_word)
        cls.uab = native_columns(3, cls.uab_word)
        cls.actual = native_columns(3, cls.displayed_word)
        cls.e0 = native_columns(3, cls.e0_word)
        cls.correction = native_columns(3, cls.correction_word)
        cls.clean = _tau(1 << 3, 3)
        cls.s = _conjugate(cls.uc, cls.clean)
        cls.r = _tau((1 | 2) << 3, 3)
        cls.plus = scale(add(eye(8), cls.s), Fraction(1, 2))
        cls.minus = scale(add(eye(8), scale(cls.s, -1)), Fraction(1, 2))

    def test_actual_native_transition_inverse_and_paid_correction(self):
        self.assertEqual(self.actual, multiply(adjoint(self.uc), self.uab))
        self.assertEqual(self.actual, native_columns(3, self.baseline_word))
        for word in (self.displayed_word, self.baseline_word):
            self.assertEqual(native_columns(3, _inverse(word)), adjoint(self.actual))
        self.assertNotEqual(self.actual, adjoint(self.actual))
        self.assertEqual(_conjugate(self.uab, self.clean), self.s)
        columns = _clifford_columns(self.s, 3)
        self.assertEqual(_rank([label ^ (1 << bit)
                                for bit, label in enumerate(columns)]), 3)
        for column in range(0, 8, 2):
            self.assertEqual([self.actual[row][column] for row in range(8)],
                             [CO if row == column else CZ for row in range(8)])
        # The correction is one actual CNOT and fixes every occupied column.
        corrected = multiply(self.correction, self.actual)
        self.assertEqual(corrected, self.e0)
        self.assertEqual(adjoint(self.correction), self.correction)
        self.assertEqual(native_columns(3, _inverse(self.e0_word)), adjoint(corrected))
        self.assertEqual(_conjugate(corrected, self.clean), self.clean)
        z2 = _tau(1 << 5, 3)
        self.assertEqual(_conjugate(corrected, z2), z2)
        self.assertNotEqual(_conjugate(corrected, _tau(1 << 4, 3)), _tau(1 << 4, 3))
        for word in (self.displayed_word, _inverse(self.displayed_word)):
            counts = Counter(gate[0] for gate in word)
            self.assertEqual(counts['T'] + counts['TDG'], 8)
            self.assertEqual(counts['H'] + counts['CX'], 13)

    def test_complete_constant_defect_subgroup_and_zero_probe_characters(self):
        # Nonempty constant-defect subgroup H=<Z_2>, with chi(Z_2)=1.
        # Its other coset has alpha=iR or -iR, a nonconstant central phase.
        constant_labels = []
        for subset in range(4):
            z = _tau((subset << 1) << 3, 3)
            c = _conjugate(self.uc, z)
            a = _conjugate(self.uab, z)
            w = multiply(c, a)
            alpha = scale(self.r, CI) if subset & 1 else eye(8)
            if subset & 2:
                alpha = scale(alpha, -1)
            self.assertEqual(w, add(self.plus, multiply(alpha, self.minus)))
            self.assertEqual(multiply(w, self.plus), self.plus)
            if w == eye(8) or w == self.s:
                constant_labels.append(subset)
            # N is the even-weight Z subgroup: check every member, and
            # every logical product, rather than only one generator.
            for mask in (0, 3, 5, 6):
                probe = _tau(mask << 3, 3)
                self.assertEqual(_conjugate(a, probe), _conjugate(c, probe))
        self.assertEqual(constant_labels, [0, 2])
        z2 = _tau(1 << 5, 3)
        self.assertEqual(_conjugate(self.actual, z2), multiply(z2, self.clean))
        self.assertEqual(_conjugate(self.correction, z2), multiply(z2, self.clean))


if __name__ == '__main__':
    unittest.main()
