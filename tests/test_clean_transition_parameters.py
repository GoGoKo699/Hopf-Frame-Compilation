"""Identity targets with unbounded native determinant-defect codimension.

One initialized wire is used. Every other wire is arbitrary, and the
checks retain actual inverse words and every physical matrix column.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _conjugate, _tau


def _word(n):
    data_h = [('H', j) for j in range(1, n + 1)]
    fanout = [('CX', 0, j) for j in range(1, n + 1)]
    return (data_h + fanout + [('TDG', j) for j in range(1, n + 1)]
            + fanout + [('H', 0)] + [('T', j) for j in range(n + 1)]
            + data_h + [('TDG', 0), ('H', 0)])


class CleanTransitionParameterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixtures = []
        for n in range(1, 4):
            q = n + 1
            word = _word(n)
            v = native_columns(q, word)
            uc = native_columns(q, [('H', 0), ('T', 0), ('H', 0)])
            s = _conjugate(uc, _tau(1 << q, q))
            plus = scale(add(eye(1 << q), s), Fraction(1, 2))
            minus = scale(add(eye(1 << q), scale(s, -1)), Fraction(1, 2))
            cls.fixtures.append((n, q, word, v, uc, s, plus, minus))

    def test_actual_native_words_literal_return_and_inverse_ledger(self):
        for n, q, word, v, uc, s, _, _ in self.fixtures:
            with self.subTest(n=n):
                # Three distinct-wire physical T bands, with both orientations.
                bands = []
                current = []
                for gate in word + [('H', 0)]:
                    if gate[0] in ('T', 'TDG'):
                        current.append(gate[1])
                    elif current:
                        bands.append(current)
                        current = []
                self.assertEqual([len(b) for b in bands], [n, n + 1, 1])
                self.assertTrue(all(len(set(b)) == len(b) for b in bands))
                counts = Counter(gate[0] for gate in word)
                self.assertEqual(counts['T'] + counts['TDG'], 2 * n + 2)
                self.assertEqual(counts['H'], 2 * n + 2)
                self.assertEqual(counts['CX'], 2 * n)
                inv = _inverse(word)
                self.assertEqual(native_columns(q, inv), adjoint(v))
                self.assertEqual(Counter(gate[0] for gate in inv)['H'], 2 * n + 2)
                self.assertEqual(Counter(gate[0] for gate in inv)['CX'], 2 * n)
                clean = _tau(1 << q, q)
                self.assertEqual(_conjugate(v, clean), clean)
                for col in range(0, 1 << q, 2):
                    self.assertEqual([v[row][col] for row in range(1 << q)],
                                     [CO if row == col else CZ
                                      for row in range(1 << q)])
                # The displayed first clean image is Pauli, the last rank one.
                first_axes = [(1 << j) | (1 << q) for j in range(1, q)]
                self.assertEqual(_rank(first_axes), n)
                columns = []
                for bit in range(2 * q):
                    label = 1 << bit
                    image_label = label | 1 if bit == q else label
                    sign = -1 if bit in (0, q) else 1
                    self.assertEqual(_conjugate(s, _tau(label, q)),
                                     scale(_tau(image_label, q), sign))
                    columns.append(image_label)
                self.assertEqual(_rank([label ^ (1 << bit)
                                        for bit, label in enumerate(columns)]), 1)
                self.assertEqual(_conjugate(multiply(uc, v), clean), s)

    def test_complete_defects_zero_probe_characters_and_trivial_kernel(self):
        for n, q, _, v, uc, s, plus, minus in self.fixtures:
            with self.subTest(n=n):
                uab = multiply(uc, v)
                constant_labels = []
                for z in range(1 << n):
                    logical_z = _tau((z << 1) << q, q)
                    c = _conjugate(uc, logical_z)
                    a = _conjugate(uab, logical_z)
                    w = multiply(c, a)
                    phase = CO
                    for _ in range(z.bit_count()):
                        phase *= CI
                    alpha = scale(_tau(z << 1, q), phase)
                    self.assertEqual(w, add(plus, multiply(alpha, minus)))
                    self.assertEqual(multiply(w, plus), plus)
                    if w == eye(1 << q) or w == s:
                        constant_labels.append(z)
                    # These generate the entire probe space N=<X_1,...,X_n>.
                    # Every logical product z is checked, not only generators.
                    for j in range(1, q):
                        probe = _tau(1 << j, q)
                        self.assertEqual(_conjugate(a, probe),
                                         _conjugate(c, probe))
                self.assertEqual(constant_labels, [0])


if __name__ == '__main__':
    unittest.main()
