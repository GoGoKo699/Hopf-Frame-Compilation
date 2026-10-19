"""An easy clean target with unbounded distortion at a one-core interface.

The physical native words retain every column and their actual inverses.
The profile is recovered from exact Pauli conjugations, independently of
the symbolic active-axis formula used in the arbitrary-width proof.
"""
from collections import Counter
import unittest

from tests.exact_matrices import CO, CZ, adjoint, eye, multiply, native_columns
from tests.test_clean_code_clifford_rigidity import _pauli_label, _rank
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _conjugate, _symplectic, _tau


def _physical_words(m):
    width = m + 2
    change_b = [gate for bit in range(2, width)
                for gate in [('CX', 0, bit), ('CX', 1, bit)]]
    change_c = [('H', 0), ('CX', 0, 1)] + change_b
    a = [('H', 1), ('T', 0), ('T', 1), ('H', 1)]
    b = _inverse(change_b) + [('TDG', 0)]
    b += [('T', bit) for bit in range(2, width)] + change_b
    c = _inverse(change_c) + [('T', 0)]
    c += [('T', bit) for bit in range(2, width)] + change_c
    word = [('H', 0)] + a + b + _inverse(c) + [('H', 0)]
    zc = 1 << width
    active = [zc | (1 << (width + 1)) | (1 << (width + bit))
              for bit in range(2, width)]
    axes = [[zc, 1 << 1], [zc] + active, active + [3]]
    return word, axes


def _observed_displacement_rank(clifford, width):
    displacements = []
    for bit in range(2 * width):
        label = _pauli_label(_conjugate(clifford, _tau(1 << bit, width)), width)
        if label is None:
            raise AssertionError('expected an actual Clifford Pauli image')
        displacements.append(label ^ (1 << bit))
    return _rank(displacements)


class RankOneCleanPhaseInterfaceTests(unittest.TestCase):
    def test_clean_native_identity_and_unbounded_one_core_profile_jump(self):
        for m in (1, 2, 3):
            with self.subTest(m=m):
                width = m + 2
                word, axes = _physical_words(m)
                actual = native_columns(width, word)
                self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
                # The independent full-space expression retains the unused
                # clean block omega^-1 S_X; it is not full-space identity.
                reduced = [('H', 1), ('T', 1), ('H', 1),
                           ('H', 1), ('CX', 0, 1), ('TDG', 1),
                           ('CX', 0, 1), ('H', 1)]
                self.assertEqual(actual, native_columns(width, reduced))
                self.assertNotEqual(actual, eye(1 << width))
                for column in range(0, 1 << width, 2):
                    for row in range(1 << width):
                        self.assertEqual(actual[row][column], CO if row == column else CZ)
                for axis_list in axes:
                    self.assertEqual(_rank(axis_list), len(axis_list))
                    self.assertTrue(all(not _symplectic(left, right, width)
                                        for left in axis_list for right in axis_list))
                self.assertEqual([len(band) for band in _native_t_bands(word)],
                                 [2, m + 1, m + 1])
                for schedule in (word, _inverse(word)):
                    counts = Counter(gate[0] for gate in schedule)
                    self.assertEqual(counts['T'] + counts['TDG'], 2 * m + 4)
                    self.assertEqual(counts['CX'], 8 * m + 2)
                    self.assertEqual(counts['H'], 6)

                logical_width = m + 1
                fanout = [('CX', 0, bit) for bit in range(1, logical_width)]
                b_word = fanout + [('T', bit) for bit in range(1, logical_width)]
                b_word += _inverse(fanout)
                b = native_columns(logical_width, b_word)
                f = native_columns(logical_width, _inverse(b_word))
                a_word = [('H', 0), ('T', 0), ('H', 0)]
                a = native_columns(logical_width, a_word)
                self.assertEqual(multiply(adjoint(a), multiply(f, multiply(b, a))),
                                 eye(1 << logical_width))
                z_label = 1 << logical_width
                y_label = z_label | 1
                j = _conjugate(a, _tau(z_label, logical_width))
                self.assertEqual(_pauli_label(_conjugate(j, _tau(z_label, logical_width)),
                                             logical_width), y_label)
                ranks = [_observed_displacement_rank(
                    _conjugate(b, _tau(label, logical_width)), logical_width)
                    for label in (z_label, y_label)]
                self.assertEqual(ranks, [0, m])
                middle_axes = [z_label | (1 << (logical_width + bit))
                               for bit in range(1, logical_width)]
                union = middle_axes + [1]
                pairing_rows = [sum(_symplectic(left, right, logical_width) << bit
                                    for bit, right in enumerate(union)) for left in union]
                self.assertEqual(_rank(pairing_rows), 2)


if __name__ == '__main__':
    unittest.main()
