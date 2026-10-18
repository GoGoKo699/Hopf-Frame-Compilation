"""Exact falsification fixtures for projected three-band proof steps.

These native words test dependent singleton cancellation, an ordered-Pauli
scalar cocycle, and fourth-degree phase-bit carries. They do not infer a
general degree bound from finite examples or discard the unused clean code.
"""
import unittest

from tests.exact_matrices import (
    CZ, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_projected_band_matching import _cq_root
from tests.test_two_layer_diagonal_rigidity import _primitive_inverse
from tests.test_two_layer_galois_return import _conjugate, _tau


def _anf(values):
    coefficients = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for mask in range(len(values)):
            if mask >> bit & 1:
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return {mask for mask, coefficient in enumerate(coefficients) if coefficient}


class PauliBoundaryThreeBandTests(unittest.TestCase):
    def test_even_relation_singleton_retains_the_full_unused_sector_sign(self):
        for relation in (4, 6):
            with self.subTest(relation=relation):
                width, clean = relation, relation - 1
                size = 1 << width
                compute = [('CX', bit, clean) for bit in range(clean)]
                e_word = compute + [('T', bit) for bit in range(width)] + list(reversed(compute))
                e = native_columns(width, e_word)
                self.assertEqual(native_columns(width, _primitive_inverse(e_word)), adjoint(e))
                # Physical axes are independent. Their projected labels have
                # exactly one relation: Z_0,...,Z_(n-1), product(data Z).
                physical_axes = [1 << bit for bit in range(clean)] + [(1 << width) - 1]
                projected_axes = [axis & ((1 << clean) - 1) for axis in physical_axes]
                self.assertEqual(_rank(physical_axes), relation)
                self.assertEqual(_rank(projected_axes), relation - 1)
                b_label = (1 << clean) - 1
                b = _tau(b_label, width)
                for axis in physical_axes:
                    axis_matrix = _tau(axis << width, width)
                    self.assertEqual(multiply(axis_matrix, b), scale(multiply(b, axis_matrix), -1))
                b_word = [gate for bit in range(clean)
                          for gate in [('H', bit), ('S', bit), ('S', bit), ('H', bit)]]
                self.assertEqual(native_columns(width, b_word), b)
                # J=E B E† is an actual one-band Pauli image. The input
                # logical Z family can choose B as its first generator;
                # every other generator is an even-weight data Z parity.
                j_word = _primitive_inverse(e_word) + b_word + e_word
                j = native_columns(width, j_word)
                self.assertEqual(j, _conjugate(e, b))
                self.assertEqual(native_columns(width, _primitive_inverse(j_word)), adjoint(j))
                actual_word = j_word + b_word + j_word
                actual = native_columns(width, actual_word)
                expected = scale(multiply(_tau(1 << (width + clean), width), b),
                                 (-1) ** (relation // 2))
                self.assertEqual(actual, expected)
                self.assertEqual(native_columns(width, _primitive_inverse(actual_word)), adjoint(actual))
                for column in range(size):
                    expected_sign = (-1) ** ((relation // 2) + ((column >> clean) & 1))
                    for row in range(size):
                        self.assertEqual(actual[row][column], b[row][column] * expected_sign)
                # The projected singleton is +B for r=4 and -B for r=6,
                # despite anticommuting with every projected native axis.
                # The full-space image has the indispensable extra Z_c.
                self.assertNotEqual(actual, scale(b, (-1) ** (relation // 2)))

    def test_ordered_pauli_blocks_retain_i_and_the_quadratic_scalar_cocycle(self):
        # Central controls a,b are arbitrary. Target t is the quantum core.
        # In operator order: S_a CZ_(b,t) CX_(a,t) H_t.
        width, a, b, target = 3, 0, 1, 2
        word = [('H', target), ('CX', a, target),
                ('H', target), ('CX', b, target), ('H', target), ('S', a)]
        actual = native_columns(width, word)
        self.assertEqual(native_columns(width, _primitive_inverse(word)), adjoint(actual))
        h = native_columns(1, [('H', 0)])
        x, z = _tau(1, 1), _tau(2, 1)
        for avalue in (0, 1):
            for bvalue in (0, 1):
                controls = avalue | (bvalue << 1)
                ordered = multiply(x if avalue else eye(2), z if bvalue else eye(2))
                ordered = multiply(ordered, h)
                # Z^b X^a=(-1)^(ab) X^a Z^b. Its scalar is literal.
                expected = scale(ordered, _cq_root(2 * avalue + 4 * avalue * bvalue))
                for column in range(2):
                    full_column = controls | (column << target)
                    for full_row in range(1 << width):
                        value = (expected[full_row >> target][column]
                                 if full_row & 3 == controls else CZ)
                        self.assertEqual(actual[full_row][full_column], value)
                if avalue and bvalue:
                    self.assertNotEqual(expected, scale(ordered, _cq_root(2 * avalue)))

    def test_native_linear_phase_has_a_fourth_degree_top_binary_bit(self):
        width = 4
        word = [('T', bit) for bit in range(width)]
        actual = native_columns(width, word)
        self.assertEqual(native_columns(width, _primitive_inverse(word)), adjoint(actual))
        roots = [_cq_root(exponent) for exponent in range(8)]
        observed_exponents = []
        for column in range(1 << width):
            for row in range(1 << width):
                if row != column:
                    self.assertEqual(actual[row][column], CZ)
            observed_exponents.append(roots.index(actual[column][column]))
        # Derive the bits from actual native amplitudes, not a predicted
        # ANF. The phase polynomial itself is linear modulo eight.
        self.assertEqual(observed_exponents, [x.bit_count() for x in range(1 << width)])
        for bit, expected_degree in enumerate((1, 2, 4)):
            coefficients = _anf([(exponent >> bit) & 1 for exponent in observed_exponents])
            self.assertEqual(coefficients,
                             {mask for mask in range(1 << width)
                              if mask.bit_count() == expected_degree})
        self.assertEqual(_anf([(exponent >> 2) & 1 for exponent in observed_exponents]), {15})


if __name__ == '__main__':
    unittest.main()
