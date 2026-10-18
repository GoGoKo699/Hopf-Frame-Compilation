"""Bounded phase/support checks for exact two-layer one-clean rigidity.

The abstract block fixtures are not claimed to be native synthesis.
In particular, a general conjugated Pauli representation can implement
literal AND while failing the physical one-T-conjugation restriction.
"""
from itertools import combinations, product
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, adjoint, eye, multiply, native_columns, scale, zeros,
)


def _hermitian_pauli(width, x, z):
    word = [('Z', j) for j in range(width) if (z >> j) & 1]
    word += [('X', j) for j in range(width) if (x >> j) & 1]
    return scale(native_columns(width, word),
                 (CO, CI, -CO, -CI)[(x & z).bit_count() % 4])


def _support(matrix):
    width = len(matrix).bit_length() - 1
    result = []
    for x, z in product(range(1 << width), repeat=2):
        pauli = _hermitian_pauli(width, x, z)
        relative = multiply(pauli, matrix)
        coefficient = sum((relative[j][j] for j in range(len(matrix))), CZ)
        if coefficient != CZ:
            result.append(pauli)
    return result


def _commute(left, right):
    return multiply(left, right) == multiply(right, left)


def _blocks(upper_left, upper_right, lower_left, lower_right):
    return [a + b for a, b in zip(upper_left, upper_right)] + [
        a + b for a, b in zip(lower_left, lower_right)]


def _product(*matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = multiply(result, matrix)
    return result


class ReadOnlyTwoLayerExactPhaseTests(unittest.TestCase):
    def test_four_term_one_t_support_cannot_be_pairwise_commuting(self):
        clifford = native_columns(2, [('H', 0), ('CX', 0, 1), ('S', 1)])
        for signs in product(('T', 'TDG'), repeat=2):
            layer = native_columns(2, [(signs[0], 0), (signs[1], 1)])
            for z in range(4):
                original = _hermitian_pauli(2, 3, z)
                dressed = _product(clifford, layer, original,
                                   adjoint(layer), adjoint(clifford))
                terms = _support(dressed)
                self.assertEqual(len(terms), 4)
                self.assertTrue(any(not _commute(left, right)
                                    for left, right in combinations(terms, 2)))

    def test_abstract_literal_and_violates_the_physical_support_condition(self):
        # t=0,d=1,c=2. A_b=Z_d CNOT(c,t) is a Hermitian Clifford
        # with balanced spectrum, hence an abstract conjugated Pauli.
        # Its four commuting support terms preclude one-T conjugation.
        a_b = native_columns(3, [('CX', 2, 0), ('Z', 1)])
        terms = _support(a_b)
        self.assertEqual(len(terms), 4)
        self.assertTrue(all(_commute(left, right)
                            for left, right in combinations(terms, 2)))
        self.assertEqual(adjoint(a_b), a_b)
        self.assertEqual(multiply(a_b, a_b), eye(8))
        self.assertEqual(sum((a_b[j][j] for j in range(8)), CZ), CZ)
        for a, b in product((0, 1), repeat=2):
            word = [('X', 2)] if a else []
            if b:
                word += [('CX', 2, 0), ('Z', 1), ('Z', 1)]
            word += [('X', 2)] if a else []
            actual = native_columns(3, word, columns=range(4))
            for column in range(4):
                row = column ^ (a * b)
                self.assertEqual([actual[j][column] for j in range(8)],
                                 [CO if j == row else CZ for j in range(8)])

    def test_distinct_anticommuting_sector_paulis_retain_the_imaginary_phase(self):
        b0 = _hermitian_pauli(1, 0, 1)  # Z
        b1 = _hermitian_pauli(1, 1, 1)  # Y
        x = _hermitian_pauli(1, 1, 0)
        result = multiply(b1, b0)
        self.assertEqual(result, scale(x, CI))
        self.assertEqual(adjoint(result), scale(result, -CO))
        self.assertNotEqual(result, x)
        self.assertNotEqual(scale(result, -CO), x)

    def test_independent_syndrome_block_identity_cancels_actual_complement(self):
        # Pure block-algebra fixture; A_b need not be a native one-T
        # conjugation. The independent-syndrome argument does not use it.
        identity, zero = eye(2), zeros(2)
        x = _hermitian_pauli(1, 1, 0)
        z = _hermitian_pauli(1, 0, 1)
        complement = native_columns(1, [('H', 0), ('T', 0)])
        w = _blocks(identity, zero, zero, complement)
        p_a = r_b = _blocks(zero, identity, identity, zero)
        p_b = _blocks(z, zero, zero, z)
        r_a = _blocks(x, zero, zero, x)
        f = scale(x, -CI)
        g_adjoint = scale(multiply(z, adjoint(complement)), CI)
        a_a = _blocks(zero, adjoint(f), f, zero)
        a_b = _blocks(zero, g_adjoint, adjoint(g_adjoint), zero)
        values = [
            w,
            scale(_product(p_a, a_a, w, r_a), CI),
            scale(_product(p_b, a_b, w, r_b), -CI),
            scale(_product(p_a, p_b, a_a, a_b, w, r_a, r_b), -CO),
        ]
        for value in values:
            self.assertEqual([row[:2] for row in value],
                             identity + zero)


if __name__ == '__main__':
    unittest.main()
