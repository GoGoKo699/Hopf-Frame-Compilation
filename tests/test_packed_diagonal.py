"""Exact small native audits of the packed, one-tail diagonal primitive.

All amplitudes lie in Q(sqrt(2), i). Actual gate words retain literal phases,
all dirty input columns, and rejected flags. Finite checks are not resource
proofs or a complete-frame compiler.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import unittest

from tests.exact_matrices import (
    CQ2, Q2, CZ, CO, CI, add, adjoint, eye, multiply, native_columns, scale,
    zeros,
)
from tests.test_one_clean_compiler import _literal_mask, _pauli_t_word
from tests.test_operator_source_compiler import _adjoint, _sign_word


ROOT_HALF = CQ2(Q2(0, F(1, 2)))


def _source_word(q):
    word = []
    for j in range(2 * q - 1):
        factors = {j // 2: 'Z'} if j % 2 == 0 else {j // 2: 'X', j // 2 + 1: 'X'}
        word += _pauli_t_word(factors)
    return word


@lru_cache(None)
def _source(q):
    word = _source_word(q)
    matrix = native_columns(q, _adjoint(word) + [('X', 0)] + word)
    gammas = []
    for j in range(q):
        prefix = [('Z', k) for k in range(j)]
        gammas.append(native_columns(q, prefix + [('X', j)]))
        gammas.append(native_columns(q, prefix + [('SDG', j), ('X', j), ('S', j)]))
    coefficients = []
    value = CO
    for j in range(2 * q):
        if j < 2 * q - 1:
            value *= ROOT_HALF
        coefficients.append(value)
    return matrix, gammas, coefficients


def _controlled_mask(q, cosine_signs, sine_signs):
    selector = q + 1
    word = []
    for value, signs in enumerate((cosine_signs, sine_signs)):
        negate = [('X', selector)] if value == 0 else []
        word += negate
        for gate, target in _literal_mask(signs):
            if gate == 'X':
                word.append(('CX', selector, target))
            else:
                word += [('H', target), ('CX', selector, target), ('H', target)]
        word += negate
    return word


def _diagonal_word(q, cosine_signs, sine_signs):
    source = _source_word(q)
    inverse = _adjoint(source)
    flag, selector = q, q + 1
    center = inverse + [('X', 0)] + source
    control1 = inverse + [('CX', flag, 0)] + source
    control0 = [('X', flag)] + control1 + [('X', flag)]
    mask = _controlled_mask(q, cosine_signs, sine_signs)
    return ([('H', selector), ('H', flag)] + control1 + _adjoint(mask)
            + center + mask + control0
            + [('H', flag), ('S', selector), ('H', selector)])


def _reflection_word(q):
    a, b = q, q + 1
    return [('X', a), ('X', b), ('H', b), ('CX', a, b),
            ('H', b), ('X', b), ('X', a)]


class PackedDiagonalTests(unittest.TestCase):
    def test_native_one_tail_source_every_mask_and_scalar_anticommutator(self):
        for q in range(1, 4):
            source, gammas, coefficients = _source(q)
            size = 1 << q
            expansion = zeros(size)
            for coefficient, gamma in zip(coefficients, gammas):
                expansion = add(expansion, scale(gamma, coefficient))
            self.assertEqual(source, expansion)
            self.assertEqual(multiply(source, source), eye(size))
            self.assertEqual(sum((a * a for a in coefficients), CZ), CO)
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in _source_word(q)), 2 * q - 1)
            for signs in product(range(2), repeat=2 * q):
                with self.subTest(q=q, signs=signs):
                    mask = native_columns(q, _literal_mask(signs))
                    for bit, gamma in zip(signs, gammas):
                        self.assertEqual(multiply(multiply(mask, gamma), adjoint(mask)),
                                         scale(gamma, (-1) ** bit))
                    other = multiply(multiply(mask, source), adjoint(mask))
                    scalar = sum((a * a * (-1) ** bit
                                  for a, bit in zip(coefficients, signs)), CZ)
                    anticommutator = add(multiply(source, other), multiply(other, source))
                    self.assertEqual(anticommutator, scale(eye(size), 2 * scalar))

    def test_literal_complex_block_oaa_and_full_initialized_isometry(self):
        target = CQ2(F(3, 5), F(4, 5))
        for q in range(1, 4):
            # Unequal independently rounded real/imaginary values in the
            # signed grid; the native address masks share the same selector.
            bits = 2 * q - 1
            denominator = 1 << bits
            integers = [round((1 - value) * denominator / 2)
                        for value in (F(3, 5), F(4, 5))]
            signs = [_sign_word(integer, bits) for integer in integers]
            cosine, sine = [F(1) - F(2 * integer, denominator) for integer in integers]
            z = CQ2(cosine, sine)
            size = 1 << q
            word = _diagonal_word(q, *signs)
            output = native_columns(q + 2, word, range(size))
            self.assertEqual(output[:size], scale(eye(size), z / 2))
            self.assertEqual(multiply(adjoint(output), output), eye(size))
            reflection = _reflection_word(q)
            # XZXZ implements the displayed leading minus sign exactly.
            minus = [('X', 0), ('Z', 0), ('X', 0), ('Z', 0)]
            amplified_word = word + reflection + _adjoint(word) + reflection + word + minus
            amplified = native_columns(q + 2, amplified_word, range(size))
            accepted = z * (3 - z.abs2()) / 2
            self.assertEqual(amplified[:size], scale(eye(size), accepted))
            self.assertEqual(multiply(adjoint(amplified), amplified), eye(size))
            error = [row[:] for row in amplified]
            for j in range(size):
                error[j][j] -= target
            h = target.conj() * (z - target)
            a, b = h.real, h.imag
            squared_error = 3 * a * a + b * b + a * (a * a + b * b)
            self.assertEqual(multiply(adjoint(error), error), scale(eye(size), CQ2(squared_error)))
            self.assertLessEqual(squared_error, (3 + abs(a) + abs(b)) * (a * a + b * b))
            grid_half_step = F(1, denominator)
            self.assertLessEqual(squared_error, 50 * grid_half_step * grid_half_step)
            # A wrong inverse changes the complete initialized output.
            wrong = native_columns(q + 2, word + reflection + word + reflection + word + minus,
                                   range(size))
            self.assertNotEqual(wrong, amplified)


if __name__ == '__main__':
    unittest.main()
