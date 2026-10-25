"""Exact complete-bank checks for Clifford-fiber projector-Gram rigidity.

These small literal Clifford words include arbitrary dirty bank columns.
They support the stabilizer and kernel identities, not an asymptotic emitter.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_two_layer_galois_return import _tau


def _rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def _xor_selected(columns, selection):
    result = 0
    for j, column in enumerate(columns):
        if selection >> j & 1:
            result ^= column
    return result


def _inverse(word):
    result = []
    for name, *wires in reversed(word):
        result.extend([(name, *wires)] * (3 if name == 'S' else 1))
    return result


def _compress(matrix, bank_width):
    size = 1 << bank_width
    return [row[:size] for row in matrix[:size]]


def _pauli_labels(operator, width):
    for label in range(1 << (2*width)):
        pauli = _tau(label, width)
        if operator == pauli or operator == scale(pauli, -1):
            return label & ((1 << width)-1), label >> width
    raise AssertionError('literal Clifford image is not a Pauli')


def _gram(width, bank_width, word):
    native = native_columns(width, word)
    accepted = _compress(native, bank_width)
    return multiply(adjoint(accepted), accepted)


class CliffordProjectorGramTests(unittest.TestCase):
    def test_complete_stabilizer_compression_with_interleaved_dirty_banks(self):
        # Every bank column is included; no dirty bank is initialized or traced.
        cases = [
            (2, [('H', 0), ('CX', 0, 2), ('S', 1), ('CX', 2, 1),
                 ('H', 2), ('CX', 1, 0)]),
            (2, [('X', 2)]),  # negative-identity stabilizer gives zero Gram
            (1, [('H', 0), ('CX', 0, 1), ('H', 2), ('S', 0),
                 ('CX', 1, 0), ('H', 1), ('S', 2), ('CX', 0, 2)]),
            (1, [('CX', 0, 1), ('H', 2)]),
            (1, [('H', 1), ('H', 2)]),
        ]
        for bank_width, word in cases:
            with self.subTest(bank_width=bank_width, word=word):
                width, k = 3, 3-bank_width
                native = native_columns(width, word)
                self.assertEqual(native_columns(width, _inverse(word)), adjoint(native))
                self.assertEqual(multiply(adjoint(native), native), eye(1 << width))
                images, clean_x = [], []
                for signal in range(bank_width, width):
                    image = multiply(multiply(adjoint(native),
                                              _tau(1 << (width+signal), width)), native)
                    images.append(image)
                    x, _ = _pauli_labels(image, width)
                    clean_x.append(x >> bank_width)
                summed = scale(eye(1 << bank_width), 0)
                for subset in range(1 << k):
                    if _xor_selected(clean_x, subset):
                        continue
                    element = eye(1 << width)
                    for j, operator in enumerate(images):
                        if subset >> j & 1:
                            element = multiply(element, operator)
                    summed = add(summed, _compress(element, bank_width))
                gram = _gram(width, bank_width, word)
                self.assertEqual(gram, scale(summed, Fraction(1, 1 << k)))
                # Either zero, or 2^-rho times a projector, on the whole bank.
                self.assertEqual(multiply(gram, gram),
                                 scale(gram, Fraction(1, 1 << _rank(clean_x))))

    def test_two_clean_exact_source_and_strict_quarter_boundary(self):
        z = _tau(2, 1)
        target = scale(add(eye(2), z), Fraction(1, 4))
        exact = _gram(3, 1, [('CX', 0, 1), ('H', 2)])
        self.assertEqual(exact, target)
        blind = _gram(3, 1, [('H', 1), ('H', 2)])
        self.assertEqual(blind, scale(eye(2), Fraction(1, 4)))
        defect = add(blind, scale(target, -1))
        self.assertEqual(multiply(defect, defect), scale(eye(2), Fraction(1, 16)))
        # A distinct Pauli half-projector has the larger sqrt(2)/4 gap.
        other = _gram(3, 1, [('H', 0), ('CX', 0, 1), ('H', 2)])
        difference = add(other, scale(target, -1))
        self.assertEqual(multiply(difference, difference),
                         scale(eye(2), Fraction(1, 8)))

    def test_one_clean_extreme_grams_retain_the_half_error_floor(self):
        target = scale(add(eye(2), _tau(2, 1)), Fraction(1, 4))
        expected = [scale(eye(2), 0), eye(2), scale(eye(2), Fraction(1, 2)),
                    scale(add(eye(2), _tau(2, 1)), Fraction(1, 2))]
        words = [[('X', 1)], [], [('H', 1)], [('CX', 0, 1)]]
        for word, value in zip(words, expected):
            gram = _gram(2, 1, word)
            self.assertEqual(gram, value)
            diagonal = [gram[j][j]-target[j][j] for j in range(2)]
            self.assertTrue(any(entry == CO/2 or entry == -CO/2 for entry in diagonal)
                            or gram == eye(2))

    def test_two_signal_column_or_recovers_every_rank_one_kernel(self):
        # Exhaust every nonzero rank-one 2x2 M and every data-label row c.
        for columns in product(range(4), repeat=2):
            if _rank(columns) != 1:
                continue
            w = [int(column != 0) for column in columns]
            kernel_vector = w[1] | (w[0] << 1)
            self.assertNotEqual(kernel_vector, 0)
            self.assertEqual(_xor_selected(columns, kernel_vector), 0)
            rowspace = {0}
            for row in range(2):
                rowspace |= {((columns[0] >> row) & 1)
                             | (((columns[1] >> row) & 1) << 1)}
            for c in range(4):
                coordinate = (c & kernel_vector).bit_count() & 1
                self.assertEqual(coordinate, int(c not in rowspace))
                # OR(a,b)=a+b+ab, so each kernel entry has degree <=2D.
                ors = []
                for column in columns:
                    a, b = column & 1, column >> 1
                    ors.append(a ^ b ^ (a & b))
                self.assertEqual(ors, w)


if __name__ == '__main__':
    unittest.main()
