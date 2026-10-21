"""Exact guards for corank-one Fourier/Bell acceptance and its failed repair.

The p=2 Walsh propagation uses complete arbitrary-bank unitaries, not a
normalized trace.  These receipts challenge phases, transposes and success
weights; they do not provide the missing one-clean central-character source.
"""
from fractions import Fraction
import unittest

import numpy as np

from tests.exact_matrices import CO, CZ, CQ2, adjoint, multiply, native_columns
from tests.test_bent_gram_source import ry_t
from tests.test_half_bilinear_gram_source import (
    _binary_image, _private_toffoli, _walsh_integer,
)
from tests.test_majorana_walsh_query import _inverse as word_inverse
from tests.test_nilpotent_pencil_observable_source import (
    _field_matrix, _power, _transpose,
)
from tests.test_phase_signal_query_limits import _counts


def _kernel(p, x):
    """Independent binary expansion of (I+xL)^T P and its left null vector."""
    m, n = 1 << p, p * (1 << p)
    if p == 1:
        pencil = (1, x | 2)
        vector = 1 | (x << 1)
    else:
        assert p == 2
        pencil = _field_matrix([
            [int(i == j) ^ (x if i == j + 1 else 0) for j in range(m)]
            for i in range(m)])
        vector = sum(_power(x, j) << (p*j) for j in range(m))
    return tuple(row & ~1 for row in _transpose(pencil)), vector, n


def _swap(a, b):
    return [('CX', a, b), ('CX', b, a), ('CX', a, b)]


def _w_word(rows, data, auxiliary):
    word = [('CX', data[j], auxiliary[i])
            for i, row in enumerate(rows) for j in range(len(rows))
            if (row >> j) & 1]
    word += [('H', wire) for wire in data + auxiliary]
    for a, b in zip(data, auxiliary):
        word += _swap(a, b)
    return word


def _reverse_bell_outcomes(p, x, a, sign):
    """Integer propagation of W^dagger on (|0>+sign|v>)/sqrt(2) tensor |a>.

    One full 2n-qubit Walsh wall has denominator q; this function retains
    integer amplitudes and divides squared norms by 2q^2 only at the end.
    The shear is applied as a permutation after the actual inverse walls.
    """
    rows, vector, n = _kernel(p, x)
    q = 1 << n
    values = np.zeros(q*q, dtype=np.int64)
    # SWAP is the first chronological operation in W^dagger.
    values[a] = 1
    values[a | (vector << n)] = sign
    _walsh_integer(values)
    images = np.array([_binary_image(rows, z) for z in range(q)], dtype=np.int64)
    labels = np.arange(q*q, dtype=np.int64)
    permutation = (labels & (q-1)) | (((labels >> n) ^ images[labels & (q-1)]) << n)
    output = np.zeros_like(values)
    output[permutation] = values
    numerators = (output.reshape(q, q) ** 2).sum(axis=1)
    return [Fraction(int(value), 2*q*q) for value in numerators]


class CorankFourierBellSourceTests(unittest.TestCase):
    def test_p1_p2_full_fourier_grams_fix_the_two_kernel_orientations(self):
        for p in (1, 2):
            for x in range(1 << p):
                rows, vector, n = _kernel(p, x)
                q = 1 << n
                images = [_binary_image(rows, z) for z in range(q)]
                self.assertEqual([z for z, image in enumerate(images) if not image], [0, 1])
                self.assertEqual(_binary_image(_transpose(rows), vector), 0)
                self.assertEqual(len(set(images)), q // 2)
                fourier = np.array([
                    [1 - 2*((b & image).bit_count() & 1) for image in images]
                    for b in range(q)], dtype=np.int64)
                # Native normalization is 1/sqrt(q), not 1/q or a trace.
                right, left = fourier.T @ fourier, fourier @ fourier.T
                expected_right = np.zeros((q, q), dtype=np.int64)
                expected_left = np.zeros((q, q), dtype=np.int64)
                for z in range(q):
                    expected_right[z, z] = expected_right[z ^ 1, z] = q
                    expected_left[z, z] = expected_left[z ^ vector, z] = q
                np.testing.assert_array_equal(right, expected_right)
                np.testing.assert_array_equal(left, expected_left)
                if x:
                    self.assertNotEqual(vector, 1)
                    self.assertFalse(np.array_equal(right, left))

    def test_literal_p1_bell_word_all_address_bank_and_dirty_columns(self):
        # W has x=0, data=1,2, auxiliary=3,4, arbitrary phase helper=5.
        query = _private_toffoli(0, 2, 3, 5) + [('CX', 2, 4)]
        word = query + [('H', wire) for wire in (1, 2, 3, 4)]
        word += _swap(1, 3) + _swap(2, 4)
        self.assertEqual(_counts(word), (8, 21))
        actual = native_columns(6, word)
        for column in range(64):
            x, data, auxiliary, helper = column & 1, (column >> 1) & 3, (column >> 3) & 3, column >> 5
            rows, _, _ = _kernel(1, x)
            shifted = auxiliary ^ _binary_image(rows, data)
            for row in range(64):
                b, a = (row >> 1) & 3, (row >> 3) & 3
                phase = -1 if ((a & data).bit_count() + (b & shifted).bit_count()) & 1 else 1
                expected = CQ2(Fraction(phase, 4)) if row & 1 == x and row >> 5 == helper else CZ
                self.assertEqual(actual[row][column], expected)
        self.assertEqual(native_columns(6, word_inverse(word)), adjoint(actual))

    def test_shared_outcome_weight_from_complete_inverse_walsh_propagation(self):
        for p in (1, 2):
            n = p * (1 << p)
            q = 1 << n
            for x in range(1 << p):
                plus = _reverse_bell_outcomes(p, x, a=q-1, sign=1)
                minus = _reverse_bell_outcomes(p, x, a=q//2, sign=-1)
                self.assertEqual(sum(plus), 1)
                self.assertEqual(sum(minus), 1)
                # Equality of the two complete outcome registers, with
                # every accepted outcome kept orthogonal, is not a trace.
                same_sector = sum(value*value for value in plus)
                wrong_sector = sum(a*b for a, b in zip(plus, minus))
                self.assertEqual(same_sector, Fraction(2, q))
                self.assertEqual(wrong_sector, 0)
                self.assertEqual(same_sector, Fraction(1, 2) if p == 1 else Fraction(1, 128))

    def test_coherent_swap_complete_native_action_and_wrong_sector(self):
        # Fixed x=1, p=1; use all basis and two-bit probes to challenge
        # the phase as well as the transpose in the conjugated swap.
        rows, vector, n = _kernel(1, 1)
        banks = (list(range(0, 2)), list(range(2, 4)),
                 list(range(4, 6)), list(range(6, 8)))
        pair = (_w_word(rows, banks[0], banks[2])
                + _w_word(rows, banks[1], banks[3]))
        auxiliary_swap = _swap(4, 6) + _swap(5, 7)
        word = word_inverse(pair) + auxiliary_swap + pair
        columns = sorted({0, vector | (vector << n)}
                         | {1 << i for i in range(8)}
                         | {(1 << i) | (1 << j) for i in range(8) for j in range(i)})
        actual = native_columns(8, word, columns=columns)
        for k, column in enumerate(columns):
            c, d, y, z = ((column >> (n*j)) & 3 for j in range(4))
            shift = _binary_image(_transpose(rows), c ^ d)
            target = d | (c << n) | ((y ^ shift) << (2*n)) | ((z ^ shift) << (3*n))
            for row in range(256):
                self.assertEqual(actual[row][k], CO if row == target else CZ)
        # The two components of the wrong-sector witness are both fixed.
        for column in (0, vector | (vector << n)):
            self.assertEqual(actual[column][columns.index(column)], CO)
        self.assertNotEqual(_binary_image(rows, 1),
                            _binary_image(_transpose(rows), 1))

    def test_full_equality_gram_retains_off_diagonal_auxiliary_blocks(self):
        # Derive the whole four-bank operator from two literal native
        # W matrices and an actual output equality projection.  Integer
        # numerators retain every auxiliary off-diagonal block.
        p, n, q = 1, 2, 4
        labels = np.arange(q**4)
        c, d, y, z = [(labels >> (n*j)) & (q-1) for j in range(4)]
        left_labels, right_labels = c | (y << n), d | (z << n)
        equal_columns = labels[y == z]
        for x in range(1 << p):
            rows, vector, _ = _kernel(p, x)
            one = native_columns(2*n, _w_word(rows, [0, 1], [2, 3]))
            signs = np.zeros((q*q, q*q), dtype=np.int64)
            for i, row in enumerate(one):
                for j, value in enumerate(row):
                    self.assertIn(value, (CQ2(Fraction(1, q)), CQ2(Fraction(-1, q))))
                    signs[i, j] = 1 if value == CQ2(Fraction(1, q)) else -1
            pair = (signs[left_labels[:, None], left_labels[None, :]]
                    * signs[right_labels[:, None], right_labels[None, :]])
            actual = pair[:, equal_columns] @ pair[:, equal_columns].T
            expected = np.zeros_like(actual)
            for column in range(q**4):
                for ell in range(q):
                    shift = _binary_image(_transpose(rows), ell)
                    output = column ^ ell ^ (ell << n) ^ (shift << (2*n)) ^ (shift << (3*n))
                    expected[output, column] += q**3
            np.testing.assert_array_equal(actual, expected)
            # Denominator q^4: rank q^3, and the desired charge is one
            # stabilizer among n, rather than the entire projector.
            self.assertEqual(int(np.trace(actual)), q**7)
            np.testing.assert_array_equal(actual,
                                          actual[labels ^ vector ^ (vector << n)])
            moved_auxiliary = next(ell for ell in range(q)
                                   if _binary_image(_transpose(rows), ell))
            shift = _binary_image(_transpose(rows), moved_auxiliary)
            output = moved_auxiliary | (moved_auxiliary << n) | (shift << (2*n)) | (shift << (3*n))
            self.assertEqual(actual[output, 0], q**3)

    def test_literal_normalized_signal_fredkin_has_correct_phase_and_inverse(self):
        # Data0,1, clean signal2, arbitrary returned mask3.  This tests
        # the actual controlled-SWAP primitive, not an abstract reflection.
        fredkin = ([('CX', 0, 1)] + _private_toffoli(2, 1, 0, 3)
                   + [('CX', 0, 1)])
        word = ry_t(2, inverse=True) + fredkin + ry_t(2, cube=True)
        self.assertEqual(_counts(word), (10, 25))
        actual = native_columns(4, word)
        self.assertEqual(native_columns(4, word_inverse(word)), adjoint(actual))
        clean = [index for index in range(16) if not index & 4]
        block = [[actual[row][column] for column in clean] for row in clean]
        gram = multiply(adjoint(block), block)
        alpha = CQ2(Fraction(1, 4), Fraction(1, 4))
        for j, column in enumerate(clean):
            swapped = (column & 8) | ((column & 1) << 1) | ((column & 2) >> 1)
            for i, row in enumerate(clean):
                coefficient = int(row == column) + int(row == swapped)
                self.assertEqual(block[i][j], coefficient * alpha)
                self.assertEqual(gram[i][j], CQ2(Fraction(coefficient, 4)))


if __name__ == '__main__':
    unittest.main()
