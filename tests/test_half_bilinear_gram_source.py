"""Exact guards for a normalized, but wrong-sector, dense half-phase trial.

The p=2 fixture uses integer Walsh numerators. It is finite evidence,
not an all-size source exclusion or a scalable native emitter.
"""
from fractions import Fraction
import unittest

import numpy as np

from tests.exact_matrices import CI, CO, CQ2, CZ, Q2, adjoint, native_columns
from tests.test_gram_four_call_readout import _partial_gram_source
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts
from tests.test_nilpotent_pencil_observable_source import _field_matrix, _power


def _private_toffoli(a, b, target, phase):
    controls = (a, b, target)
    vertices = (0, 1, 3, 2, 6, 7, 5, 4)
    word, old = [('H', target)], 0
    for vertex in vertices:
        difference = old ^ vertex
        if difference:
            word.append(('CX', controls[difference.bit_length() - 1], phase))
        word.append(('TDG' if vertex.bit_count() % 2 else 'T', phase))
        old = vertex
    word.extend([('CX', target, phase), ('H', target)])
    return word


def _ccs_word():
    a, b, c, d, e, phase = range(6)
    td = _private_toffoli(a, b, d, phase)
    te = _private_toffoli(a, b, e, phase)
    l = [('T', d), ('CX', c, d), ('TDG', d), ('CX', c, d)]
    ccz = [('H', e)] + _private_toffoli(c, d, e, phase) + [('H', e)]
    return td + l + td + _inverse(l) + te + ccz + te + ccz


def _binary_image(rows, value):
    return sum(((row & value).bit_count() % 2) << bit
               for bit, row in enumerate(rows))


def _walsh_integer(values):
    distance = 1
    while distance < len(values):
        view = values.reshape(-1, 2 * distance)
        left, right = view[:, :distance].copy(), view[:, distance:].copy()
        view[:, :distance] = left + right
        view[:, distance:] = left - right
        distance *= 2


def _probe_expectation(p, left_sign=1, right_sign=1):
    """Exact p<=2 carrier expectation at field x=1, z=e_p, t=0.

    Four full Hadamard walls give denominator 2^(4n). Intermediate
    numerators have magnitude <=2^(4n)<=2^32, safely inside int64.
    Convert to Python integers BEFORE squaring for the exact receipt.
    """
    m, n = 1 << p, p * (1 << p)
    size, dimension = 1 << n, 1 << (2 * n)
    basis = ((3, 2) if p == 1 else
             _field_matrix([[_power(u, j) for u in range(m)] for j in range(m)]))
    images = np.array([_binary_image(basis, z) for z in range(size)], dtype=np.int64)
    labels = np.arange(dimension, dtype=np.int64)
    z, t = labels & (size - 1), labels >> n
    permutation = (images[z] ^ t) | (z << n)
    exponent = np.zeros(dimension, dtype=np.int64)
    # At x=1, multiplication by x is the identity in the binary basis.
    for block in range(m - 1):
        for bit in range(p):
            exponent += (((z >> (block * p + bit)) & 1)
                         * ((t >> ((block + 1) * p + bit)) & 1))
    real, imag = np.zeros(dimension, dtype=np.int64), np.zeros(dimension, dtype=np.int64)
    real[1 << p] = 1

    def phase(sign):
        old_real, old_imag = real.copy(), imag.copy()
        powers = (sign * exponent) % 4
        for power in (1, 2, 3):
            mask = powers == power
            if power == 1:
                real[mask], imag[mask] = -old_imag[mask], old_real[mask]
            elif power == 2:
                real[mask], imag[mask] = -old_real[mask], -old_imag[mask]
            else:
                real[mask], imag[mask] = old_imag[mask], -old_real[mask]

    def mixer(inverse=False):
        nonlocal real, imag
        if not inverse:
            new_real, new_imag = np.zeros_like(real), np.zeros_like(imag)
            new_real[permutation], new_imag[permutation] = real, imag
            real, imag = new_real, new_imag
        _walsh_integer(real)
        _walsh_integer(imag)
        if inverse:
            real, imag = real[permutation], imag[permutation]

    # V=D^-left E^-right D^left E^right; E=R† D R.
    mixer()
    phase(right_sign)
    mixer(True)
    phase(left_sign)
    mixer()
    phase(-right_sign)
    mixer(True)
    phase(-left_sign)
    squares = [int(a) ** 2 + int(b) ** 2 for a, b in zip(real, imag)]
    denominator = 1 << (8 * n)
    norm = Fraction(sum(squares), denominator)
    numerator = sum((-1 if (int(zi) ^ int(ti)) & 1 else 1) * square
                    for zi, ti, square in zip(z, t, squares))
    return norm, Fraction(numerator, denominator)


class HalfBilinearGramSourceTests(unittest.TestCase):
    def test_base_source_is_literal_h_or_x_on_both_signal_inputs(self):
        word = _partial_gram_source()
        self.assertEqual(_counts(word), (2, 19))
        self.assertEqual(_counts(_inverse(word)), (2, 19))
        actual = native_columns(3, word)
        inverse_root = CQ2(Q2(0, Fraction(1, 2)))
        for column in range(8):
            equal = (column & 1) == ((column >> 1) & 1)
            signal = (column >> 2) & 1
            for row in range(8):
                expected = CZ
                if row & 3 == column & 3:
                    if equal:
                        expected = inverse_root * (-1) ** (signal * ((row >> 2) & 1))
                    elif row == column ^ 4:
                        expected = CO
                self.assertEqual(actual[row][column], expected)

    def test_literal_52_t_ccs_returns_every_dirty_input(self):
        word = _ccs_word()
        self.assertEqual(_counts(word), (52, 68))
        self.assertTrue(all(gate[1] in (3, 4, 5) for gate in word
                            if gate[0] in ('T', 'TDG')))
        actual = native_columns(6, word)
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))
        for column in range(64):
            scalar = CI if column & 7 == 7 else CO
            for row in range(64):
                self.assertEqual(actual[row][column], scalar if row == column else CZ)

    def test_binary_affine_control_keeps_arbitrary_phase_helper_return(self):
        word = []
        for gate in _ccs_word():
            if gate[0] == 'CX' and gate[1] == 0:
                word.extend([gate, ('CX', 6, gate[2])])
            else:
                word.append(gate)
        self.assertEqual(_counts(word), (52, 84))
        actual = native_columns(7, word)
        self.assertEqual(native_columns(7, _inverse(word)), adjoint(actual))
        for column in range(128):
            affine = ((column >> 0) ^ (column >> 6)) & 1
            scalar = CI if affine and column & 6 == 6 else CO
            for row in range(128):
                self.assertEqual(actual[row][column], scalar if row == column else CZ)

    def test_dense_commutator_p1_exact_sector_failure_and_orientations(self):
        orientations = ((1, 1), (1, -1), (-1, 1), (-1, -1),
                        (2, 1), (1, 2), (2, 2))
        for left, right in orientations:
            norm, expectation = _probe_expectation(1, left, right)
            self.assertEqual(norm, 1)
            self.assertEqual(expectation, 1)
            # At x=1,z=e_1,t=0 the desired paired lookup character is -1.
            self.assertEqual(abs(expectation + 1) / 4, Fraction(1, 2))

    def test_dense_commutator_p2_exact_walsh_receipts(self):
        expected = {(1, 1): Fraction(19, 64), (1, -1): Fraction(39, 128),
                    (-1, 1): Fraction(39, 128), (-1, -1): Fraction(19, 64),
                    (2, 1): Fraction(3, 32), (1, 2): Fraction(1, 4),
                    (2, 2): Fraction(0)}
        for orientations, value in expected.items():
            norm, expectation = _probe_expectation(2, *orientations)
            self.assertEqual(norm, 1)
            self.assertEqual(expectation, value)
            self.assertGreaterEqual(abs(expectation + 1) / 4, Fraction(1, 4))
        self.assertEqual((1 + expected[(1, 1)]) / 4, Fraction(83, 256))


if __name__ == '__main__':
    unittest.main()
