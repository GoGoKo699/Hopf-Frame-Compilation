"""Focused checks for the scoped literal-permutation transport proof.

The non-C3 source below uses four native bands, not the three-band prefix
required in the actual six-layer comparison. No returned six-layer family
or unrestricted source theorem is inferred from these finite fixtures.
"""
from __future__ import annotations

from collections import Counter
import math
import unittest

import numpy as np

from tests.test_dirty_indicator_depth import _toffoli_two_layers
from tests.test_five_layer_native_defect import _anf, _degree
from tests.test_one_band_axis_matching import _projective_distance
from tests.test_operator_source_compiler import _word_matrix
from tests.test_outer_axis_normalizer import _image, _tableau, _tau
from tests.test_read_only_three_layer_rigidity import H, S, Z
from tests.test_t_depth import _inverse, _word


def _permutation(images):
    matrix = np.zeros((len(images), len(images)), dtype=complex)
    matrix[images, np.arange(len(images))] = 1
    return matrix


def _candidate_image(z):
    a, b, c = [(z >> j) & 1 for j in range(3)]
    return ((z & ~7) | a | ((b ^ (a & c) ^ (a & b)) << 1)
            | ((c ^ (a & b)) << 2))


def _linear_image(columns, value):
    result = 0
    for j, column in enumerate(columns):
        if value >> j & 1:
            result ^= column
    return result


class PermutationSourceTransport(unittest.TestCase):
    def test_rounding_radius_is_sharp_for_both_boundary_mechanisms(self):
        delta = math.sqrt(2 - math.sqrt(2))
        # A nonmonomial Clifford column meets the bound exactly.
        self.assertAlmostEqual(_projective_distance(Z, H), delta)
        self.assertTrue(np.allclose(abs(H), 1 / math.sqrt(2)))
        # The same monomial permutation with two distinct fourth-root
        # ratios meets the bound after the best common phase adjustment.
        self.assertAlmostEqual(_projective_distance(np.eye(2), S), delta)
        np.testing.assert_allclose(
            np.linalg.norm(np.eye(2) - np.exp(-1j * math.pi / 4) * S, 2),
            delta, atol=2e-14)
        # Different computational permutations have orthogonal columns.
        swapped = np.array([[0, 1], [1, 0]], dtype=complex)
        self.assertAlmostEqual(_projective_distance(np.eye(2), swapped),
                               math.sqrt(2))

    def test_four_band_non_c3_candidate_full_operator_and_paid_inverse(self):
        width = 4
        schedule = (_toffoli_two_layers(0, 1, 2, 3)
                    + _toffoli_two_layers(0, 2, 1, 3))
        expected = _permutation([_candidate_image(z) for z in range(16)])
        np.testing.assert_allclose(_word_matrix(width, _word(schedule)),
                                   expected, atol=5e-11, rtol=0)
        np.testing.assert_allclose(_word_matrix(width, _word(_inverse(schedule))),
                                   expected.conj().T, atol=5e-11, rtol=0)
        for native in (schedule, _inverse(schedule)):
            counts = Counter(gate[0] for gate in _word(native))
            self.assertEqual(counts['T'] + counts['TDG'], 16)
            self.assertEqual(counts['CX'] + counts['H'], 92)
            self.assertEqual(sum(kind == 'T' for kind, _ in native), 4)
            self.assertEqual({wire for gate in _word(native) for wire in gate[1:]},
                             set(range(width)))
        # The map is quadratic and returns the arbitrary helper on every
        # column, so the full matrix also retains helper coherences.
        for bit in range(width):
            table = [(_candidate_image(z) >> bit) & 1 for z in range(16)]
            self.assertLessEqual(_degree(_anf(table)), 2)
        self.assertTrue(all((_candidate_image(z) ^ z) & 8 == 0
                            for z in range(16)))

        # A native Pauli's image is not Clifford, so this source is not C3.
        image_xa = expected @ _tau(1, width) @ expected.conj().T
        with self.assertRaises(ValueError):
            _tableau(image_xa, width)
        nonaffine = []
        for z in range(16):
            a, b, c = [(z >> j) & 1 for j in range(3)]
            nonaffine.append((z & ~7) | (a ^ 1)
                             | ((c ^ (a & b)) << 1)
                             | ((b ^ c ^ (a & c)) << 2))
        np.testing.assert_allclose(image_xa, _permutation(nonaffine), atol=2e-14)

        # These are axes in the two displayed middle bands. Their
        # anticommutation forbids a direct common-Clifford band merge.
        left_axis = _tau((8 | 1 | 2) << width, width)  # Z_d Z_a Z_b
        right_axis = _tau(2 | (8 << width), width)    # Z_d X_b
        np.testing.assert_allclose(left_axis @ right_axis,
                                   -right_axis @ left_axis, atol=2e-14)

    def test_non_c3_source_has_compatible_monomial_pauli_inputs(self):
        # This checks transport compatibility only. F has the four-band
        # implementation above and is not a qualifying three-band source.
        width = 3
        source = _permutation([_candidate_image(z) for z in range(8)])
        labels = [2, 4, 2 << width, 4 << width]
        for label in labels:
            endpoint = source @ _tau(label, width) @ source.conj().T
            _tableau(endpoint, width)
            self.assertTrue(np.all(np.sum(abs(endpoint) > 1e-10, axis=0) == 1))

    def test_relative_phase_and_fixed_probe_sign_retain_integer_carries(self):
        width = 3
        # An invertible affine map with a nonzero displacement and a
        # diagonal Clifford phase exercises both formulas independently
        # of any C3 assumption about an external source.
        columns, translation = [3, 6, 4], 5
        linear = [_linear_image(columns, z) for z in range(8)]
        inverse_linear = [linear.index(z) for z in range(8)]
        inverse_displacement = inverse_linear[translation]
        ell, pairs = (1, 3, 2), ((0, 1), (0, 2), (1, 2))

        def phase(z):
            bits = [(z >> j) & 1 for j in range(width)]
            return (2 * sum(coefficient * bit for coefficient, bit in zip(ell, bits))
                    + 4 * sum(bits[i] * bits[j] for i, j in pairs)) % 8

        unitary = (_permutation([z ^ translation for z in linear])
                   @ np.diag([np.exp(1j * math.pi * phase(z) / 4)
                              for z in range(8)]))
        matrix, _ = _tableau(unitary, width)
        theta_carry_seen = False
        for z in range(8):
            image, sign = _image(unitary, z, width)
            a, c = image & 7, image >> width
            predicted = (4 * sign + 2 * (a & c).bit_count()
                         + 4 * (c & translation).bit_count()) % 8
            self.assertEqual(predicted, phase(z))
            without_integer_carry = (4 * sign + 2 * ((a & c).bit_count() % 2)
                                     + 4 * (c & translation).bit_count()) % 8
            theta_carry_seen |= without_integer_carry != phase(z)
        self.assertTrue(theta_carry_seen)

        sign_carry_seen = False
        for label in range(64):
            a, b = label & 7, label >> width
            incoming = np.array([(label >> j) & 1 for j in range(6)])
            out = matrix @ incoming % 2
            output = sum(int(bit) << j for j, bit in enumerate(out))
            a_out, c_out = output & 7, output >> width
            residual = (phase(inverse_displacement ^ a) - phase(inverse_displacement)
                        + 2 * (a & b).bit_count()
                        + 4 * (b & inverse_displacement).bit_count()
                        - 2 * (a_out & c_out).bit_count()) % 8
            self.assertEqual(residual % 4, 0)
            sign = residual // 4
            np.testing.assert_allclose(unitary @ _tau(label, width) @ unitary.conj().T,
                                       (-1) ** sign * _tau(output, width),
                                       atol=3e-14)
            wrong = (phase(inverse_displacement ^ a) - phase(inverse_displacement)
                     + 2 * (a & b).bit_count()
                     + 4 * (b & inverse_displacement).bit_count()
                     - 2 * ((a_out & c_out).bit_count() % 2)) % 8
            sign_carry_seen |= wrong != residual
        self.assertTrue(sign_carry_seen)


if __name__ == '__main__':
    unittest.main()
