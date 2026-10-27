"""Exact quotient-label and interference checks for the copy-axis span bound."""
from itertools import product
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, adjoint, multiply, native_columns,
)


def _span(generators):
    result = {0}
    for generator in generators:
        result |= {value ^ generator for value in tuple(result)}
    return result


def _pauli(width, x, z):
    size = 1 << width
    phase = (CO, CI, -CO, -CI)[(x & z).bit_count() % 4]
    return [[phase * (-1 if (z & column).bit_count() % 2 else 1)
             if row == column ^ x else CZ
             for column in range(size)] for row in range(size)]


def _coefficients(matrix, width):
    size = 1 << width
    result = {}
    for x, z in product(range(size), repeat=2):
        pauli = _pauli(width, x, z)
        result[x, z] = sum((pauli[row][column].conj() * matrix[row][column]
                           for row in range(size) for column in range(size)),
                          CZ) / size
    return result


def _sector_weights(coefficients, output, span):
    weights = {}
    for (x, z), coefficient in coefficients.items():
        label = ((x >> output) & 1) | (((z >> output) & 1) << 1)
        key = min(label ^ value for value in span)
        weights[key] = weights.get(key, 0) + coefficient.abs2()
    return weights


class CopyAxisSpanTests(unittest.TestCase):
    def test_all_two_bit_affine_x_maps_and_one_axis_quotients(self):
        # 65,536 exact cases. The unquotiented maximum is three;
        # a single quotient axis can make all four addresses match.
        maxima = {0: 0, 1: 0}
        for axis in range(16):
            span = _span([axis])
            rank = len(span).bit_length() - 1
            for offset, a, b in product(range(16), repeat=3):
                values = (offset, offset ^ a, offset ^ b, offset ^ a ^ b)
                hits = sum((value ^ (1 << x)) in span
                           for x, value in enumerate(values))
                self.assertLessEqual(hits, min(4, 3 + rank))
                maxima[rank] = max(maxima[rank], hits)
        self.assertEqual(maxima, {0: 3, 1: 4})

    def test_mixed_axes_use_intersection_not_full_projected_rank(self):
        # X labels occupy the low four bits, Z labels the high four.
        # Neither mixed generator is X-only, but their sum is.
        mixed = _span([0x11, 0x12])
        self.assertEqual(mixed & set(range(16)), {0, 3})
        cases = [(_span([0x11]), 0), (mixed, 1), (_span([0x10, 0x20]), 0)]
        for span, rank_x in cases:
            self.assertEqual(len(span & set(range(16))), 1 << rank_x)
            # Include non-X affine coefficient vectors and dependent
            # quotient targets; do not presume quotient independence.
            coefficients = (0, 1, 2, 4, 8, 0x11, 0x12, 0x24, 0x48)
            for offset, a, b in product(coefficients, repeat=3):
                values = (offset, offset ^ a, offset ^ b, offset ^ a ^ b)
                hits = sum((value ^ (1 << x)) in span
                           for x, value in enumerate(values))
                self.assertLessEqual(hits, min(4, 3 + rank_x))

    def test_exact_noncommuting_rotations_preserve_coset_norm(self):
        # Carrier 0, output 1. Two native axes X0 Z1 and Z0 Z1
        # anticommute, but their output span contains only Z1.
        zz = [('CX', 0, 1), ('T', 1), ('CX', 0, 1)]
        xz = [('H', 0)] + zz + [('H', 0)]
        word = xz + zz
        unitary = native_columns(2, word)
        span = {0, 2}
        for x, z in product(range(4), repeat=2):
            pauli = _pauli(2, x, z)
            image = multiply(multiply(unitary, pauli), adjoint(unitary))
            coefficients = _coefficients(image, 2)
            input_coset = (x >> 1) & 1
            for (out_x, _), coefficient in coefficients.items():
                if coefficient != CZ:
                    self.assertEqual((out_x >> 1) & 1, input_coset)
        # A full-space Hermitian unitary has several coefficients in
        # one sector. Those coefficients change/interfere, while its
        # quotient-sector squared norms remain exactly unchanged.
        reflection = native_columns(2, [('H', 0), ('H', 1)])
        before = _coefficients(reflection, 2)
        image = multiply(multiply(unitary, reflection), adjoint(unitary))
        after = _coefficients(image, 2)
        self.assertNotEqual(before, after)
        self.assertEqual(_sector_weights(before, 1, span),
                         _sector_weights(after, 1, span))
        restored = multiply(multiply(adjoint(unitary), image), unitary)
        self.assertEqual(restored, reflection)


if __name__ == '__main__':
    unittest.main()
