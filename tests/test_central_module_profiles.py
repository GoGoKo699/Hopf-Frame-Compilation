"""Small exact native-band guards for central-module coefficient profiles.

These check literal full-space norm identities, not a growing original-pair
Boolean family or an asymptotic matching theorem. No input is initialized.
"""
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _conjugate, _tau


def _coefficient(basis, matrix):
    size = len(matrix)
    return sum((basis[row][column].conj() * matrix[row][column]
                for row in range(size) for column in range(size)), CZ) / size


def _grouped_native_image(width, central_generators, orientations):
    """Return exact central-character coefficients in a fixed Pauli module.

    The independent native axes are Z_i and the probe is X on every wire.
    Representatives are P Z^a. Central factors are moved to their left,
    retaining the sign even when P anticommutes with R.
    """
    size = 1 << width
    probe = _tau(size - 1, width)
    band = [('T' if sign == 1 else 'TDG', wire)
            for wire, sign in enumerate(orientations)]
    x_word = [gate for wire in range(width)
              for gate in [('H', wire), ('S', wire), ('S', wire), ('H', wire)]]
    word = _inverse(band) + x_word + band
    actual = native_columns(width, word)

    coefficients = []
    reconstruction = zeros(size)
    for mask in range(size):
        basis = multiply(probe, _tau(mask << width, width))
        coefficient = _coefficient(basis, actual)
        coefficients.append(coefficient)
        reconstruction = add(reconstruction, scale(basis, coefficient))

    central_coordinates = {}
    for coordinates in range(1 << len(central_generators)):
        mask = 0
        for bit, generator in enumerate(central_generators):
            if coordinates >> bit & 1:
                mask ^= generator
        if mask in central_coordinates:
            raise ValueError('central generators must be independent')
        central_coordinates[mask] = coordinates

    grouped = []
    remaining = set(range(size))
    while remaining:
        representative = min(remaining)
        coset = {representative ^ central for central in central_coordinates}
        remaining -= coset
        values = []
        for character in range(1 << len(central_generators)):
            value = CZ
            for mask in coset:
                central = mask ^ representative
                coordinates = central_coordinates[central]
                # Z^central P=(-1)^|central| P Z^central.
                parity = central.bit_count() + (coordinates & character).bit_count()
                value += coefficients[mask] * (-1 if parity % 2 else 1)
            values.append(value)
        grouped.append(values)
    return actual, word, reconstruction, coefficients, grouped


class CentralModuleProfileTests(unittest.TestCase):
    def test_four_axes_distinguish_all_three_profiles(self):
        width = 4
        generators = (0b0011, 0b0110, 0b1100)
        gram = [sum(((left & right).bit_count() % 2) << bit
                    for bit, right in enumerate(generators))
                for left in generators]
        self.assertEqual(_rank(generators), 3)
        self.assertEqual(_rank(gram), 2)
        # k=4,h=3,r=2: nu=2, delta=1, mu=0, not raw rank k=4.
        self.assertEqual((4 - 2, 4 - 3, 4 - 2 * 3 + 2), (2, 1, 0))
        actual, word, reconstruction, coefficients, grouped = _grouped_native_image(
            width, generators, (1, 1, 1, 1))
        self.assertEqual(actual, reconstruction)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(multiply(adjoint(actual), actual), eye(16))
        self.assertTrue(all(value != CZ for value in coefficients))
        displacements = []
        for bit in range(width):
            x_label = 1 << bit
            y_label = x_label | (1 << (width + bit))
            self.assertEqual(_conjugate(actual, _tau(x_label, width)),
                             _tau(y_label, width))
            displacements.append(x_label ^ y_label)
        self.assertEqual(_rank(displacements), 4)
        self.assertEqual(len(grouped), 2)
        squared = [[value.conj() * value for value in row] for row in grouped]
        # Each coset has supremum magnitude one; at each actual R character
        # exactly one survives. Thus sum-sup=2, sup-sum=1, support=2.
        for row in squared:
            self.assertEqual(set(row), {CZ, CO})
        for character in range(8):
            self.assertEqual(sum((row[character] for row in squared), CZ), CO)
            self.assertEqual(sum(row[character] != CZ for row in squared), 1)

    def test_odd_central_relation_and_negative_orientation(self):
        # R=<Z>, P=X: K_p=<1> is not even-weight. The dot rank is one.
        actual, word, reconstruction, coefficients, grouped = _grouped_native_image(
            1, (1,), (-1,))
        self.assertEqual(_rank((1,)), 1)
        self.assertEqual((1 - 1, 1 - 1, 1 - 2 + 1), (0, 0, 0))
        self.assertEqual(actual, reconstruction)
        self.assertEqual(native_columns(1, _inverse(word)), adjoint(actual))
        self.assertEqual(multiply(adjoint(actual), actual), eye(2))
        self.assertTrue(all(value != CZ for value in coefficients))
        self.assertEqual(len(grouped), 1)
        self.assertEqual([value.conj() * value for value in grouped[0]], [CO, CO])
        self.assertNotEqual(grouped[0][0], grouped[0][1])


if __name__ == '__main__':
    unittest.main()
