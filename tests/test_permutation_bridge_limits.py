"""Finite witnesses for two carefully scoped permutation-bridge claims.

The four-shear fixture exhausts the permitted clean/dirty basis inputs.
The encoded-Clifford fixture uses a dense spectral isometry whose native
implementation is not supplied.  These are algebra checks, not general
proofs, resource estimates, or scalable permutation synthesis routines.
"""
from __future__ import annotations

import unittest

import numpy as np


ATOL = 3e-12


def _four_shears(data, work):
    """Chronological A_f, B_g, A_h, B_k; first data bit is least significant."""
    a, b, c = data
    d, e = work
    # A_f: only Y changes, with f(X)=((1+a)(1+b), 1+a).
    d ^= (1 ^ a) & (1 ^ b)
    e ^= 1 ^ a
    # B_g: only X changes, with g(Y)=(1, Y_2, Y_1).
    a ^= 1
    b ^= e
    c ^= d
    # A_h reads the current data, with h(X)=(ab, a).
    d ^= a & b
    e ^= a
    # B_k has k(Y)=(0, 0, Y_1).
    c ^= d
    return (a, b, c), (d, e)


def _word(bits):
    return sum(bit << index for index, bit in enumerate(bits))


def _cycle_permutation(size):
    """Independent target: ordinary addition by one on computational labels."""
    result = np.zeros((size, size), dtype=complex)
    for old in range(size):
        result[(old + 1) % size, old] = 1
    return result


def _binary_multiplier(width, polynomial):
    """Basis permutation for multiplication by x modulo a binary polynomial."""
    dimension = 1 << width
    columns = []
    for bit in range(width):
        image = 1 << (bit + 1)
        if image & dimension:
            image ^= polynomial
        columns.append(image)
    mapping = []
    for old in range(dimension):
        new = 0
        for bit, image in enumerate(columns):
            if (old >> bit) & 1:
                new ^= image
        mapping.append(new)
    result = np.zeros((dimension, dimension), dtype=complex)
    for old, new in enumerate(mapping):
        result[new, old] = 1
    return result, mapping


class PermutationBridgeLimitsTests(unittest.TestCase):
    def test_one_clean_four_shear_order_eight_is_sharp(self):
        actual = np.zeros((32, 16), dtype=complex)
        expected = np.zeros_like(actual)
        for value in range(8):
            for dirty in range(2):
                data = tuple((value >> bit) & 1 for bit in range(3))
                returned_data, returned_work = _four_shears(data, (dirty, 0))
                self.assertEqual(returned_work, (dirty, 0))
                # Arithmetic decrement is independent of the shear formulas.
                target = (value - 1) % 8
                self.assertEqual(_word(returned_data), target)
                column = 2 * value + dirty
                actual[4 * _word(returned_data) + _word(returned_work), column] = 1
                expected[4 * target + dirty, column] = 1
        # All permitted columns, including both dirty values, agree literally.
        np.testing.assert_array_equal(actual, expected)

        value, visited = 0, set()
        for _ in range(8):
            visited.add(value)
            bits = tuple((value >> bit) & 1 for bit in range(3))
            value = _word(_four_shears(bits, (1, 0))[0])
        self.assertEqual(value, 0)
        self.assertEqual(visited, set(range(8)))

        # A four-bit triangular increment has order sixteen, exceeding eight.
        increment = [(value + 1) % 16 for value in range(16)]
        for value in range(16):
            for bit in range(4):
                lower = (1 << bit) - 1
                carry = int(value & lower == lower)
                self.assertEqual((increment[value] >> bit) & 1,
                                 ((value >> bit) & 1) ^ carry)
            moved = value
            for _ in range(8):
                moved = increment[moved]
            self.assertNotEqual(moved, value)
            for _ in range(8):
                moved = increment[moved]
            self.assertEqual(moved, value)

    def test_approximate_cycle_embedding_in_exact_linear_clifford(self):
        # x^3+x+1, x^4+x+1, and x^5+x^2+1 are checked by their full orbits.
        polynomials = {3: 0b1011, 4: 0b10011, 5: 0b100101}
        for logical_width in (2, 3, 4):
            with self.subTest(logical_width=logical_width):
                size = 1 << logical_width
                physical_width = logical_width + 1
                period = 2 * size - 1
                clifford, mapping = _binary_multiplier(
                    physical_width, polynomials[physical_width])
                self.assertEqual(mapping[0], 0)
                self.assertEqual(len(set(mapping)), 2 * size)
                orbit, value = [], 1
                for _ in range(period):
                    orbit.append(value)
                    value = mapping[value]
                self.assertEqual(value, 1)
                self.assertEqual(set(orbit), set(range(1, 2 * size)))
                # The explicit invertible binary map is a CNOT/SWAP Clifford.
                for left in range(2 * size):
                    for right in range(2 * size):
                        self.assertEqual(mapping[left ^ right],
                                         mapping[left] ^ mapping[right])

                target = _cycle_permutation(size)
                indices = np.arange(size)
                logical_fourier = np.exp(
                    -2j * np.pi * np.outer(indices, indices) / size) / np.sqrt(size)
                # Integer nearest-root rounding, including the half-turn tie.
                roots = [(2 * period * k + size) // (2 * size) % period
                         for k in range(size)]
                self.assertEqual(len(set(roots)), size)
                physical_fourier = np.zeros((2 * size, size), dtype=complex)
                for column, root in enumerate(roots):
                    physical_fourier[orbit, column] = np.exp(
                        -2j * np.pi * root * np.arange(period) / period) / np.sqrt(period)
                embedding = physical_fourier @ logical_fourier.conj().T
                np.testing.assert_allclose(embedding.conj().T @ embedding,
                                           np.eye(size), atol=ATOL, rtol=0)
                error = np.linalg.norm(clifford @ embedding - embedding @ target, 2)
                self.assertLessEqual(error, np.pi / period + ATOL)
                self.assertAlmostEqual(error, 2 * np.sin(np.pi / (2 * period)),
                                       delta=ATOL)
                # Approximate intertwining does not give exact finite order.
                np.testing.assert_allclose(np.linalg.matrix_power(clifford, period),
                                           np.eye(2 * size), atol=ATOL, rtol=0)
                self.assertGreater(np.linalg.norm(
                    np.linalg.matrix_power(target, period) - np.eye(size), 2), 1)


if __name__ == "__main__":
    unittest.main()
