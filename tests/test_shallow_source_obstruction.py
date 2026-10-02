"""Bounded checks of the full-input two-T-layer source obstruction.

Native matrices audit all Pauli-transfer entries on two and three wires.
Exact subset sums and high-precision gap enumeration independently check
the diagonal witnesses; these finite checks do not prove a depth bound.
"""
from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import product
import unittest

import numpy as np

from tests.test_operator_source_compiler import (
    ATOL, I2, X, Y, Z, _adjoint, _optimized_source_word, _pauli,
    _source_word, _word_matrix,
)


def _pauli_basis(width):
    return np.array([_pauli(width, dict(enumerate(factors)))
                     for factors in product((I2, X, Y, Z), repeat=width)])


def _clifford_word(width, rng):
    # Include unrestricted H/CX blocks, rather than Majorana permutations.
    gates = [('H', 0), ('CX', 0, width - 1)]
    for _ in range(5 * width):
        if rng.integers(3) == 0:
            control, target = rng.choice(width, size=2, replace=False)
            gates.append(('CX', int(control), int(target)))
        else:
            gates.append((('H', 'S', 'SDG', 'X')[rng.integers(4)],
                          int(rng.integers(width))))
    return gates


def _weights(width):
    return [Fraction(1, 2 ** min(j + 1, width - 1))
            for j in range(width)]


def _decimal(value):
    return Decimal(value.numerator) / Decimal(value.denominator)


def _enumerated_gap(grid, precision_bits):
    # The alphabet includes the grid's smallest positive point. Smaller
    # alphabet points cannot improve its distance at any positive grid point.
    root_two = Decimal(2).sqrt()
    alphabet = [Decimal(0)] + [root_two ** (-j)
                               for j in range(2 * precision_bits + 3)]
    return max(min(abs(_decimal(abs(value)) - point) for point in alphabet)
               for value in grid)


def _midpoint_gap(spacing):
    inverse_root = 1 / Decimal(2).sqrt()
    midpoint = (1 + inverse_root) / 2
    step = _decimal(spacing)
    nearest = step * (midpoint / step).to_integral_value()
    return (1 - inverse_root) / 2 - abs(nearest - midpoint)


class ShallowSourceObstructionTests(unittest.TestCase):
    def test_native_two_layer_full_transfer_alphabet(self):
        rng = np.random.default_rng(17029)
        for width in (2, 3):
            dimension = 1 << width
            paulis = _pauli_basis(width)
            alphabet = np.array([0.0] + [sign * 2 ** (-j / 2)
                                for sign in (-1, 1)
                                for j in range(2 * width + 1)])
            for depth in (0, 1, 2):
                for sample in range(8):
                    word = _clifford_word(width, rng)
                    for _ in range(depth):
                        exponents = rng.integers(-1, 2, size=width)
                        if not np.any(exponents):
                            exponents[0] = 1
                        layer = [('T' if exponent == 1 else 'TDG', wire)
                                 for wire, exponent in enumerate(exponents)
                                 if exponent]
                        self.assertEqual(len(layer), len({gate[1] for gate in layer}))
                        word += layer + _clifford_word(width, rng)
                    unitary = _word_matrix(width, word)
                    transformed = unitary @ paulis @ unitary.conj().T
                    transfer = np.einsum('aij,bji->ab', paulis, transformed) / dimension
                    with self.subTest(width=width, depth=depth, sample=sample):
                        self.assertLess(np.max(np.abs(transfer.imag)), ATOL)
                        distances = np.min(np.abs(transfer.real[..., None] - alphabet),
                                           axis=-1)
                        self.assertLess(np.max(distances), ATOL)
                        np.testing.assert_allclose(transfer.real @ transfer.real.T,
                                                   np.eye(len(paulis)),
                                                   atol=ATOL, rtol=0)

    def test_exact_subset_grids_and_optimized_gap(self):
        with localcontext() as context:
            context.prec = 80
            tolerance = Decimal('1e-65')
            for width in range(2, 10):
                sums = {Fraction(0)}
                for weight in _weights(width):
                    sums |= {value + weight for value in sums}
                source_grid = {1 - 2 * value for value in sums}
                source_spacing = Fraction(1, 2 ** (width - 2))
                expected_source = {-1 + index * source_spacing
                                   for index in range(2 ** (width - 1) + 1)}
                self.assertEqual(source_grid, expected_source)
                controlled_grid = {(1 + value) / 2 for value in source_grid}
                controlled_spacing = source_spacing / 2
                self.assertEqual(controlled_grid,
                                 {index * controlled_spacing
                                  for index in range(2 ** (width - 1) + 1)})
                for controlled, grid, spacing in (
                        (False, source_grid, source_spacing),
                        (True, controlled_grid, controlled_spacing)):
                    actual = _enumerated_gap(grid, width)
                    predicted = _midpoint_gap(spacing)
                    with self.subTest(width=width, controlled=controlled):
                        self.assertLess(abs(actual - predicted), tolerance)
                        if width >= (4 if controlled else 5):
                            self.assertGreaterEqual(actual + tolerance, Decimal(1) / 8)
                            limit = (1 - 1 / Decimal(2).sqrt()) / 4
                            self.assertLessEqual(actual / 2, limit + tolerance)
                            self.assertLessEqual(limit - actual / 2,
                                                 _decimal(spacing) / 4 + tolerance)

    def test_native_source_witnesses_with_arbitrary_dirty_extensions(self):
        for width, witness, coefficient in ((3, 1, Fraction(1, 2)),
                                             (4, 2, Fraction(3, 4)),
                                             (5, 3, Fraction(7, 8))):
            gammas = [_pauli(width, {**{k: Z for k in range(j)}, j: X})
                      for j in range(width)]
            expected = sum(np.sqrt(float(weight)) * gamma
                           for weight, gamma in zip(_weights(width), gammas))
            word = _optimized_source_word(width)
            if width == 3:
                # Two singleton T layers attain the first width below exclusion.
                self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 2)
            actual = _word_matrix(width + 1, word)
            np.testing.assert_allclose(actual, np.kron(I2, expected),
                                       atol=ATOL, rtol=0)
            # Include off-diagonal helper Paulis, not only its zero-input sector.
            for helper_pauli in (I2, X, Y, Z):
                pauli = np.kron(helper_pauli, _pauli(width, {witness: Z}))
                transfer = np.trace(pauli @ actual @ pauli @ actual.conj().T) / len(actual)
                self.assertLess(abs(transfer - float(coefficient)), ATOL)

        for width, witness, coefficient in ((2, 0, Fraction(1, 2)),
                                             (3, 1, Fraction(3, 4)),
                                             (4, 2, Fraction(7, 8))):
            loader = _source_word(width)
            word = _adjoint(loader) + [('CX', width, 0)] + loader
            if width == 2:
                self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 2)
            # Core, arbitrary control, and one arbitrary helper: at most 64 columns.
            actual = _word_matrix(width + 2, word)
            source = _word_matrix(width, _optimized_source_word(width))
            identity, zero = np.eye(1 << width), np.zeros_like(source)
            expected = np.block([[identity, zero], [zero, source]])
            np.testing.assert_allclose(actual, np.kron(I2, expected),
                                       atol=ATOL, rtol=0)
            pauli = _pauli(width + 2, {witness: Z, width + 1: X})
            transfer = np.trace(pauli @ actual @ pauli @ actual.conj().T) / len(actual)
            self.assertLess(abs(transfer - float(coefficient)), ATOL)


if __name__ == '__main__':
    unittest.main()
