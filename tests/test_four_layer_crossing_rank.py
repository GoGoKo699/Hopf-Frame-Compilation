"""Finite algebra audits for the four-layer crossing-rank obstruction.

These check the symplectic reduction, literal native boundary matrices, and
the abstract scaled-cubic-bias counterexample. They neither certify the
imported polynomial-limit theorem nor realize the abstract bias as a native
four-layer transfer coefficient.
"""
from __future__ import annotations

from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction
import unittest

import numpy as np

from tests.test_operator_source_compiler import ATOL, X, Z, _adjoint, _pauli, _word_matrix


def _binary_basis(rows):
    """Exact row reduction on bit-packed vectors over F2."""
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return list(pivots.values())


def _pairing(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((left >> width) & (right & mask)).bit_count()) % 2


def _pairing_rank(left, right, width):
    rows = [sum(_pairing(first, second, width) << index
                for index, second in enumerate(right)) for first in left]
    return len(_binary_basis(rows))


def _clifford_vector(vector, width, word):
    """Binary Clifford action; Pauli signs do not enter crossing rank."""
    mask = (1 << width) - 1
    x, z = vector & mask, vector >> width
    for gate in word:
        name, wire, *rest = gate
        if name == 'H':
            difference = ((x >> wire) ^ (z >> wire)) & 1
            x ^= difference << wire
            z ^= difference << wire
        elif name == 'S':
            z ^= ((x >> wire) & 1) << wire
        elif name == 'CX':
            target = rest[0]
            x ^= ((x >> wire) & 1) << target
            z ^= ((z >> target) & 1) << wire
        else:
            raise ValueError(name)
    return x | (z << width)


def _simultaneous_coordinates(left, right, width):
    """Construct hyperbolic pairs, then complete their isotropic remainder.

    Exhaustive dual-vector selection is deliberately restricted to these
    four-qubit algebra fixtures; it is not a claimed compiler algorithm.
    """
    left, right = list(left), list(right)
    z_axes, x_axes = [], []
    while True:
        pair = next(((i, j) for i, first in enumerate(left)
                     for j, second in enumerate(right)
                     if _pairing(first, second, width)), None)
        if pair is None:
            break
        i, j = pair
        z_axis, x_axis = left.pop(i), right.pop(j)
        left = [row ^ (z_axis if _pairing(row, x_axis, width) else 0)
                for row in left]
        right = [row ^ (x_axis if _pairing(z_axis, row, width) else 0)
                 for row in right]
        z_axes.append(z_axis)
        x_axes.append(x_axis)
    crossing_rank = len(z_axes)
    z_axes.extend(_binary_basis(left + right))
    for vector in range(1, 1 << (2 * width)):
        if len(z_axes) == width:
            break
        if (all(_pairing(vector, row, width) == 0 for row in z_axes + x_axes)
                and len(_binary_basis(z_axes + [vector])) > len(z_axes)):
            z_axes.append(vector)
    if len(z_axes) != width:
        raise AssertionError('Failed to extend the isotropic subspace.')
    for index in range(crossing_rank, width):
        dual = next(vector for vector in range(1, 1 << (2 * width))
                    if all(_pairing(row, vector, width) == int(i == index)
                           for i, row in enumerate(z_axes))
                    and all(_pairing(vector, row, width) == 0 for row in x_axes))
        x_axes.append(dual)

    def coordinates(vector):
        x = sum(_pairing(row, vector, width) << i for i, row in enumerate(z_axes))
        z = sum(_pairing(row, vector, width) << i for i, row in enumerate(x_axes))
        return x | (z << width)

    return crossing_rank, z_axes, x_axes, coordinates


def _cubic_bias(pair_count, blocks):
    """Enumerate disjoint copies of a sum_j x_j y_j, exactly over F2."""
    block_width = 1 + 2 * pair_count
    monomials = [(1 << (block * block_width))
                 | (1 << (block * block_width + 1 + j))
                 | (1 << (block * block_width + 1 + pair_count + j))
                 for block in range(blocks) for j in range(pair_count)]
    width = block_width * blocks
    sign_sum = sum(1 - 2 * (sum((assignment & term) == term for term in monomials) % 2)
                   for assignment in range(1 << width))
    return Fraction(sign_sum, 1 << width)


class FourLayerCrossingRankTests(unittest.TestCase):
    def test_exact_simultaneous_normal_form_including_common_axes(self):
        width = 4
        x = [1 << i for i in range(width)]
        z = [1 << (width + i) for i in range(width)]
        # (L, M, crossing rank, intersection dimension).
        cases = [([], [z[0]], 0, 0),
                 ([z[0], z[1]], [z[1], z[2]], 0, 1),
                 ([z[0], z[1], z[2]], [x[0], z[1], z[3]], 1, 1),
                 ([z[0], z[1], z[2]], [x[0], x[1], z[2]], 2, 1),
                 ([z[0], z[1], z[2]], [x[0], x[1], x[2]], 3, 0)]
        scramble = [('H', 0), ('S', 1), ('CX', 0, 2), ('H', 3),
                    ('CX', 2, 1), ('S', 0), ('CX', 3, 0), ('H', 1)]
        for original_left, original_right, expected_rank, intersection in cases:
            left, right = [[_clifford_vector(row, width, scramble) for row in rows]
                           for rows in (original_left, original_right)]
            for rows in (left, right):
                if len(rows) > 1:
                    rows[0] ^= rows[-1]
                if len(rows) > 2:
                    rows[-1] ^= rows[1]
                self.assertEqual(len(_binary_basis(rows)), len(rows))
                self.assertEqual(_pairing_rank(rows, rows, width), 0)
            self.assertEqual(_pairing_rank(left, right, width), expected_rank)
            self.assertEqual(len(left) + len(right) - len(_binary_basis(left + right)),
                             intersection)
            rank, z_axes, x_axes, coordinates = _simultaneous_coordinates(left, right, width)
            with self.subTest(rank=expected_rank, intersection=intersection):
                self.assertEqual(rank, expected_rank)
                self.assertEqual(len(_binary_basis(z_axes + x_axes)), 2 * width)
                for i in range(width):
                    for j in range(width):
                        self.assertEqual(_pairing(z_axes[i], z_axes[j], width), 0)
                        self.assertEqual(_pairing(x_axes[i], x_axes[j], width), 0)
                        self.assertEqual(_pairing(z_axes[i], x_axes[j], width), int(i == j))
                transformed_left = [coordinates(row) for row in left]
                transformed_right = [coordinates(row) for row in right]
                x_mask = (1 << width) - 1
                for row in transformed_left:
                    self.assertEqual(row & x_mask, 0)
                for row in transformed_right:
                    self.assertEqual((row & x_mask) >> rank, 0)
                    self.assertEqual((row >> width) & ((1 << rank) - 1), 0)
                self.assertEqual(_pairing_rank(transformed_left, transformed_right, width), rank)
                self.assertEqual(len(left) + len(right)
                                 - len(_binary_basis(transformed_left + transformed_right)),
                                 intersection)

    def test_active_native_family_attains_beta_equal_to_crossing_rank(self):
        for width in range(1, 4):
            layer = [('T', wire) for wire in range(width)]
            hadamards = [('H', wire) for wire in range(width)]
            # H S H maps Y to Z and preserves X, so this one-T-layer
            # word sends X^{tensor q} to the Hermitian Clifford H^{tensor q}.
            source_word = layer + [gate for wire in range(width)
                                   for gate in [('H', wire), ('S', wire), ('H', wire)]]
            source = _word_matrix(width, source_word)
            boundary = _word_matrix(width, hadamards)
            phase_layer = _word_matrix(width, layer)
            pauli = _pauli(width, {wire: X for wire in range(width)})
            np.testing.assert_allclose(source @ pauli @ source.conj().T,
                                       boundary, atol=ATOL, rtol=0)
            word = source_word + layer + hadamards + layer + _adjoint(source_word)
            unitary = _word_matrix(width, word)
            coefficient = np.trace(pauli @ unitary @ pauli @ unitary.conj().T) / (1 << width)
            expected = ((1 + 2 * np.sqrt(2)) / 4) ** width
            self.assertLess(abs(coefficient - expected), ATOL)
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 4 * width)

            for wire in range(width):
                np.testing.assert_allclose(boundary @ _pauli(width, {wire: Z}) @ boundary,
                                           _pauli(width, {wire: X}), atol=ATOL, rtol=0)
            rank = _pairing_rank([1 << (width + wire) for wire in range(width)],
                                 [1 << wire for wire in range(width)], width)
            self.assertEqual(rank, width)
            vector = boundary.reshape(-1) / np.sqrt(1 << width)
            kernel = np.kron(boundary, boundary.conj())
            diagonal = np.diag(np.kron(phase_layer, phase_layer.conj()))
            terms = (vector.conj()[:, None] * diagonal[:, None] * kernel
                     * diagonal[None, :] * vector[None, :])
            support_rank = int(np.log2(np.count_nonzero(np.abs(vector) > ATOL)))
            compatible_rank = int(np.log2(np.count_nonzero(np.abs(terms) > ATOL)))
            beta = Fraction(compatible_rank - support_rank - rank)
            self.assertEqual(support_rank, 2 * width)
            self.assertEqual(compatible_rank, 4 * width)
            self.assertEqual(beta, rank)
            phases = terms / (2 ** (-3 * width))
            np.testing.assert_allclose(np.abs(phases), 1, atol=ATOL, rtol=0)
            np.testing.assert_allclose(phases ** 8, 1, atol=ATOL, rtol=0)
            self.assertLess(abs(2 ** int(beta) * np.mean(phases) - coefficient), ATOL)

    def test_abstract_cubic_bias_and_disjoint_products_are_exact(self):
        for pair_count in range(1, 5):
            expected = (1 + Fraction(1, 1 << pair_count)) / 2
            self.assertEqual(_cubic_bias(pair_count, blocks=1), expected)
        for pair_count, blocks in ((1, 2), (1, 3), (2, 2), (3, 2)):
            bias = _cubic_bias(pair_count, blocks)
            expected = ((1 + Fraction(1, 1 << pair_count)) / 2) ** blocks
            self.assertEqual(bias, expected)
            scaled = Fraction(2) ** (blocks - 2) * bias
            self.assertEqual(scaled, Fraction(1, 4)
                             * (1 + Fraction(1, 1 << pair_count)) ** blocks)

    def test_growing_scaled_cubic_family_approaches_one_third(self):
        # Evaluate the factored exact formula, never a large truth table.
        # These are abstract polynomial biases, not native-circuit fixtures.
        with localcontext() as context:
            context.prec = 90
            target = Decimal(1) / 3
            logarithm = (Decimal(4) / 3).ln()
            previous_error = Decimal(1)
            previous_exponent = -1
            for pair_count in (4, 8, 12, 16, 24, 28):
                delta = Decimal(2) ** (-pair_count)
                blocks = int((logarithm / delta).to_integral_value(rounding=ROUND_FLOOR))
                exponent = blocks - 2
                value = (1 + delta) ** blocks / 4
                error = target - value
                # log(1+d) is between d-d^2/2 and d, and the floor
                # loses less than one block. This is the analytic O(d) bound.
                error_bound = (1 + logarithm / 2) * delta / 3
                self.assertGreater(exponent, previous_exponent)
                self.assertGreater(error, 0)
                self.assertLess(error, error_bound)
                self.assertLess(error, previous_error)
                previous_exponent, previous_error = exponent, error
            self.assertLess(previous_error, Decimal('2e-9'))


if __name__ == '__main__':
    unittest.main()
