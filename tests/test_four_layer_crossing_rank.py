"""Finite algebra audits for the four-layer crossing-rank obstruction.

These check the symplectic reduction, literal native boundary matrices,
weighted phase digits, initialized-observable normalization, and the
quantitative cluster constants. They neither certify an infinite-width
theorem numerically nor realize the abstract bias as a native coefficient.
"""
from __future__ import annotations

from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction
from itertools import combinations
import unittest

import numpy as np

from tests.test_clean_aware_depth_obstruction import (
    _degree, _f2_value, _integer_mobius, _phase_digits, _phase_events,
)
from tests.test_operator_source_compiler import ATOL, X, Y, Z, _adjoint, _pauli, _word_matrix


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


def _affine_phase_table(coefficients, width, rng):
    """Evaluate an actual invertible affine pullback, before interpolation."""
    rows = [1 << bit for bit in range(width)]
    for _ in range(4 * width):
        first, second = rng.choice(width, 2, replace=False)
        rows[first] ^= rows[second]
    offset = int(rng.integers(1 << width))
    values = []
    for assignment in range(1 << width):
        original = offset ^ sum(((assignment & row).bit_count() % 2) << bit
                                for bit, row in enumerate(rows))
        values.append(sum(value for mask, value in coefficients.items()
                          if original & mask == mask) % 8)
    return values


def _weighted_phase_fixture(width, rng):
    coefficients = {0: int(rng.integers(8))}
    for degree in (1, 2, 3):
        for support in combinations(range(width), degree):
            mask = sum(1 << bit for bit in support)
            coefficients[mask] = (int(rng.integers(8)) << (degree - 1)) % 8
    return coefficients


class FourLayerCrossingRankTests(unittest.TestCase):
    def test_weighted_phase_digits_and_residues_after_affine_pullback(self):
        width = 8
        rng = np.random.default_rng(740751)
        # The linear weight sum attains all three digit bounds and degree
        # seven for a residue event; affine changes preserve those degrees.
        fixtures = [{1 << bit: 1 for bit in range(width)}]
        fixtures += [_weighted_phase_fixture(width, rng) for _ in range(5)]
        for sample, coefficients in enumerate(fixtures):
            values = _affine_phase_table(coefficients, width, rng)
            pulled = _integer_mobius(values, width)
            for mask, coefficient in pulled.items():
                degree = mask.bit_count()
                self.assertLessEqual(degree, 3)
                if degree:
                    self.assertEqual(coefficient % (1 << (degree - 1)), 0)
            digits = _phase_digits(pulled)
            events = _phase_events(digits)
            self.assertTrue(all(_degree(digit) <= bound
                                for digit, bound in zip(digits, (1, 2, 4))))
            self.assertLessEqual(max(map(_degree, events)), 7)
            if sample == 0:
                self.assertEqual(list(map(_degree, digits)), [1, 2, 4])
                self.assertEqual(max(map(_degree, events)), 7)
            for assignment, value in enumerate(values):
                self.assertEqual([_f2_value(digit, assignment) for digit in digits],
                                 [(value >> bit) & 1 for bit in range(3)])
                self.assertEqual([_f2_value(event, assignment) for event in events],
                                 [int(value == residue) for residue in range(8)])

    def test_fiber_phase_coefficient_events_have_degrees_seven_three_one(self):
        # Three untouched variables t and three active fiber variables w.
        # Integer interpolation follows an affine substitution mixing both.
        rng = np.random.default_rng(740752)
        attained = [0, 0, 0, 0]
        for _ in range(10):
            values = _affine_phase_table(_weighted_phase_fixture(6, rng), 6, rng)
            pulled = _integer_mobius(values, 6)
            for w_mask in range(8):
                w_degree = w_mask.bit_count()
                coefficient = {mask & 7: value for mask, value in pulled.items()
                               if mask >> 3 == w_mask}
                for t_mask, value in coefficient.items():
                    total_degree = t_mask.bit_count() + w_degree
                    self.assertLessEqual(total_degree, 3)
                    if total_degree:
                        self.assertEqual(value % (1 << (total_degree - 1)), 0)
                events = _phase_events(_phase_digits(coefficient))
                event_degree = max(map(_degree, events))
                self.assertLessEqual(event_degree, (7, 3, 1, 0)[w_degree])
                attained[w_degree] = max(attained[w_degree], event_degree)
                for t in range(8):
                    value = sum(c for mask, c in coefficient.items() if t & mask == mask) % 8
                    self.assertEqual([_f2_value(event, t) for event in events],
                                     [int(value == residue) for residue in range(8)])
        self.assertEqual(attained[1:], [3, 1, 0])

    def test_initialized_two_sided_observable_identity_includes_native_leakage(self):
        # Logical wire 0, dirty wire 1, and clean wire 2. The four physical
        # T layers form a commutator of rotations about entangling axes.
        first_basis = [('CX', 1, 2), ('H', 2)]
        second_basis = [('CX', 0, 2)]
        first = _adjoint(first_basis) + [('T', 2)] + first_basis
        second = _adjoint(second_basis) + [('T', 2)] + second_basis
        word = first + second + _adjoint(first) + _adjoint(second)
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 4)
        unitary = _word_matrix(3, word)
        embedding = np.eye(8)[:, :4]
        projector = embedding @ embedding.conj().T
        isometry = unitary @ embedding
        eta = np.linalg.norm(isometry - embedding, ord=2)
        self.assertGreater(np.linalg.norm((np.eye(8) - projector) @ isometry, ord=2), 0.01)
        self.assertLess(eta, 0.5)
        clean_z = _pauli(3, {2: Z})
        for factors in ({}, {0: Y}, {1: X}):
            pauli = _pauli(3, factors)
            observable = np.sqrt(2) * projector @ pauli
            evolved = unitary @ observable @ unitary.conj().T
            coefficient = np.trace(observable @ evolved) / 8
            discrepancy = np.linalg.norm(evolved - observable, ord='fro') ** 2 / 8
            self.assertLess(abs(coefficient.imag), ATOL)
            self.assertAlmostEqual(np.trace(observable @ observable).real / 8, 1)
            self.assertAlmostEqual(1 - coefficient.real, discrepancy / 2, places=12)
            expanded = sum(np.trace(left @ unitary @ right @ unitary.conj().T) / 8
                           for left in (pauli, pauli @ clean_z)
                           for right in (pauli, pauli @ clean_z)) / 2
            self.assertLess(abs(coefficient - expanded), ATOL)
            self.assertLessEqual(1 - coefficient.real, 2 * eta ** 2 + ATOL)
            if not factors:
                # The one-sided identity witness is always one; this
                # two-sided projector witness instead detects clean leakage.
                self.assertAlmostEqual(np.trace(isometry.conj().T @ isometry).real / 4, 1)
                self.assertLess(coefficient.real, 0.99)

    def test_one_clean_scale_recurrence_pays_for_cross_cluster_interference(self):
        clean, count = 1, 4
        multiplier = 2 * (count ** 2 - 1)
        for rank in (0, 1, 3, 8):
            active = 2 * rank
            degree = 7 + 3 * active + active * (active - 1) // 2
            constant = count * degree + clean + 4 * active + 4 * 2 + 1
            threshold = 0
            for _ in range(count):
                following = multiplier * threshold + constant
                hull_codimension = 2 * (count - 1) * threshold
                log_gap = (-(count + 1) * hull_codimension - count * degree
                           - clean - 4 * active - 2 * 2)
                log_cross_bound = 2 * 2 - following
                self.assertEqual(log_cross_bound, log_gap - 1)
                self.assertGreaterEqual(log_gap - 1, -following)
                threshold = following
            self.assertEqual(threshold, 27931 * (8 * rank ** 2 + 28 * rank + 38))
            self.assertGreaterEqual(threshold, 27931 * 8 * rank ** 2)
            self.assertLessEqual(threshold, 27931 * 38 * (rank + 1) ** 2)

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
