"""Finite fixtures for the four-band scalar gap and linear interface bounds."""
from itertools import combinations, product
from math import comb, cos, inf, pi, sin, sqrt
import unittest

import numpy as np

from tests.test_four_layer_crossing_rank import _binary_basis, _weighted_phase_fixture
from tests.test_one_band_walsh_blocks import _affine_spaces
from tests.test_operator_source_compiler import I2, X, Y, Z, _pauli


def _valuation(vector):
    """Exponent of 1-zeta in the integral basis 1,zeta,zeta²,zeta³."""
    if not any(vector):
        return inf
    order = 0
    while sum(vector) % 2 == 0:
        first = (vector[0] - sum(vector[1:])) // 2
        vector = [first, first + vector[1], first + sum(vector[1:3]),
                  first + sum(vector[1:])]
        order += 1
    return order


def _mobius_table(values):
    coefficients = values.copy()
    for bit in range((len(values) - 1).bit_length()):
        step = 1 << bit
        for offset in range(0, len(values), 2 * step):
            coefficients[offset + step:offset + 2 * step] -= coefficients[offset:offset + step]
    return coefficients


class NativeGapRefinementTests(unittest.TestCase):
    def test_affine_support_incidence_and_scalar_threshold(self):
        spaces = [set(points) for _, _, points in _affine_spaces(3)]
        for support in spaces:
            near = [other for other in spaces
                    if 4 * len(support & other) ** 2 > len(support) * len(other)]
            for other in near:
                self.assertTrue(support <= other or other <= support)
                self.assertIn(len(other) / len(support), (0.5, 1, 2))
            for first, second in combinations(near, 2):
                if len(first) == len(second) == 2 * len(support):
                    self.assertFalse((first - support) & (second - support))
                if 2 * len(first) == len(support) == 2 * len(second) and first | second != support:
                    self.assertGreaterEqual(4 * len(support - (first | second)), len(support))
        bound = 2 * sin(pi / 24)
        block_loss = (1 / sqrt(2) - 0.5)**2 + (1 / sqrt(2) - sqrt(3) / 2)**2
        self.assertAlmostEqual(bound**2, block_loss)
        self.assertAlmostEqual(bound**2 / 16, (1 - cos(pi / 12)) / 8)
        self.assertLess(8 / 128, bound**2)
        self.assertGreater(8 / 64, bound**2)

    def test_native_sparse_probe_has_two_orthogonal_clifford_terms(self):
        phase = np.diag([1, np.exp(1j*pi/4)])
        band = np.kron(phase, phase.conj())
        hadamard = (X + Z) / sqrt(2)
        cnot = np.eye(4)[[0, 1, 3, 2]]
        core = cnot @ np.kron(hadamard, Z) @ cnot
        image = band @ core @ band.conj().T
        axes = (I2, X, Z, Y)
        paulis = {}
        for x, z in product(range(4), repeat=2):
            factors = {w: axes[((x >> w) & 1) + 2 * ((z >> w) & 1)] for w in range(2)}
            paulis[x, z] = _pauli(2, factors)
        for shift in (0, 1, 2):
            for mask in range(4):
                probe = paulis[shift, mask]
                pulled = band.conj().T @ probe @ band
                coefficients = [np.trace(p @ pulled)/4 for p in paulis.values()]
                terms = [(c, band @ core @ p @ core @ band.conj().T)
                         for c, p in zip(coefficients, paulis.values()) if abs(c) > 1e-10]
                self.assertEqual(len(terms), 1 + shift.bit_count())
                np.testing.assert_allclose(sum(c * m for c, m in terms),
                                           image @ probe @ image, atol=2e-12)
                for _, matrix in terms:
                    for p in paulis.values():
                        output = matrix @ p @ matrix.conj().T
                        self.assertAlmostEqual(max(abs(np.trace(output @ q)) for q in paulis.values()), 4)
                if len(terms) == 2:
                    self.assertAlmostEqual(abs(np.trace(terms[0][1] @ terms[1][1])), 0)

    def test_exponent_digits_after_quadratic_carry(self):
        width = 7
        rng = np.random.default_rng(75103)
        points = range(1 << width)
        for _ in range(12):
            parities = rng.integers(1, 1 << width, size=12)
            quadratics = rng.choice([m for m in points if m.bit_count() <= 2], size=7)
            values = []
            for x in points:
                total = sum((int(mask) & x).bit_count() % 2 for mask in parities)
                for digit in range(3):
                    self.assertEqual((total >> digit) & 1, comb(total, 1 << digit) % 2)
                correction = sum((x & int(mask)) == mask for mask in quadratics) % 2
                values.append((total - 2 * correction) % 8)
            for digit in range(3):
                bits = (np.array(values) >> digit) & 1
                coefficients = _mobius_table(bits) % 2
                self.assertTrue(all(x.bit_count() <= 1 << digit
                                    for x in points if coefficients[x]))

    def test_cyclotomic_filtration_and_minimum_stratum(self):
        rng = np.random.default_rng(74025)
        roots = np.vstack((np.eye(4, dtype=int), -np.eye(4, dtype=int)))
        for width in range(2, 7):
            size = 1 << width
            marked = np.zeros((size, 4), dtype=int)
            pi_power = roots[0].copy()
            for _ in range(width):
                shifted = np.roll(pi_power, 1)
                shifted[0] *= -1
                pi_power = pi_power - shifted
            marked[-1] = pi_power
            self.assertEqual(_valuation(marked[-1]), width)
            for _ in range(6):
                samples = marked.copy()
                for _ in range(3):
                    phase = _weighted_phase_fixture(width, rng)
                    for x in range(size):
                        exp = sum(c for mask, c in phase.items() if (mask & x) == mask)
                        samples[x] += roots[exp % 8]
                coefficients = _mobius_table(samples)
                self.assertTrue(all(_valuation(v) >= mask.bit_count()
                                    for mask, v in enumerate(coefficients)))
                orders = [_valuation(v) for v in samples]
                lowest = min(orders)
                self.assertEqual(lowest, min(map(_valuation, coefficients)))
                self.assertGreaterEqual(orders.count(lowest) * (1 << lowest), size)

    def test_occupied_quadratic_support_and_outer_crossings(self):
        width = 4
        monomials = [mask for mask in range(1 << width) if mask.bit_count() <= 2]
        for coeffs in range(1, 1 << len(monomials)):
            terms = [m for i, m in enumerate(monomials) if (coeffs >> i) & 1]
            support = {x for x in range(1 << width)
                       if sum(x & m == m for m in terms) % 2}
            reached = {next(iter(support))}
            while True:
                more = {x for x in support if any((x ^ y).bit_count() <= 3 for y in reached)}
                if more == reached:
                    break
                reached = more
            self.assertEqual(reached, support)
        for rows in product(range(8), repeat=3):
            form = [sum((((x & y).bit_count() % 2) << j) for j, y in enumerate(rows))
                    for x in rows]
            self.assertLessEqual(len(_binary_basis(form)), len(_binary_basis(rows)))
        self.assertEqual(sum(6**i for i in range(4)), 259)
        self.assertEqual(sum(126**i for i in range(64)), (126**64 - 1) // 125)
