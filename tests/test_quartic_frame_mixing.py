"""Quartic cube correlation, full-register boundaries, and a native dirty echo."""

from functools import reduce
from itertools import permutations, product
from operator import xor
from random import Random
import unittest

import numpy as np

from tests.test_clifford_orbit_query_rank import _rank
from tests.test_joint_lookup_pauli_source import _inverse
from tests.test_projective_field_recoupler import equal_ring, ring_matrix
from tests.test_recursive_field_shear_decoder import _tower_product, _tower_trace
from tests.test_returned_diagonal_echo_degree import _ccz


def _differences(table, base, directions):
    points = [base]
    for h in directions:
        points += [x ^ h for x in points]
    return sum((-1)**i.bit_count() * table[x] for i, x in enumerate(points))


class QuarticFrameMixingTests(unittest.TestCase):
    def test_quartic_polarization_and_weighted_cubic_derivatives(self):
        quartic = [int(x == 15) for x in range(16)]
        full_rank = 0
        for directions in product(range(16), repeat=4):
            polarization = _differences(quartic, 0, directions) % 2
            self.assertEqual(polarization, _rank(directions) == 4)
            full_rank += polarization
        self.assertEqual(full_rank, 20160)
        rng = Random(721634)
        for _ in range(20):
            coefs = [rng.randrange(8) * 2**max(0, mask.bit_count()-1) for mask in range(32)]
            table = [sum(c for mask, c in enumerate(coefs) if x & mask == mask) % 8
                     for x in range(32)]
            for _ in range(50):
                self.assertEqual(_differences(table, rng.randrange(32),
                                              rng.choices(range(32), k=4)) % 8, 0)

    def test_exact_native_dirty_quartic_all_32_columns(self):
        def flatten(schedule):
            return [gate for _, band in schedule for gate in band]
        phase = flatten(_ccz(4, 2, 3))
        toggle = [('H', 4)] + flatten(_ccz(0, 1, 4)) + [('H', 4)]
        native = phase + toggle + phase + _inverse(toggle)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in native), 28)
        expected = np.zeros((4, 32, 32), dtype=np.int64)
        expected[0] = np.diag([-1 if x & 15 == 15 else 1 for x in range(32)])
        self.assertTrue(equal_ring(ring_matrix(native, 5), (expected, 0)))

    def test_affine_support_intersection_and_path_sum_with_both_boundary_cliffords(self):
        # The left and right Clifford supports have different dimensions.
        h = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        ident = np.eye(2)
        cx = np.eye(4)[:, [0, 3, 2, 1]]
        left = np.kron(h, ident) @ cx
        right = cx @ np.kron(h, h)
        diag = np.diag(np.exp(1j*np.pi/4*np.array([0, 2, 3, 1])))
        middle = np.kron(ident, h)
        reconstructed = np.zeros(4, dtype=complex)
        for shift in (0, 1):
            for x in range(4):
                path = [left[x, v ^ shift] * right[v, x] * diag[v, v]
                        * (-1)**((v & (v ^ shift)) & 1) for v in range(4)]
                common = sum(abs(value) > 1e-10 for value in path)
                self.assertLessEqual(abs(sum(path)), common / np.sqrt(2*4) + 1e-12)
                reconstructed[x] += sum(path) / np.sqrt(2)
        np.testing.assert_allclose(reconstructed, np.diag(left @ middle @ diag @ right),
                                   atol=1e-12)

    def test_moore_adjoint_matches_explicit_third_derivative(self):
        def mul(a, b):
            return _tower_product(level, a, b)
        def power(a, degree):
            return reduce(mul, [a] * degree, 1)
        def trace(a):
            return _tower_trace(level, a)
        rng = Random(431)
        for level, size in [(3, 16), (4, 256)]:
            f = [trace(power(x, 15)) for x in range(size)]
            for _ in range(45):
                directions = rng.choices(range(size), k=3)
                minors = []
                for omitted in range(4):
                    orders = [j for j in range(4) if j != omitted]
                    minor = reduce(xor, (reduce(mul, (power(v, 1 << j)
                                    for v, j in zip(directions, ordering)), 1)
                                    for ordering in permutations(orders)), 0)
                    minors.append(minor)
                coeff = reduce(xor, (power(c, 1 << (3-j)) for j, c in enumerate(minors)))
                constant = _differences(f, 0, directions) % 2
                vanishes = True
                for x in range(size):
                    bit = (constant + _differences(f, x, directions)) % 2
                    expected = trace(reduce(xor, (mul(c, power(x, 1 << j))
                                            for j, c in enumerate(minors))))
                    self.assertEqual(bit, expected)
                    vanishes &= bit == 0
                self.assertEqual(vanishes, coeff == 0)


if __name__ == '__main__':
    unittest.main()
