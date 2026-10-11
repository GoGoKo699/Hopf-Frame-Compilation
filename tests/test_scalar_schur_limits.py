"""Independent complete-matrix and exact-completion checks for the Schur gauge.

These are finite algebraic fixtures, not a native source implementation or
a test of the unrestricted endpoint.
"""
from fractions import Fraction as F
import math
import unittest

import numpy as np

from tests.exact_matrices import CQ2, CO, CZ


def _word(n, coins):
    result = np.eye(1 << n, dtype=complex)
    for depth in range(n):
        half = 1 << (n - depth - 1)
        for prefix in range(1 << depth):
            a, b = coins[(1 << depth) - 1 + prefix]
            left, right = 2 * prefix * half, (2 * prefix + 1) * half
            x, y = result[left].copy(), result[right].copy()
            result[left] = a * x - b.conjugate() * y
            result[right] = b * x + a.conjugate() * y
    return result


def _gauge(n, angles, coins):
    size = 1 << n
    node = np.ones(2 * size, dtype=complex)
    columns = np.ones(size, dtype=complex)
    edges = np.zeros(size)
    for depth in range(n):
        half = 1 << (n - depth - 1)
        for prefix in range(1 << depth):
            v = (1 << depth) + prefix
            a, b = coins[v - 1]
            c, s = np.cos(angles[v - 1]), np.sin(angles[v - 1])
            lam = c / abs(a)
            node[2*v] = node[v] * lam * abs(a) / a
            node[2*v+1] = node[v] * lam * abs(b) / b
            j = (2 * prefix + 1) * half
            columns[j] = a * b / (abs(a) * abs(b) * node[v])
            edges[j] = s * abs(a) / (c * abs(b))
    gram = np.zeros((size, size))
    gram[0, 0] = 1
    for j in range(1, size):
        parent = j - (j & -j)
        gram[j] = edges[j] * gram[parent]
        gram[j, j] = np.sqrt(max(0., 1 - edges[j]**2))
    return node[size:], columns, gram @ gram.T


def _determinant(matrix):
    """Exact Gaussian elimination, independent of bidiagonal proof."""
    work = [row[:] for row in matrix]
    determinant = CO
    for column in range(len(work)):
        pivot = next((j for j in range(column, len(work))
                      if work[j][column] != CZ), None)
        if pivot is None:
            return CZ
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            determinant = -determinant
        diagonal = work[column][column]
        determinant *= diagonal
        for row in range(column + 1, len(work)):
            ratio = work[row][column] / diagonal
            for j in range(column + 1, len(work)):
                work[row][j] -= ratio * work[column][j]
            work[row][column] = CZ
    return determinant


class ScalarSchurLimitsTests(unittest.TestCase):
    def test_complex_coarse_gauge_all_columns_including_zero_target(self):
        rng = np.random.default_rng(161803)
        for n in range(1, 6):
            size = (1 << n) - 1
            beta = rng.uniform(.15, .75, size)
            theta = beta * rng.uniform(0, 1, size)
            theta[::4] = 0
            coins = list(zip(np.cos(beta) * np.exp(1j*rng.normal(size=size)),
                             np.sin(beta) * np.exp(1j*rng.normal(size=size))))
            target = _word(n, list(zip(np.cos(theta), np.sin(theta))))
            coarse = _word(n, coins)
            rows, columns, kernel = _gauge(n, theta, coins)
            actual = rows[:, None] * coarse * kernel * columns[None, :]
            np.testing.assert_allclose(actual, target, atol=2e-14, rtol=0)
            np.testing.assert_allclose(np.diag(kernel), 1, atol=2e-14)

    def test_all_free_complex_upper_entries_leave_exact_determinant(self):
        for m in range(2, 9):
            qs = [F(j+1, j+3) for j in range(m-1)]
            matrix = []
            for i in range(m):
                row = []
                for j in range(m):
                    if i >= j:
                        value = CQ2(math.prod(qs[j:i]))
                    elif j == i+1:
                        value = CQ2(qs[i])
                    else:
                        value = CQ2(F(2*i-3*j, i+j+1), F(i+1, j+2))
                    row.append(value)
                matrix.append(row)
            wanted = CQ2(math.prod(1-q*q for q in qs))
            self.assertEqual(_determinant(matrix), wanted)

    def test_chain_support_and_minor_weight_lower_bound(self):
        for n in range(2, 8):
            beta = 1/(10*n)
            theta = math.atan(math.tan(beta)/2)
            size = 1 << n
            coarse = _word(n, [(math.cos(beta), math.sin(beta))]*(size-1))
            target = _word(n, [(math.cos(theta), math.sin(theta))]*(size-1))
            rows, columns, _ = _gauge(n, np.full(size-1, theta),
                                     [(math.cos(beta), math.sin(beta))]*(size-1))
            weighted = rows[:, None] * coarse * columns[None, :]
            labels = [size-(1 << (n-j)) for j in range(n+1)]
            for j in range(1, n+1):
                support = np.flatnonzero(coarse[:, labels[j]])
                np.testing.assert_array_equal(
                    support, np.arange(labels[j-1],
                                       labels[j-1]+(1 << (n-j+1))))
            r = min(4, n)
            for i in range(r+1):
                for j in range(r+1):
                    if i >= j or j == i+1:
                        x, y = labels[i], labels[j]
                        self.assertGreaterEqual(abs(weighted[x, y]),
                                                1/(2*(10*n)**r))
                        self.assertAlmostEqual(target[x, y]/weighted[x, y],
                                               .5**abs(i-j), places=13)
            if n >= 3:
                self.assertEqual(coarse[labels[1], labels[3]], 0)

    def test_bounded_rank_approximation_obeys_robust_gap(self):
        for m in range(2, 8):
            kernel = np.fromfunction(lambda i, j: .5**np.abs(i-j), (m, m))
            u, singular, vh = np.linalg.svd(kernel)
            singular[-1] = 0
            approximation = (u * singular) @ vh
            pattern = np.fromfunction(lambda i, j: i >= j-1, (m, m))
            delta = np.max(np.abs(approximation-kernel)[pattern])
            bound = (.75**(m-1) /
                     (m**(1+m/2)*max(1, np.max(abs(approximation)))**(m-1)))
            self.assertGreaterEqual(delta, bound)


if __name__ == "__main__":
    unittest.main()
