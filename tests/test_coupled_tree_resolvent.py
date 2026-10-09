"""Bounded checks for graded tree inverses and their scattering dilation.

The operator and all physical rejected columns are retained. These checks
exercise the analytic estimates; they do not assign a native circuit cost.
"""
import math
import unittest

import numpy as np


ATOL = 3e-13


def _norm(matrix):
    return np.linalg.norm(matrix, 2)


def _resolvent(shift, height):
    result = np.eye(len(shift), dtype=complex)
    power = result.copy()
    for _ in range(height):
        power = power @ shift
        result += power
    return result


def _tree_pair(height, delta, seed):
    rng = np.random.default_rng(seed)
    size = (2 << height) - 1
    base = np.zeros((size, size), dtype=complex)
    defect = np.zeros_like(base)
    for vertex in range(1, 1 << height):
        depth = vertex.bit_length() - 1
        b = rng.normal(size=2) + 1j * rng.normal(size=2)
        e = rng.normal(size=2) + 1j * rng.normal(size=2)
        if vertex % 4 == 0:
            b[1] = 0
        if vertex % 3 == 0:
            e[0] = 0
        if vertex % 7 == 0:
            e[:] = 0
        b *= rng.uniform(.1, 1) / np.linalg.norm(b)
        if np.linalg.norm(e):
            e *= delta * 2 ** (depth - height) * rng.uniform(.1, 1) / np.linalg.norm(e)
        base[2 * vertex - 1:2 * vertex + 1, vertex - 1] = b
        defect[2 * vertex - 1:2 * vertex + 1, vertex - 1] = e
    return base, defect


def _scattering(base, defect):
    """Assemble the literal root-to-leaf local completion, including rejection."""
    target = base + defect
    size, sigma = len(base), 3
    rho, beta, gamma = {}, {}, {}

    def children(vertex):
        return [child for child in (2 * vertex, 2 * vertex + 1) if child <= size]

    for vertex in range(size, 0, -1):
        descendants = children(vertex)
        beta[vertex] = sigma ** 2 - 1 - sum(
            rho[c] * abs(defect[c - 1, vertex - 1]) ** 2 for c in descendants)
        gamma[vertex] = 1 + sum(
            rho[c] * target[c - 1, vertex - 1].conjugate() * defect[c - 1, vertex - 1]
            for c in descendants)
        rho[vertex] = 1 + sum(
            rho[c] * abs(target[c - 1, vertex - 1]) ** 2 for c in descendants
        ) + abs(gamma[vertex]) ** 2 / beta[vertex]
    unitary = np.eye(2 * size, dtype=complex)
    local_errors = []
    for vertex in range(1, size + 1):
        descendants = children(vertex)
        p = np.array([1] + [math.sqrt(rho[c]) * target[c - 1, vertex - 1]
                           for c in descendants]
                     + [-gamma[vertex].conjugate() / math.sqrt(beta[vertex])]) / math.sqrt(rho[vertex])
        q = np.array([1] + [math.sqrt(rho[c]) * defect[c - 1, vertex - 1]
                           for c in descendants] + [math.sqrt(beta[vertex])]) / sigma
        columns = [p, q]
        for basis in np.eye(len(p), dtype=complex):
            candidate = basis.copy()
            for col in columns:
                candidate -= np.vdot(col, candidate) * col
            if np.linalg.norm(candidate) > 1e-10:
                columns.append(candidate / np.linalg.norm(candidate))
            if len(columns) == len(p):
                break
        local = np.column_stack(columns)
        local_errors.append(_norm(local.conj().T @ local - np.eye(len(p))))
        if descendants:
            local = local[[0, 3, 1, 2]]
        ids = [size + vertex - 1, vertex - 1] + [size + c - 1 for c in descendants]
        unitary[ids] = local @ unitary[ids]
    return unitary[list(range(size, 2 * size)) + list(range(size))], beta, local_errors


class CoupledTreeResolventTests(unittest.TestCase):
    def test_graded_powers_tails_and_input_support(self):
        for height in range(1, 6):
            for seed, delta in ((173, .2), (937, .7), (281, .95)):
                with self.subTest(height=height, delta=delta):
                    base, defect = _tree_pair(height, delta, seed + height)
                    shift = _resolvent(base, height) @ defect
                    coupled = _resolvent(shift, height)
                    bound = 1 + delta + delta ** 2 / (3 * (1 - delta / 7))
                    self.assertLessEqual(_norm(coupled), bound + ATOL)
                    partial = np.eye(len(base), dtype=complex)
                    power = partial.copy()
                    for order in range(1, height + 1):
                        power = power @ shift
                        power_bound = delta ** order / math.prod(2 ** j - 1 for j in range(1, order + 1))
                        self.assertLessEqual(_norm(power), power_bound + ATOL)
                        cap = (1 << (height - order + 1)) - 1
                        np.testing.assert_array_equal(power[:, cap:], 0)
                        self.assertLessEqual(np.linalg.matrix_rank(power, tol=1e-18), cap)
                        tail = coupled - partial
                        tail_bound = power_bound / (1 - delta / (2 ** (order + 1) - 1))
                        self.assertLessEqual(_norm(tail), tail_bound + ATOL)
                        np.testing.assert_allclose(tail[:, cap:], 0, atol=ATOL, rtol=0)
                        partial += power
                    np.testing.assert_array_equal(power @ shift, 0)

    def test_rank_one_tail_keeps_independent_descendant_parameters(self):
        height, delta = 4, .8
        size = (2 << height) - 1
        base, defect = np.zeros((size, size)), np.zeros((size, size))
        angles = {}
        for vertex in range(1, 1 << height):
            depth = vertex.bit_length() - 1
            theta = delta * 2 ** (depth - height) * (.3 + .2 * (vertex % 5) / 4)
            angles[vertex] = theta
            base[2 * vertex - 1, vertex - 1] = 1
            defect[2 * vertex - 1, vertex - 1] = -2 * math.sin(theta / 2) ** 2
            defect[2 * vertex, vertex - 1] = math.sin(theta)
        target, internal = base + defect, (1 << height) - 1
        np.testing.assert_allclose(target[:, :internal].T @ target[:, :internal], np.eye(internal), atol=ATOL, rtol=0)
        shift = _resolvent(base, height) @ defect
        last = np.linalg.matrix_power(shift, height)
        np.testing.assert_allclose(last, np.linalg.matrix_power(defect, height), atol=1e-22, rtol=0)
        np.testing.assert_array_equal(last[:, 1:], 0)
        self.assertEqual(np.count_nonzero(last[internal:, 0]), 1 << height)
        derivatives = []
        for vertex in range(1 << (height - 1), 1 << height):
            ancestor, prefactor = vertex, 1
            while ancestor > 1:
                parent = ancestor // 2
                prefactor *= defect[ancestor - 1, parent - 1]
                ancestor = parent
            derivative = np.zeros(size)
            derivative[2 * vertex - 1] = -prefactor * math.sin(angles[vertex])
            derivative[2 * vertex] = prefactor * math.cos(angles[vertex])
            self.assertGreater(np.linalg.norm(derivative), 0)
            derivatives.append(derivative / np.linalg.norm(derivative))
        jacobian = np.column_stack(derivatives)
        np.testing.assert_allclose(jacobian.T @ jacobian, np.eye(1 << (height - 1)), atol=ATOL, rtol=0)

    def test_complex_thirty_mode_scattering_is_unitary_with_correct_block(self):
        base, defect = _tree_pair(3, .7, 319)
        unitary, beta, local_errors = _scattering(base, defect)
        self.assertEqual(unitary.shape, (30, 30))
        self.assertGreaterEqual(min(beta.values()), 9 - (43 / 18) ** 2)
        self.assertLess(max(local_errors), ATOL)
        expected = _resolvent(base + defect, 3) @ (np.eye(15) - base) / 3
        np.testing.assert_allclose(unitary[:15, :15], expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(unitary.conj().T @ unitary, np.eye(30), atol=ATOL, rtol=0)

    def test_zero_defect_rejected_block_contains_prefix_history(self):
        base, defect = _tree_pair(3, .7, 319)
        for column in range(7):
            base[:, column] /= np.linalg.norm(base[:, column])
        defect[:] = 0
        unitary, _, _ = _scattering(base, defect)
        fixed = np.kron(np.array([[1, -math.sqrt(8)], [math.sqrt(8), 1]]) / 3, np.eye(15))
        split = fixed.conj().T @ unitary
        path = np.zeros(15, complex)
        path[0] = 1
        for vertex in range(2, 16):
            path[vertex - 1] = base[vertex - 1, vertex // 2 - 1] * path[vertex // 2 - 1]
        np.testing.assert_allclose(split[:15, :15], np.eye(15), atol=ATOL, rtol=0)
        np.testing.assert_allclose(split[15:, :15], 0, atol=ATOL, rtol=0)
        np.testing.assert_allclose(split[15:, 15], -path / 2, atol=ATOL, rtol=0)
        for depth in range(4):
            self.assertAlmostEqual(np.linalg.norm(path[(1 << depth) - 1:(2 << depth) - 1]), 1)


if __name__ == '__main__':
    unittest.main()
