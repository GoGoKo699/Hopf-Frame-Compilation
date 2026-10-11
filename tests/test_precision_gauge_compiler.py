"""Independent all-column checks for the precision-gauge source identities.

Exact fixtures retain complex coarse phases and singular target corners.
Numerical completions exercise arbitrary dirty rejections, not just scalar
padding. Neither fixture assigns a native cost to a matrix factorization.
"""
from fractions import Fraction as F
import unittest

import numpy as np

from tests.exact_matrices import (
    CQ2, Q2, CO, CZ, CI, add, adjoint, eye, multiply, zeros,
)


INV_SQRT2 = CQ2(Q2(0, F(1, 2)))
PHASES = (CQ2(F(3, 5), F(4, 5)), CQ2(F(5, 13), -F(12, 13)),
          CQ2(-F(7, 25), F(24, 25)))
TARGETS = ((CQ2(F(3, 5)), CQ2(F(4, 5))),
           (CQ2(-F(5, 13)), CQ2(F(12, 13))),
           (CQ2(F(15, 17)), CQ2(-F(8, 17))))


def _node(n, vertex):
    depth = vertex.bit_length() - 1
    prefix = vertex - (1 << depth)
    left = prefix << (n - depth)
    right = left + (1 << (n - depth - 1))
    return depth, left, right


def _su(a, b):
    return [[a, -b.conj()], [b, a.conj()]]


def _apply_pair(matrix, first, second, local):
    a, b = matrix[first][:], matrix[second][:]
    matrix[first] = [local[0][0] * x + local[0][1] * y for x, y in zip(a, b)]
    matrix[second] = [local[1][0] * x + local[1][1] * y for x, y in zip(a, b)]


def _word(n, gates):
    result = eye(1 << n)
    for vertex, gate in enumerate(gates, 1):
        _, first, second = _node(n, vertex)
        _apply_pair(result, first, second, gate)
    return result


def _subtrees(n, gates):
    """Bottom-up vectors, independent of the chronological row-action word."""
    size = 1 << n
    unit = eye(size)
    root = {size + leaf: unit[leaf] for leaf in range(size)}
    injected = zeros(size, size - 1)
    for vertex in range(size - 1, 0, -1):
        a, b = gates[vertex - 1][0][0], gates[vertex - 1][1][0]
        left, right = root[2 * vertex], root[2 * vertex + 1]
        root[vertex] = [a * x + b * y for x, y in zip(left, right)]
        for row, (x, y) in enumerate(zip(left, right)):
            injected[row][vertex - 1] = INV_SQRT2 * (x - CI * y)
    internal = [[root[v][row] for v in range(1, size)] for row in range(size)]
    return internal, injected


def _diagonal(values):
    matrix = zeros(len(values))
    for j, value in enumerate(values):
        matrix[j][j] = value
    return matrix


def _unitary(rng, size):
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    unitary, triangular = np.linalg.qr(raw)
    return unitary @ np.diag(np.diag(triangular) / np.abs(np.diag(triangular)))


class PrecisionGaugeCompilerTests(unittest.TestCase):
    def test_regular_gauge_exact_all_columns_and_complex_phases(self):
        for n in range(1, 4):
            size = 1 << n
            target, coarse, d = [], [], {1: CO}
            ctable, atable = [CO] * size, zeros(size - 1, size)
            for vertex in range(1, size):
                c, s = TARGETS[vertex % len(TARGETS)]
                # The exact identity also covers corners with literal signs.
                if vertex % 5 == 0:
                    c, s = CZ, -CO
                target.append(_su(c, s))
                _, _, detail = _node(n, vertex)
                if c == CZ or s == CZ:
                    a, b, alpha, beta = c, s, CZ, CO
                    d[2 * vertex] = d[2 * vertex + 1] = d[vertex]
                else:
                    a = CQ2(F(3, 5)) * PHASES[vertex % 3]
                    b = CQ2(F(4, 5)) * PHASES[(vertex + 1) % 3]
                    d[2 * vertex] = d[vertex] * c / a
                    d[2 * vertex + 1] = d[vertex] * s / b
                    beta = c * s / (a * b)
                    alpha = -(c * c * b.conj() / a) + s * s * a.conj() / b
                coarse.append(_su(a, b))
                ctable[detail] = CO / (d[vertex] * beta)
                atable[vertex - 1][detail] = -alpha / (d[vertex] * beta)
            internal, _ = _subtrees(n, coarse)
            reconstructed = multiply(_diagonal([d[size + j] for j in range(size)]),
                                     add(multiply(_word(n, coarse), _diagonal(ctable)),
                                         multiply(internal, atable)))
            self.assertEqual(reconstructed, _word(n, target))

    def test_complex_parseval_source_exact_recurrence(self):
        for n in range(1, 4):
            size = 1 << n
            gates = [_su(CQ2(F(3, 5)) * PHASES[v % 3],
                         CQ2(F(4, 5)) * PHASES[(v + 1) % 3])
                     for v in range(1, size)]
            split = [(CQ2(F(12, 13)) * PHASES[d % 3],
                      CQ2(F(5, 13)) * PHASES[(d + 1) % 3]) for d in range(n)]
            physical = eye(size) + zeros(size, size)
            gamma, omega = [CO] * size, [CO] * (size - 1)
            deeper = CO
            for depth in range(n - 1, -1, -1):
                r, t = split[depth]
                for vertex in range(1 << depth, 2 << depth):
                    _, carrier, detail = _node(n, vertex)
                    _apply_pair(physical, carrier, detail, adjoint(gates[vertex - 1]))
                    _apply_pair(physical, carrier, size + detail, _su(r, t))
                    gamma[detail], omega[vertex - 1] = deeper, t * deeper
                deeper *= r
            gamma[0] = deeper
            internal, _ = _subtrees(n, gates)
            formula = multiply(_diagonal(gamma), adjoint(_word(n, gates)))
            features = multiply(_diagonal(omega), adjoint(internal))
            arranged = zeros(size, size)
            for vertex in range(1, size):
                arranged[_node(n, vertex)[2]] = features[vertex - 1]
            self.assertEqual(physical, formula + arranged)
            self.assertEqual(multiply(adjoint(physical), physical), eye(size))
            self.assertEqual(physical[size], [CZ] * size)

    def test_lossy_gauge_exact_full_frame_and_local_cancellation(self):
        for n in range(1, 4):
            size = 1 << n
            target, coarse, d = [], [], {1: CO}
            ctable, atable = [CO] * size, zeros(size - 1, size)
            for vertex in range(1, size):
                c, s = TARGETS[vertex % len(TARGETS)]
                z = F(3, 5) * PHASES[vertex % 3]
                t = (c + CI * s) * (c + CI * s)
                a, b = (CO + z) / 2, (z - CO) / (2 * CI)
                gate = [[a, -b], [b, a]]
                target.append(_su(c, s))
                coarse.append(gate)
                d[2 * vertex] = d[vertex] * c / a
                d[2 * vertex + 1] = d[vertex] * s / b
                cvalue = (z - (c * c - s * s)) / (2 * CI * d[vertex] * c * s)
                avalue = -(z - t) * (z - t.conj()) * INV_SQRT2 / (2 * d[vertex] * c * s)
                alternative = -(2 * INV_SQRT2 * CI * (z - t) / d[vertex])
                alternative -= (z - t) * (z - t) * INV_SQRT2 / (2 * d[vertex] * c * s)
                self.assertEqual(avalue, alternative)
                detail = _node(n, vertex)[2]
                ctable[detail], atable[vertex - 1][detail] = cvalue, avalue
                ell = [[INV_SQRT2], [-CI * INV_SQRT2]]
                self.assertEqual(multiply(gate, ell), [[z * row[0]] for row in ell])
            _, injections = _subtrees(n, coarse)
            reconstructed = multiply(_diagonal([d[size + j] for j in range(size)]),
                                     add(multiply(_word(n, coarse), _diagonal(ctable)),
                                         multiply(injections, atable)))
            self.assertEqual(reconstructed, _word(n, target))

    def test_uneven_source_keeps_initial_attenuation_failure(self):
        # Three physical leaves: the left root child terminates at depth 1,
        # while the right child has two leaves at depth 2.
        parent = _su(CQ2(F(3, 5)) * PHASES[0], CQ2(F(4, 5)) * PHASES[1])
        child = _su(CQ2(F(5, 13)) * PHASES[2], CQ2(F(12, 13)) * PHASES[0])
        split0 = _su(CQ2(F(12, 13)) * PHASES[1], CQ2(F(5, 13)) * PHASES[2])
        split1 = _su(CQ2(F(15, 17)) * PHASES[0], CQ2(F(8, 17)) * PHASES[1])
        r0, t0, r1, t1 = split0[0][0], split0[1][0], split1[0][0], split1[1][0]
        word = eye(3)
        _apply_pair(word, 0, 1, parent)
        _apply_pair(word, 1, 2, child)
        root = [word[row][0] for row in range(3)]
        subtree = [CZ, child[0][0], child[1][0]]
        internal = [[root[row], subtree[row]] for row in range(3)]
        physical = _diagonal([r1, CO, CO]) + zeros(3, 3)
        _apply_pair(physical, 1, 2, adjoint(child))
        _apply_pair(physical, 1, 5, split1)
        _apply_pair(physical, 0, 1, adjoint(parent))
        _apply_pair(physical, 0, 4, split0)
        formula = multiply(_diagonal([r0 * r1, r1, CO]), adjoint(word))
        features = multiply(_diagonal([t0 * r1, t1]), adjoint(internal))
        self.assertEqual(physical, formula + [[CZ] * 3] + features)
        gram = multiply(adjoint(physical), physical)
        self.assertEqual(gram, _diagonal([CQ2(r1.abs2()), CO, CO]))
        self.assertNotEqual(gram, eye(3))

    def test_two_query_replay_preserves_the_actual_dirty_rejection(self):
        rng = np.random.default_rng(18373)
        for dirty_size in (2, 3, 5):
            identity = np.eye(dirty_size)
            for angle in (.3, -1.7, 2.2):
                z = .63 * np.exp(1j * angle)
                loss = 1 - abs(z) ** 2
                rejection = np.sqrt(loss) * _unitary(rng, dirty_size)
                other = _unitary(rng, dirty_size)
                word = np.block([[z * identity, -np.sqrt(loss) * other],
                                 [rejection, z.conjugate() * rejection @ other / np.sqrt(loss)]])
                np.testing.assert_allclose(word.conj().T @ word, np.eye(2 * dirty_size), atol=2e-14)
                phase = np.exp(1j * .83)
                signal = np.diag([phase] * dirty_size + [1] * dirty_size)
                actual = (word @ signal @ word.conj().T)[:, :dirty_size]
                wanted = np.vstack(((1 + (phase - 1) * abs(z) ** 2) * identity,
                                    (phase - 1) * z.conjugate() * rejection))
                np.testing.assert_allclose(actual, wanted, atol=2e-14)
                reversed_queries = (word.conj().T @ signal @ word)[:, :dirty_size]
                self.assertGreater(np.linalg.norm(reversed_queries - wanted), .1)

    def test_lossy_tree_and_four_bank_overlap_on_all_dirty_columns(self):
        rng = np.random.default_rng(9173)
        dirty_size = 3
        identity = np.eye(dirty_size)
        ell = np.array([1, -1j]) / np.sqrt(2)
        for n in (2, 3):
            size = 1 << n
            width = size * dirty_size
            target = np.eye(size, dtype=complex)
            coarse = target.copy()
            gates, rejected, losses = {}, {}, {}
            gauges = {1: 1 + 0j}
            coefficients, corrections = np.ones(size, dtype=complex), np.zeros(size, dtype=complex)
            for vertex in range(1, size):
                theta = (.13 + .11 * vertex) * (-1 if vertex % 3 == 0 else 1)
                c, s = np.cos(theta), np.sin(theta)
                tau = 1e-4 * min(abs(c), abs(s))
                z = (1 - tau) * np.exp(2j * theta)
                a, b = (1 + z) / 2, (z - 1) / (2j)
                gates[vertex] = np.array([[a, -b], [b, a]])
                losses[vertex] = 1 - abs(z) ** 2
                rejected[vertex] = np.sqrt(losses[vertex]) * _unitary(rng, dirty_size)
                _, first, detail = _node(n, vertex)
                target[[first, detail]] = np.array([[c, -s], [s, c]]) @ target[[first, detail]]
                coarse[[first, detail]] = gates[vertex] @ coarse[[first, detail]]
                gauges[2 * vertex] = gauges[vertex] * c / a
                gauges[2 * vertex + 1] = gauges[vertex] * s / b
                coefficients[detail] = (z - np.cos(2 * theta)) / (2j * gauges[vertex] * c * s)
                t = np.exp(2j * theta)
                corrections[detail] = -(z - t) * (z - t.conjugate()) / (2 * np.sqrt(2) * gauges[vertex] * c * s)

            # Apply fresh-feature three-port merges directly to all physical
            # logical/dirty input columns; no scalar replacement of F occurs.
            physical = np.vstack((np.eye(width), np.zeros((width, width)))).astype(complex)
            roots = {size + leaf: np.eye(size, dtype=complex)[:, leaf] for leaf in range(size)}
            injections = {}
            for vertex in range(size - 1, 0, -1):
                _, first, detail = _node(n, vertex)
                ids = list(range(first * dirty_size, (first + 1) * dirty_size))
                ids += list(range(detail * dirty_size, (detail + 1) * dirty_size))
                children = physical[ids].copy()
                physical[ids] = np.kron(gates[vertex].conj().T, identity) @ children
                slot = slice((size + detail) * dirty_size, (size + detail + 1) * dirty_size)
                np.testing.assert_array_equal(physical[slot], 0)
                physical[slot] = rejected[vertex] @ np.kron(ell.conj()[None, :], identity) @ children
                left, right = roots[2 * vertex], roots[2 * vertex + 1]
                roots[vertex] = gates[vertex][0, 0] * left + gates[vertex][1, 0] * right
                injections[vertex] = ell[0] * left + ell[1] * right
            np.testing.assert_allclose(physical.conj().T @ physical, np.eye(width), atol=3e-13)
            np.testing.assert_allclose(physical[:width], np.kron(coarse.conj().T, identity), atol=3e-13)
            for vertex in range(1, size):
                detail = _node(n, vertex)[2]
                slot = slice((size + detail) * dirty_size, (size + detail + 1) * dirty_size)
                expected = rejected[vertex] @ np.kron(injections[vertex].conj()[None, :], identity)
                np.testing.assert_allclose(physical[slot], expected, atol=3e-13)

            leaves = np.array([gauges[size + j] for j in range(size)])
            qsource, psource = np.zeros((4 * width, width), complex), np.zeros((4 * width, width), complex)
            qsource[:2 * width] = physical @ np.kron(np.diag(leaves.conj() / 1.25), identity)
            qsource[3 * width:] = np.kron(np.diag(np.sqrt(1 - abs(leaves / 1.25) ** 2)), identity)
            psource[:width] = np.kron(np.diag(coefficients / 1.6), identity)
            radial = abs(coefficients / 1.6) ** 2
            for vertex in range(1, size):
                detail = _node(n, vertex)[2]
                slot = slice((size + detail) * dirty_size, (size + detail + 1) * dirty_size)
                column = slice(detail * dirty_size, (detail + 1) * dirty_size)
                psource[slot, column] = corrections[detail] * rejected[vertex] / (1.6 * losses[vertex])
                radial[detail] += abs(corrections[detail]) ** 2 / (1.6 ** 2 * losses[vertex])
            self.assertLess(np.max(radial), 1)
            psource[2 * width:3 * width] = np.kron(np.diag(np.sqrt(1 - radial)), identity)
            np.testing.assert_allclose(psource.conj().T @ psource, np.eye(width), atol=3e-13)
            np.testing.assert_allclose(qsource.conj().T @ qsource, np.eye(width), atol=3e-13)
            np.testing.assert_allclose(qsource.conj().T @ psource, np.kron(target / 2, identity), atol=3e-13)

    def test_recovery_word_full_rejected_output_and_linear_error(self):
        rng = np.random.default_rng(2371)
        size = 4
        target = _unitary(rng, size)
        left, right = _unitary(rng, size), _unitary(rng, size)
        reflector = np.diag([-1] * size + [1] * size)
        for delta in (0, 1e-5, 1e-3):
            singular = .5 + delta * np.linspace(-1, 1, size)
            a = target @ right @ np.diag(singular) @ right.conj().T
            top_failure = target @ right @ np.diag(np.sqrt(1 - singular ** 2)) @ right.conj().T
            bottom_failure = left @ right @ np.diag(np.sqrt(1 - singular ** 2)) @ right.conj().T
            lower = left @ right @ np.diag(singular) @ right.conj().T
            unitary = np.block([[a, -top_failure], [bottom_failure, lower]])
            np.testing.assert_allclose(unitary.conj().T @ unitary, np.eye(2 * size), atol=2e-14)
            recovered = (-unitary @ reflector @ unitary.conj().T @ reflector @ unitary)[:, :size]
            formula = np.vstack((3 * a - 4 * a @ a.conj().T @ a,
                                 bottom_failure @ (np.eye(size) - 4 * a.conj().T @ a)))
            np.testing.assert_allclose(recovered, formula, atol=2e-14)
            error = np.linalg.norm(recovered - np.vstack((target, np.zeros((size, size)))), 2)
            self.assertLessEqual(error, 12 * delta + 2e-14)


if __name__ == '__main__':
    unittest.main()
