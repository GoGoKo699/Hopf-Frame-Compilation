"""Full-matrix diagnostics for DISJOINT_BLOCK_COMPILER.md.

These tests check exact source/recovery algebra and the nontrivial QR
schedule. They do not certify asymptotic native T counts.
"""
from __future__ import annotations

import unittest

import numpy as np

ATOL = 4e-12


def _su(size, rng):
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    unitary, _ = np.linalg.qr(raw)
    unitary[:, 0] *= np.linalg.det(unitary).conjugate()
    return unitary


def _near_identity_su(size, rng):
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    hermitian = (raw + raw.conjugate().T) / 2
    hermitian -= np.trace(hermitian) * np.eye(size) / size
    values, vectors = np.linalg.eigh(hermitian)
    return (vectors * np.exp(0.03j * values)) @ vectors.conjugate().T


def _reflection_to(center, target):
    axis = -target.copy()
    axis[center] += 1
    axis /= np.linalg.norm(axis)
    return np.eye(len(target)) - 2 * np.outer(axis, axis.conjugate())


class DisjointBlockCompilerTests(unittest.TestCase):
    def test_local_two_diagonal_chart_retains_literal_complex_state(self):
        rng = np.random.default_rng(20261015)
        for size in (2, 4, 8):
            walsh = np.ones((1, 1))
            while len(walsh) < size:
                walsh = np.block([[walsh, walsh], [walsh, -walsh]]) / np.sqrt(2)
            inverse = np.linalg.inv(walsh @ np.diag(walsh[0]))
            center_phase = np.eye(size, dtype=complex)
            center_phase[0, 0] = 1j

            def chart(parameters):
                real_b = parameters[:size]
                real_c = np.r_[0.0, parameters[size:]]
                left = (walsh * np.exp(1j * (inverse @ real_b))) @ walsh.T
                right = (walsh * np.exp(1j * (inverse @ real_c))) @ walsh.T
                return left @ center_phase @ right @ center_phase.conjugate().T

            def coordinates(state):
                return np.r_[state.imag, state[1:].real]

            perturbation = rng.normal(size=size) + 1j * rng.normal(size=size)
            perturbation /= 128 * size * np.linalg.norm(perturbation)
            target = perturbation.copy()
            target[0] += 1
            target /= np.linalg.norm(target)
            parameters = np.zeros(2 * size - 1)
            for _ in range(24):
                parameters += coordinates(target) - coordinates(chart(parameters)[:, 0])
            word = chart(parameters)
            with self.subTest(size=size):
                # This comparison deliberately does not optimize global phase.
                np.testing.assert_allclose(word[:, 0], target, atol=ATOL, rtol=0)
                center_reflection = np.eye(size)
                center_reflection[0, 0] = -1
                np.testing.assert_allclose(
                    word @ center_reflection @ word.conjugate().T,
                    np.eye(size) - 2 * np.outer(target, target.conjugate()),
                    atol=ATOL, rtol=0,
                )

    def test_dirty_boolean_clock_echo_on_every_input_column(self):
        values = [0, 1, 1, 0]
        dimension = 4 * len(values)
        increment = np.zeros((dimension, dimension), dtype=complex)
        clock = np.zeros(dimension, dtype=complex)
        target = np.zeros(dimension, dtype=complex)
        for address, bit in enumerate(values):
            for dirty in range(4):
                column = 4 * address + dirty
                # Conditional carry precedes the low-bit toggle.
                low, high = dirty & 1, dirty >> 1
                high ^= bit & low
                low ^= bit
                increment[4 * address + low + 2 * high, column] = 1
                clock[column] = (1j) ** dirty
                target[column] = (1j) ** bit
        phase = np.diag(clock)
        echo = phase @ increment @ phase.conjugate().T @ increment.conjugate().T
        np.testing.assert_allclose(echo, np.diag(target), atol=ATOL, rtol=0)

    def test_adjacent_qr_mesh_preserves_padding_and_prior_zeros(self):
        rng = np.random.default_rng(20261013)
        for size in range(2, 9):
            padded = 1 << (size - 1).bit_length()
            matrix = np.eye(padded, dtype=complex)
            matrix[:size, :size] = _su(size, rng)
            schedule = {}
            for column in range(size - 1):
                for row in range(column + 1, size):
                    schedule.setdefault(padded - row + 2 * column, []).append((column, row))
            self.assertLessEqual(max(schedule), 2 * padded - 3)
            for time in sorted(schedule):
                pairs = schedule[time]
                touched = [index for _, row in pairs for index in (row - 1, row)]
                self.assertEqual(len(touched), len(set(touched)))
                self.assertEqual(len({row % 2 for _, row in pairs}), 1)
                for column, row in pairs:
                    if column:
                        np.testing.assert_allclose(
                            matrix[[row - 1, row], :column], 0, atol=ATOL, rtol=0
                        )
                    x, y = matrix[row - 1, column], matrix[row, column]
                    radius = np.sqrt(abs(x) ** 2 + abs(y) ** 2)
                    if radius == 0:
                        continue
                    gate = np.array([[x.conjugate(), y.conjugate()], [-y, x]]) / radius
                    self.assertAlmostEqual(np.linalg.det(gate), 1, delta=ATOL)
                    matrix[[row - 1, row], :] = gate @ matrix[[row - 1, row], :]
            with self.subTest(size=size):
                np.testing.assert_allclose(matrix, np.diag(np.diag(matrix)), atol=ATOL, rtol=0)
                np.testing.assert_allclose(
                    matrix[size:, size:], np.eye(padded - size), atol=ATOL, rtol=0
                )
                self.assertAlmostEqual(np.linalg.det(matrix), 1, delta=ATOL)

    def test_complex_pair_sources_and_complete_recovery(self):
        rng = np.random.default_rng(20261014)
        logical, modes = 32, 64
        active_blocks = ([1, 7], [3, 20])
        ideal = np.eye(logical, dtype=complex)
        accepted_target = np.eye(logical, dtype=complex)
        columns = np.eye(modes, dtype=complex)
        rows = np.eye(modes, dtype=complex)
        pair_phase = np.ones(modes, dtype=complex)
        permutation = np.arange(modes)
        free, eta = logical, 0.0
        active = []
        for labels in active_blocks:
            size = len(labels)
            target = _near_identity_su(size, rng)
            noise = 1e-7 * (
                rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
            )
            residual = target + noise
            eta = max(eta, np.linalg.norm(noise, 2))
            ideal[np.ix_(labels, labels)] = target
            accepted_target[np.ix_(labels, labels)] = residual / 2
            pairs = np.arange(free, free + size * size).reshape(size, size)
            free += size * size
            row_anchors = np.arange(free, free + size)
            free += size
            for column, anchor in enumerate(labels):
                state = np.zeros(modes)
                state[pairs[:, column]] = np.sqrt(np.abs(residual[:, column]) / 2)
                state[anchor] = np.sqrt(1 - np.sum(np.abs(residual[:, column])) / 2)
                columns = _reflection_to(anchor, state) @ columns
            for row, anchor in enumerate(row_anchors):
                state = np.zeros(modes)
                state[pairs[row, :]] = np.sqrt(np.abs(residual[row, :]) / 2)
                state[anchor] = np.sqrt(1 - np.sum(np.abs(residual[row, :])) / 2)
                rows = _reflection_to(anchor, state) @ rows
                permutation[labels[row]], permutation[anchor] = anchor, labels[row]
            pair_phase[pairs.ravel()] = np.exp(-1j * np.angle(residual.ravel()))
            active.extend(labels)
        self.assertLessEqual(free - logical, logical / 2)
        swap = np.eye(modes)[permutation]
        overlap = swap.conjugate().T @ rows.conjugate().T @ np.diag(
            pair_phase.conjugate()
        ) @ columns
        np.testing.assert_allclose(
            overlap.conjugate().T @ overlap, np.eye(modes), atol=ATOL, rtol=0
        )
        accepted = overlap[:logical, :logical]
        rejected = overlap[logical:, :logical]
        np.testing.assert_allclose(accepted, accepted_target, atol=ATOL, rtol=0)
        sign = np.diag(np.r_[np.ones(logical), -np.ones(logical)])
        recovery = -overlap @ sign @ overlap.conjugate().T @ sign @ overlap
        np.testing.assert_allclose(
            recovery[:logical, :logical],
            3 * accepted - 4 * accepted @ accepted.conjugate().T @ accepted,
            atol=ATOL, rtol=0,
        )
        np.testing.assert_allclose(
            recovery[logical:, :logical],
            rejected @ (np.eye(logical) - 4 * accepted.conjugate().T @ accepted),
            atol=ATOL, rtol=0,
        )
        inactive = sorted(set(range(logical)) - set(active))
        recovery[inactive, :] *= -1
        expected = np.vstack([ideal, np.zeros((logical, logical))])
        self.assertLess(np.linalg.norm(recovery[:, :logical] - expected, 2), 6 * eta)
        np.testing.assert_allclose(
            recovery[:, inactive], expected[:, inactive], atol=ATOL, rtol=0
        )


if __name__ == '__main__':
    unittest.main()
