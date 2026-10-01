"""Small ideal-algebra checks of state-only reference-interference QBP.

The largest preparation matrix has dimension 64. Coins and coarse words
are ideal unitary fixtures: these checks do not emit a native compiler,
verify its asymptotic cost, or approximate the prescribed complete frame.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.frames import direct_real_frame, real_tree_data


ATOL = 4e-12
H = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2.)


def _fixtures(height):
    size = 1 << height
    return (
        np.linspace(.17, .93, size - 1),
        np.linspace(-.71, 1.13, size - 1),
        np.zeros(size - 1),
        np.full(size - 1, np.pi / 2),
        np.array([0 if node % 2 else np.pi / 2
                  for node in range(1, size)]),
        np.full(size - 1, np.pi / 4),
    )


def _envelope(angles):
    """One visit per internal node; M tracks the largest ancestor derivative."""
    size = len(angles) + 1
    incoming, maximum = np.zeros(2 * size), np.zeros(2 * size)
    incoming[1] = 1.
    for node, theta in enumerate(angles, start=1):
        c2, s2 = np.cos(theta) ** 2, np.sin(theta) ** 2
        a, m = incoming[node], maximum[node]
        incoming[2 * node:2 * node + 2] = a * np.array([c2, s2])
        maximum[2 * node] = max(m * c2, a * s2)
        maximum[2 * node + 1] = max(m * s2, a * c2)
    weights = maximum[size:]
    total = weights.sum()
    return weights, total, np.sqrt((weights + 1 / size) / (total + 1))


def _observable(size):
    axis = np.linspace(.2, .9, size) + 1j * np.linspace(.7, -.3, size)
    axis /= np.linalg.norm(axis)
    return np.eye(size) - 2 * np.outer(axis, axis.conj())


def _positive_completion(state):
    """An ideal Hopf frame with the supplied positive reference state column."""
    size = len(state)
    mass = np.zeros(2 * size)
    mass[size:] = state ** 2
    angles = np.zeros(size - 1)
    for node in range(size - 1, 0, -1):
        mass[node] = mass[2 * node] + mass[2 * node + 1]
        angles[node - 1] = np.arctan2(np.sqrt(mass[2 * node + 1]),
                                      np.sqrt(mass[2 * node]))
    return direct_real_frame(size.bit_length() - 1, angles)


def _coarse_fixture(frame, scale=1.):
    """A complex coarse completion, close on the state column including phase."""
    size = len(frame)
    error = np.eye(size, dtype=complex)
    error[0, 0] = np.exp(.013j * scale / np.sqrt(size))
    for label, phase, angle in ((1, .41, .027), (size - 1, -.63, -.019)):
        theta = angle * scale / np.sqrt(size)
        c, s = np.cos(theta), np.sin(theta)
        rotation = np.eye(size, dtype=complex)
        rotation[np.ix_([0, label], [0, label])] = [
            [c, -np.exp(1j * phase) * s], [np.exp(-1j * phase) * s, c]]
        error = rotation @ error
    return frame @ error


def _state_lcu(coarse, target):
    """Q in order (s,t,x); its accepted 00 block prepares phi/2."""
    size = len(target)
    phi = coarse.conj().T @ target
    walsh = np.array([[1.]])
    for _ in range(size.bit_length() - 1):
        walsh = np.kron(walsh, H)
    h_signal = np.kron(H, np.eye(2 * size))
    controlled_walsh = np.eye(4 * size, dtype=complex)
    controlled_walsh[2 * size:, 2 * size:] = np.kron(np.eye(2), walsh)
    coins = np.zeros((4 * size, 4 * size), dtype=complex)
    for signal in (0, 1):
        for x in range(size):
            z = phi[0] if signal == 0 else (0 if x == 0 else np.sqrt(size) * phi[x])
            assert abs(z) <= 1 + ATOL
            off = np.sqrt(max(0., 1 - abs(z) ** 2))
            indices = [2 * signal * size + x, (2 * signal + 1) * size + x]
            coins[np.ix_(indices, indices)] = [[z, -off], [off, np.conj(z)]]
    return h_signal @ coins @ controlled_walsh @ h_signal, phi


def _amplify(q, size, branches=1):
    initial = np.ones(4 * size)
    initial[0] = -1
    accepted = np.ones(4 * size)
    accepted[:size] = -1
    r_initial = np.kron(np.eye(branches), np.diag(initial))
    r_accepted = np.kron(np.eye(branches), np.diag(accepted))
    # The actual adjoint and leading literal minus are both essential.
    return -q @ r_initial @ q.conj().T @ r_accepted @ q


class ReferenceStateQBPTests(unittest.TestCase):
    def test_linear_envelope_matches_jacobian_and_extremal_totals(self):
        for height in (2, 3):
            size = 1 << height
            for angles in _fixtures(height):
                data = real_tree_data(angles)
                derivatives = np.asarray(data.derivatives)
                weights, total, reference = _envelope(angles)
                np.testing.assert_allclose(weights, np.max(derivatives ** 2, axis=0),
                                           atol=ATOL, rtol=0)
                self.assertGreaterEqual(total + ATOL, 1)
                self.assertLessEqual(total, height + ATOL)
                self.assertAlmostEqual(np.linalg.norm(reference), 1., delta=ATOL)
                self.assertTrue(np.all(reference > 0))
                scores = 2 * derivatives / reference
                for depth in range(height):
                    rows = scores[(1 << depth) - 1:(1 << (depth + 1)) - 1]
                    self.assertLessEqual(np.linalg.norm(rows, axis=0).max(),
                                         2 * np.sqrt(total + 1) + ATOL)
            self.assertAlmostEqual(_envelope(np.full(size - 1, np.pi / 4))[1],
                                   1., delta=ATOL)
            self.assertAlmostEqual(_envelope(np.zeros(size - 1))[1],
                                   height, delta=ATOL)

    def test_reference_measurement_decodes_every_raw_coordinate(self):
        for height in (2, 3):
            size = 1 << height
            observable = _observable(size)
            np.testing.assert_allclose(observable @ observable, np.eye(size), atol=ATOL)
            self.assertGreater(np.linalg.norm(observable.imag), .1)
            for angles in _fixtures(height):
                data = real_tree_data(angles)
                derivatives = np.asarray(data.derivatives)
                reference = _envelope(angles)[2]
                prepared = np.concatenate((reference, data.state)) / np.sqrt(2)
                controlled = np.eye(2 * size, dtype=complex)
                controlled[size:, size:] = observable
                output = np.kron(H, np.eye(size)) @ controlled @ prepared
                probabilities = (abs(output) ** 2).reshape(2, size)
                self.assertAlmostEqual(probabilities.sum(), 1., delta=ATOL)
                ratio = derivatives / reference
                decoded = 2 * ratio @ (probabilities[0] - probabilities[1])
                analytic = 2 * np.real(derivatives.conj() @ observable @ data.state)
                np.testing.assert_allclose(decoded, analytic, atol=ATOL, rtol=0)
                # No divisions by incoming amplitudes: vanishing raw
                # derivatives and negative oriented amplitudes are retained.
                tau = .003
                rounded = tau * np.round(ratio / tau)
                rounded_gradient = 2 * rounded @ (probabilities[0] - probabilities[1])
                self.assertLessEqual(abs(rounded_gradient - decoded).max(), 2 * tau + ATOL)

            # A signed state can itself be the reference when each local
            # branch probability is at least 1/4; the score bound is 2*sqrt(3).
            angles = np.linspace(-.9, -.6, size - 1)
            data = real_tree_data(angles)
            self.assertTrue(np.any(data.state < 0))
            derivatives = np.asarray(data.derivatives)
            ratio = derivatives / data.state
            self.assertLessEqual(abs(ratio).max(), np.sqrt(3) + ATOL)
            joint = np.concatenate((data.state, observable @ data.state)) / np.sqrt(2)
            probabilities = abs(np.kron(H, np.eye(size)) @ joint).reshape(2, size) ** 2
            decoded = 2 * ratio @ (probabilities[0] - probabilities[1])
            analytic = 2 * np.real(derivatives.conj() @ observable @ data.state)
            np.testing.assert_allclose(decoded, analytic, atol=ATOL, rtol=0)

    def test_complex_residual_lcu_and_one_step_amplification_literal_phase(self):
        for height in (2, 3):
            size = 1 << height
            frame = direct_real_frame(height, np.linspace(-.37, .81, size - 1))
            target = frame[:, 0]
            for scale in (0., 1.):
                coarse = _coarse_fixture(frame, scale)
                self.assertLessEqual(np.linalg.norm(coarse[:, 0] - target), 1 / (4 * np.sqrt(size)))
                q, phi = _state_lcu(coarse, target)
                np.testing.assert_allclose(q.conj().T @ q, np.eye(4 * size), atol=ATOL)
                np.testing.assert_allclose(q[:size, 0], phi / 2, atol=ATOL, rtol=0)
                amplified = _amplify(q, size)
                expected = np.zeros(4 * size, dtype=complex)
                expected[:size] = phi
                np.testing.assert_allclose(amplified[:, 0], expected, atol=ATOL, rtol=0)
                self.assertAlmostEqual(np.linalg.norm(-amplified[:, 0] - expected), 2., delta=ATOL)
                if scale:
                    self.assertGreater(np.linalg.norm(phi.imag), 1e-3)
                # Perturb a whole Q, use its literal adjoint, and test the
                # three-occurrence telescoping bound, including flag leakage.
                phase = np.ones(4 * size, dtype=complex)
                phase[size] = np.exp(.021j)
                changed = phase[:, None] * q
                error = np.linalg.norm(changed - q, ord=2)
                self.assertLessEqual(np.linalg.norm(_amplify(changed, size) - amplified, ord=2),
                                     3 * error + ATOL)

    def test_coherent_two_branch_preparation_returns_dirty_spectator(self):
        for height in (2, 3):
            size = 1 << height
            angles = np.linspace(-.49, .87, size - 1)
            frame = direct_real_frame(height, angles)
            reference = _envelope(angles)[2]
            reference_frame = _positive_completion(reference)
            targets = (reference, frame[:, 0])
            coarse = (_coarse_fixture(reference_frame), _coarse_fixture(frame, -.8))
            q = np.zeros((8 * size, 8 * size), dtype=complex)
            restore = np.zeros_like(q)
            initial = np.zeros((8 * size, 2), dtype=complex)
            expected = np.zeros_like(initial)
            for branch in (0, 1):
                start = branch * 4 * size
                q[start:start + 4 * size, start:start + 4 * size] = _state_lcu(coarse[branch], targets[branch])[0]
                restore[start:start + 4 * size, start:start + 4 * size] = np.kron(np.eye(4), coarse[branch])
                initial[start, branch] = 1
                expected[start:start + size, branch] = targets[branch]
            prepared = restore @ _amplify(q, size, branches=2)
            np.testing.assert_allclose(prepared @ initial, expected, atol=ATOL, rtol=0)
            # Equality on both branch columns and both dirty inputs includes
            # their coherent superpositions and arbitrary reference extensions.
            np.testing.assert_allclose(np.kron(prepared, np.eye(2)) @ np.kron(initial, np.eye(2)),
                                       np.kron(expected, np.eye(2)), atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
