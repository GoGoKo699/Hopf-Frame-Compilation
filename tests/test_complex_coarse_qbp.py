"""Small complex-chart checks with literal local native coarse blocks.

Four/eight-mode matrices test the phase gauge, all-suffix application,
complex residual state algebra, and both gradient streams. Local coarse
blocks are evaluated from elementary words; the residual LCU coins remain
ideal. This is not a full native phase-table emitter or a scaling test.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf import native_coarse_fixture as native
from compiler_robust_hopf.complex_analysis import (
    complex_magnitude_gradient, complex_phase_gradient,
)
from compiler_robust_hopf.complex_coarse_decoder import (
    apply_phase_blocks, decode_complex_coarse_frame_histograms, phase_tree_angles,
)
from compiler_robust_hopf.frames import direct_real_frame, real_tree_data
from tests.test_native_coarse_qbp import _coarse_block
from tests.test_operator_source_compiler import _adjoint, _word_matrix
from tests.test_reference_state_qbp import H, _amplify, _observable, _state_lcu


ATOL = 8e-11
Y_READOUT = np.array([[1., -1j], [1., 1j]]) / np.sqrt(2.)
PHASE_UNITS = (9, -5, 7, 3, -9, 1, 6)


def _phase_fixture(height):
    """An unwrapped phase tree with nontrivial global phase and winding."""
    size, mean = 1 << height, .37 + 4 * np.pi
    angles = np.asarray(PHASE_UNITS[:size - 1]) * np.pi / 4
    phases = np.full(size, mean)
    for leaf in range(size):
        node = 1
        for depth in range(height):
            side = (leaf >> (height - depth - 1)) & 1
            phases[leaf] += (2 * side - 1) * angles[node - 1]
            node = 2 * node + side
    return mean, angles, phases


def _angle_cases(height):
    size = 1 << height
    return (np.ones(size - 1, dtype=int),
            np.array([1, -1, 0, 2, 0, 1, -1][:size - 1]),
            np.array([0, 1, -1, 0, 2, -1, 1][:size - 1]))


def _native_blocks(height, units):
    """Return literal native phase rows and actual coarse real-tree rows."""
    k_word = [(name, *qubits) for name, qubits in native.commutator_word(target=0)]
    rz_step = [('X', 0), ('TDG', 0), ('X', 0), ('T', 0)]
    ry_step = [('Z', 0), ('H', 0)]
    phase_blocks, real_blocks = [], []
    for value in PHASE_UNITS[:(1 << height) - 1]:
        rotation = (rz_step if value >= 0 else _adjoint(rz_step)) * abs(value)
        phase_blocks.append(_word_matrix(1, k_word + rotation))
    for node, value in enumerate(units):
        rotation = (ry_step if value >= 0 else _adjoint(ry_step)) * abs(int(value))
        real_blocks.append(_word_matrix(1, (k_word if node == 0 else []) + rotation))
    return np.asarray(phase_blocks), np.asarray(real_blocks)


def _assemble_blocks(height, blocks, *, all_suffixes):
    """Independent full logical matrix; selectors are checked as binary prefixes."""
    size = 1 << height
    result = np.eye(size, dtype=complex)
    for node, block in enumerate(blocks, start=1):
        depth = node.bit_length() - 1
        prefix, target = node - (1 << depth), height - depth - 1
        layer = np.eye(size, dtype=complex)
        for label in range(size):
            if label >> (target + 1) != prefix or (label >> target) & 1:
                continue
            if not all_suffixes and label & ((1 << target) - 1):
                continue
            pair = [label, label | (1 << target)]
            layer[np.ix_(pair, pair)] = block
        result = layer @ result
    return result


def _walsh(height):
    result = np.array([[1.]])
    for _ in range(height):
        result = np.kron(result, H)
    return result


def _fixture(height, units):
    size = 1 << height
    mean, phase_angles, phases = _phase_fixture(height)
    angles = np.asarray(units) * np.pi / 4
    phase_blocks, real_blocks = _native_blocks(height, units)
    phase_coarse = _assemble_blocks(height, phase_blocks, all_suffixes=True)
    real_coarse = _assemble_blocks(height, real_blocks, all_suffixes=False)
    coarse = phase_coarse @ real_coarse
    gauge = np.diag(np.exp(1j * (phases - mean)))
    data = real_tree_data(angles)
    return dict(size=size, mean=mean, phase_angles=phase_angles, phases=phases,
                angles=angles, phase_blocks=phase_blocks, real_blocks=real_blocks,
                phase_coarse=phase_coarse, real_coarse=real_coarse, coarse=coarse,
                gauge=gauge, data=data, state=gauge @ data.state)


class ComplexCoarseQBPTests(unittest.TestCase):
    def test_phase_input_contract_overflow_and_decoder_nonmutation(self):
        for phases in (np.array([1 + 0j, 0], dtype=object),
                       np.array([np.inf, 0.]), np.array([np.nan, 0.])):
            with self.subTest(phases=phases), self.assertRaises(ValueError):
                phase_tree_angles(phases)
        largest = np.finfo(float).max
        with np.errstate(over='raise', invalid='raise'):
            for phases, expected_mean, expected_angle in (
                    ([largest, largest], largest, 0.),
                    ([-largest, largest], 0., largest)):
                mean, angles = phase_tree_angles(phases)
                self.assertTrue(np.isfinite(mean) and np.all(np.isfinite(angles)))
                self.assertEqual(mean, expected_mean)
                np.testing.assert_array_equal(angles, [expected_angle])

        fixture = _fixture(2, _angle_cases(2)[1])
        inputs = [fixture[key] for key in ('angles', 'real_blocks', 'phase_blocks', 'phases')]
        inputs += [np.array([1, 0, -1, 0]), np.array([0, 1, 0, -1])]
        saved = [value.copy() for value in inputs]
        for value in inputs:
            value.flags.writeable = False
        result = decode_complex_coarse_frame_histograms(*inputs, 4)
        self.assertTrue(np.all(np.isfinite(result)))
        for value, original in zip(inputs, saved):
            np.testing.assert_array_equal(value, original)
            self.assertFalse(value.flags.writeable)
        for index, invalid in ((3, inputs[3][:2]),
                               (1, inputs[1][:-1]), (2, inputs[2][:-1])):
            bad_inputs = list(inputs)
            bad_inputs[index] = invalid
            with self.subTest(input_index=index), self.assertRaises(ValueError):
                decode_complex_coarse_frame_histograms(*bad_inputs, 4)

    def test_subtree_phase_factorization_keeps_winding_and_global_gauge(self):
        for height in (2, 3):
            size = 1 << height
            expected_mean, expected_angles, phases = _phase_fixture(height)
            saved = phases.copy()
            mean, angles = phase_tree_angles(phases)
            self.assertAlmostEqual(mean, expected_mean, delta=ATOL)
            np.testing.assert_allclose(angles, expected_angles, atol=ATOL, rtol=0)
            np.testing.assert_array_equal(phases, saved)
            self.assertGreater(np.max(abs(phases)), 2 * np.pi)
            for variant in (phases, phases + .61 + 3 * np.pi,
                            phases + 2 * np.pi * np.arange(size)):
                common, splits = phase_tree_angles(variant)
                rows = [np.diag([np.exp(-1j * delta), np.exp(1j * delta)]) for delta in splits]
                product = _assemble_blocks(height, rows, all_suffixes=True)
                np.testing.assert_allclose(np.exp(1j * common) * product,
                                           np.diag(np.exp(1j * variant)), atol=ATOL, rtol=0)
            shifted_mean, shifted_angles = phase_tree_angles(phases + .61 + 3 * np.pi)
            self.assertAlmostEqual(shifted_mean - mean, .61 + 3 * np.pi, delta=ATOL)
            np.testing.assert_allclose(shifted_angles, angles, atol=ATOL, rtol=0)

    def test_literal_phase_rows_all_suffixes_and_complex_residual_amplification(self):
        for height in (2, 3):
            fixture = _fixture(height, _angle_cases(height)[1])
            size, coarse = fixture['size'], fixture['coarse']
            k = _coarse_block()
            for delta, actual in zip(fixture['phase_angles'], fixture['phase_blocks'], strict=True):
                rz = np.diag([np.exp(-1j * delta), np.exp(1j * delta)])
                np.testing.assert_allclose(actual, rz @ k, atol=ATOL, rtol=0)
                np.testing.assert_allclose(actual.conj().T @ actual, np.eye(2), atol=ATOL, rtol=0)
            self.assertGreater(np.max(abs(fixture['phase_blocks'][:, 0, 1])), 1e-3)
            for leaf in range(size):
                vector = np.eye(size, dtype=complex)[:, leaf]
                saved = vector.copy()
                actual = apply_phase_blocks(vector, fixture['phase_blocks'])
                np.testing.assert_allclose(actual, fixture['phase_coarse'] @ vector, atol=ATOL, rtol=0)
                np.testing.assert_array_equal(vector, saved)
            # A Hopf anchor/marker-only application omits required suffix copies.
            incorrect = _assemble_blocks(height, fixture['phase_blocks'], all_suffixes=False)
            self.assertGreater(np.linalg.norm(incorrect - fixture['phase_coarse'], ord=2), .5)
            target_frame = fixture['gauge'] @ direct_real_frame(height, fixture['angles'])
            error = np.linalg.norm(coarse - target_frame, ord=2)
            self.assertLess(error, min(1 / 64, 1 / (4 * np.sqrt(size))))
            self.assertLessEqual(error, (height + 1) * np.linalg.norm(k - np.eye(2), ord=2) + ATOL)
            q, residual = _state_lcu(coarse, fixture['state'])
            np.testing.assert_allclose(q[:size, 0], residual / 2, atol=ATOL, rtol=0)
            output = np.kron(np.eye(4), coarse) @ _amplify(q, size)[:, 0]
            expected = np.zeros(4 * size, dtype=complex)
            expected[:size] = fixture['state']
            np.testing.assert_allclose(output, expected, atol=ATOL, rtol=0)

    def test_complex_magnitude_xy_means_and_actual_block_histogram(self):
        for height in (2, 3):
            size, walsh = 1 << height, _walsh(height)
            observable = _observable(size)
            records = [(setting, leaf, (-1 if (leaf + setting) % 3 == 0 else 1))
                       for leaf in range(size) for setting in (0, 1)]
            records += [(0, 0, 1), (1, size - 1, -1)]
            hx, hy = np.zeros(size, dtype=int), np.zeros(size, dtype=int)
            for setting, leaf, sign in records:
                (hx if setting == 0 else hy)[leaf] += sign
            for units in _angle_cases(height):
                f = _fixture(height, units)
                transform = walsh @ f['coarse'].conj().T
                derivatives = (f['gauge'] @ np.asarray(f['data'].derivatives).T).T
                transformed = (transform @ derivatives.T).T
                joint = np.concatenate((transform @ f['coarse'][:, 0],
                                        transform @ observable @ f['state'])) / np.sqrt(2)
                probabilities = [abs(np.kron(readout, np.eye(size)) @ joint).reshape(2, size) ** 2
                                 for readout in (H, Y_READOUT)]
                decoded = 2 * np.sqrt(size) * (
                    transformed.real @ (probabilities[0][0] - probabilities[0][1])
                    + transformed.imag @ (probabilities[1][0] - probabilities[1][1]))
                analytic = complex_magnitude_gradient(f['angles'], f['phases'], observable)
                np.testing.assert_allclose(decoded, analytic, atol=ATOL, rtol=0)
                scale = 1.3
                expected = sum(4 * scale * np.sqrt(size) * sign *
                               (transformed[:, leaf].real if setting == 0 else transformed[:, leaf].imag)
                               for setting, leaf, sign in records) / len(records)
                actual = decode_complex_coarse_frame_histograms(
                    f['angles'], f['real_blocks'], f['phase_blocks'], f['phases'],
                    hx, hy, len(records), coefficient_scale=scale)
                np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                shifted = decode_complex_coarse_frame_histograms(
                    f['angles'], f['real_blocks'], f['phase_blocks'], f['phases'] + .61 + 3 * np.pi,
                    hx, hy, len(records), coefficient_scale=scale)
                np.testing.assert_allclose(shifted, actual, atol=ATOL, rtol=0)

    def test_direct_phase_y_sign_singular_leaves_and_wrong_gauge_witnesses(self):
        for height in (2, 3):
            size, observable = 1 << height, _observable(1 << height)
            for units in _angle_cases(height):
                f = _fixture(height, units)
                response = observable @ f['state']
                joint = np.concatenate((f['state'], response)) / np.sqrt(2)
                probabilities = abs(np.kron(Y_READOUT, np.eye(size)) @ joint).reshape(2, size) ** 2
                decoded = 2 * (probabilities[0] - probabilities[1])
                analytic = complex_phase_gradient(f['angles'], f['phases'], observable)
                np.testing.assert_allclose(decoded, analytic, atol=ATOL, rtol=0)
                self.assertAlmostEqual(decoded.sum(), 0., delta=ATOL)
                inactive = abs(f['state']) < 1e-12
                np.testing.assert_allclose(decoded[inactive], 0., atol=ATOL, rtol=0)
                np.testing.assert_allclose(complex_phase_gradient(f['angles'], f['phases'] + 5.2, observable),
                                           analytic, atol=ATOL, rtol=0)
            # Inconsistent gauge correction is different from a legitimate
            # common phase change in both the target and its derivatives.
            f = _fixture(height, _angle_cases(height)[1])
            response = observable @ f['state']
            derivatives = np.asarray(f['data'].derivatives)
            exact = complex_magnitude_gradient(f['angles'], f['phases'], observable)
            omitted_phase = 2 * np.real(derivatives @ response)
            wrong_common = 2 * np.real((np.exp(1j * f['phases']) * derivatives).conj() @ response)
            self.assertGreater(np.max(abs(omitted_phase - exact)), .1)
            self.assertGreater(np.max(abs(wrong_common - exact)), .01)


if __name__ == '__main__':
    unittest.main()
