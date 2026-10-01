"""Complete small complex-chart words with arbitrary dirty-helper inputs.

Only 64x4 preparation/operator batches and 64x2 readout columns are
propagated. Independent logical oracles use the analytic Hopf chart,
matrix-defined K3, and explicit all-suffix prefix pairs. This is a finite
exact-preparation fallback, not an emitted general residual-table compiler.
"""
from __future__ import annotations

from collections import Counter
import unittest

import numpy as np

from compiler_robust_hopf import native_complex_coarse_fixture as native
from compiler_robust_hopf.complex_analysis import (
    complex_magnitude_gradient, complex_phase_gradient,
)
from compiler_robust_hopf.complex_coarse_decoder import decode_complex_coarse_frame_histograms
from compiler_robust_hopf.conventions import marker_label
from compiler_robust_hopf.frames import direct_real_frame, hopf_ry, real_tree_data
from tests.test_native_coarse_qbp import (
    ALPHABET, ATOL, CASES, H, T, X, WALSH, _coarse_block, _embedding,
    _operator_error, _propagate, _score_operator,
)


OBSERVABLES = ('tilted_x', 'hadamard_low')
PHASE_UNITS = (9, -5, 7)
PHASE_MEAN = .37 + 4 * np.pi


def _phase_oracle():
    """Matrix-order definition independent of emitted words/logical_blocks."""
    root, left, right = np.asarray(PHASE_UNITS) * np.pi / 4
    offsets = np.array([-root - left, -root + left, root - right, root + right])
    phases = PHASE_MEAN + offsets
    ideal, coarse = np.eye(4, dtype=complex), np.eye(4, dtype=complex)
    native_rows = []
    for angle, pairs in zip((root, left, right),
                            (((0, 2), (1, 3)), ((0, 1),), ((2, 3),)), strict=True):
        row = np.diag([np.exp(-1j * angle), np.exp(1j * angle)])
        actual = row @ _coarse_block()
        native_rows.append(actual)
        for block, accumulator in ((row, ideal), (actual, coarse)):
            layer = np.eye(4, dtype=complex)
            for pair in pairs:
                layer[np.ix_(pair, pair)] = block
            accumulator[:] = layer @ accumulator
    return phases, ideal, coarse, np.asarray(native_rows)


def _oracle(units, observable='tilted_x'):
    angles = np.asarray(units) * np.pi / 4
    real_frame = direct_real_frame(2, angles)
    phases, gauge, phase_coarse, phase_blocks = _phase_oracle()
    error = np.eye(4, dtype=complex)
    error[np.ix_((0, 2), (0, 2))] = _coarse_block()
    frame, coarse = gauge @ real_frame, phase_coarse @ real_frame @ error
    data = real_tree_data(angles)
    state = gauge @ data.state
    derivatives = (gauge @ np.asarray(data.derivatives).T).T
    center = np.kron(T @ X @ T.conj().T, np.eye(2)) if observable == 'tilted_x' else np.kron(np.eye(2), H)
    operator = real_frame @ center @ real_frame.conj().T
    transformed = (WALSH @ coarse.conj().T @ derivatives.T).T
    # These analytic routines use the original physical phase tuple, so
    # their answer independently checks cancellation of the common gauge.
    magnitude = complex_magnitude_gradient(angles, phases, operator)
    phase = complex_phase_gradient(angles, phases, operator)
    return dict(angles=angles, real_frame=real_frame, frame=frame, coarse=coarse,
                gauge=gauge, phases=phases, phase_blocks=phase_blocks, data=data,
                state=state, derivatives=derivatives, operator=operator,
                transformed=transformed, magnitude=magnitude, phase=phase)


def _full_port_error(word, active, branch_value=None):
    errors = []
    for dirty in (0, 1):
        for branch in (0, 1):
            rows = np.arange(4) | (branch << 2) | (dirty << 5)
            columns, expected = np.zeros((64, 4), complex), np.zeros((64, 4), complex)
            columns[rows] = np.eye(4)
            expected[rows] = active if branch_value is None or branch == branch_value else np.eye(4)
            errors.append(_propagate(word, columns) - expected)
    # Larger assembled metrics do not increase the propagated batch size.
    return _operator_error(np.concatenate(errors, axis=1))


class NativeComplexCoarseQBPTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.initial = _embedding()
        cls.zero_branch = cls.initial[:, [0, 2]]
        cls.outputs, cls.words = {}, {}
        for units in CASES:
            for observable in OBSERVABLES:
                words = {
                    'X': native.magnitude_protocol_word(units, observable, 'X'),
                    'Y': native.magnitude_protocol_word(units, observable, 'Y'),
                    'phase': native.phase_protocol_word(units, observable),
                    'original_mag': native.original_magnitude_protocol_word(units, observable),
                    'original_phase': native.original_phase_protocol_word(units, observable),
                }
                seen = {}
                for mode, word in words.items():
                    if word not in seen:
                        seen[word] = _propagate(word, cls.zero_branch)
                    cls.outputs[units, observable, mode] = seen[word]
                    cls.words[units, observable, mode] = word

    def test_complete_prefix_selection_and_observable_phase_on_every_input(self):
        _, ideal, coarse, _ = _phase_oracle()
        np.testing.assert_allclose(ideal, np.diag([-1, 1j, 1j, 1]), atol=ATOL, rtol=0)
        self.assertLess(np.linalg.norm(_coarse_block() - np.eye(2), ord=2), 1 / 200)
        for is_coarse, target in ((False, ideal), (True, coarse)):
            for branch in (None, 0, 1):
                word = native.phase_table_word(coarse=is_coarse, branch_value=branch)
                self.assertLess(_full_port_error(word, target, branch), ATOL)
                self.assertFalse(any(q in (3, 4) for _, wires in word for q in wires))
        for observable in OBSERVABLES:
            oracle = _oracle(CASES[0], observable)
            self.assertLess(_full_port_error(native.controlled_observable_word(CASES[0], observable),
                                             oracle['operator'], 1), ATOL)

    def test_joint_preparation_all_dirty_columns_and_literal_resource_words(self):
        for units in CASES:
            oracle = _oracle(units)
            self.assertLess(np.linalg.norm(oracle['coarse'] - oracle['frame'], ord=2), 3 / 200)
            word = native.preparation_word(units)
            expected = np.zeros_like(self.initial)
            for dirty in (0, 1):
                for branch in (0, 1):
                    rows = np.arange(4) | (branch << 2) | (dirty << 5)
                    expected[rows, branch + 2 * dirty] = oracle['coarse' if branch == 0 else 'frame'][:, 0]
            output = _propagate(word, self.initial)
            self.assertLess(_operator_error(output - expected), ATOL)
            np.testing.assert_allclose(output.conj().T @ output, np.eye(4), atol=ATOL, rtol=0)
            first_query = next(i for i, gate in enumerate(word) if gate == ('CX', (0, 5)))
            intermediate = _propagate(word[:first_query + 1], self.initial)
            moved = sum(np.linalg.norm(intermediate[((np.arange(64) >> 5) & 1) != column // 2, column]) ** 2
                        for column in range(4))
            self.assertGreater(moved, 3.9)
            self.assertFalse(any(q in (3, 4) for _, wires in word for q in wires))

        for word in self.words.values():
            self.assertIsInstance(word, tuple)
            counts = Counter(name for name, _ in word)
            for name, wires in word:
                self.assertIn(name, ALPHABET)
                self.assertEqual(len(wires), 2 if name == 'CX' else 1)
                self.assertEqual(len(set(wires)), len(wires))
                self.assertTrue(all(0 <= q < 6 and q not in (3, 4) for q in wires))
            ledger = native.gate_counts(word)
            self.assertEqual(ledger['total'], len(word))
            self.assertEqual(ledger['T_count'], counts['T'] + counts['TDG'])
            self.assertEqual(ledger['Clifford_count'], len(word) - ledger['T_count'])
            for name in ALPHABET:
                self.assertEqual(ledger[name], counts[name])
        # Both phase executions use the same exact finite-size preparation.
        self.assertEqual(self.words[CASES[0], OBSERVABLES[0], 'phase'],
                         self.words[CASES[0], OBSERVABLES[0], 'original_phase'])
        # Literal supplied-word counts, not optimized synthesis estimates.
        for mode, expected in (('X', (25922, 74677)), ('original_mag', (66, 597)),
                               ('phase', (40, 350))):
            counts = native.gate_counts(self.words[CASES[0], 'tilted_x', mode])
            self.assertEqual((counts['T_count'], counts['Clifford_count']), expected)
        preparation = native.gate_counts(native.preparation_word())
        self.assertEqual((preparation['T_count'], preparation['Clifford_count']), (20506, 48630))

    def test_both_streams_and_original_protocol_have_full_dirty_score_operators(self):
        for units in CASES:
            for observable in OBSERVABLES:
                oracle = _oracle(units, observable)
                response = oracle['operator'] @ oracle['state']
                np.testing.assert_allclose(2 * np.real(oracle['derivatives'].conj() @ response),
                                           oracle['magnitude'], atol=ATOL, rtol=0)
                np.testing.assert_allclose(2 * np.imag(oracle['state'].conj() * response),
                                           oracle['phase'], atol=ATOL, rtol=0)
                for node in range(1, 4):
                    corrected = sum(_score_operator(self.outputs[units, observable, setting],
                                                    8 * (oracle['transformed'][node - 1].real if setting == 'X'
                                                         else oracle['transformed'][node - 1].imag))
                                    for setting in ('X', 'Y')) / 2
                    wanted = oracle['magnitude'][node - 1] * np.eye(2)
                    np.testing.assert_allclose(corrected, wanted, atol=ATOL, rtol=0)
                    marker = marker_label(node, 2)
                    character = np.array([(-1) ** ((leaf & marker).bit_count()) for leaf in range(4)])
                    weights = 2 * oracle['data'].incoming_amplitude[node - 1] * character
                    original = _score_operator(self.outputs[units, observable, 'original_mag'], weights)
                    np.testing.assert_allclose(original, wanted, atol=ATOL, rtol=0)
                for leaf in range(4):
                    for mode in ('phase', 'original_phase'):
                        operator = _score_operator(self.outputs[units, observable, mode], 2 * np.eye(4)[leaf])
                        np.testing.assert_allclose(operator, oracle['phase'][leaf] * np.eye(2), atol=ATOL, rtol=0)
                self.assertAlmostEqual(float(np.sum(oracle['phase'])), 0., delta=ATOL)
        singular = _oracle((0, 1, -1), 'hadamard_low')
        np.testing.assert_allclose(singular['magnitude'], [0, np.sqrt(2), 0], atol=ATOL, rtol=0)
        np.testing.assert_allclose(singular['phase'], [-1 / np.sqrt(2), 1 / np.sqrt(2), 0, 0], atol=ATOL, rtol=0)

    def test_integer_histogram_uses_actual_prefix_blocks_and_physical_phases(self):
        records = [('X', 0, 1), ('Y', 0, -1), ('X', 1, -1), ('Y', 1, 1),
                   ('X', 2, 1), ('Y', 3, 1), ('X', 3, 1), ('Y', 1, -1), ('X', 0, 1)]
        hx, hy = np.zeros(4, dtype=int), np.zeros(4, dtype=int)
        for setting, leaf, sign in records:
            (hx if setting == 'X' else hy)[leaf] += sign
        scale = 1.7
        for units in CASES:
            oracle = _oracle(units)
            a = oracle['angles']
            real_blocks = np.array([hopf_ry(a[0]) @ _coarse_block(), hopf_ry(a[1]), hopf_ry(a[2])])
            expected = sum(8 * scale * sign * (oracle['transformed'][:, leaf].real if setting == 'X'
                                              else oracle['transformed'][:, leaf].imag)
                           for setting, leaf, sign in records) / len(records)
            actual = decode_complex_coarse_frame_histograms(
                a, real_blocks, oracle['phase_blocks'], oracle['phases'], hx, hy,
                len(records), coefficient_scale=scale)
            np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
            blocks = native.logical_blocks(units)
            np.testing.assert_allclose(blocks['C'], oracle['coarse'], atol=ATOL, rtol=0)
            np.testing.assert_allclose(blocks['V'], oracle['frame'], atol=ATOL, rtol=0)
            np.testing.assert_allclose(blocks['coarse_blocks'], real_blocks, atol=ATOL, rtol=0)
            np.testing.assert_allclose(blocks['phase_blocks'], oracle['phase_blocks'], atol=ATOL, rtol=0)

    def test_relative_branch_phase_and_wrong_classical_gauge_are_detected(self):
        units, observable = CASES[0], 'tilted_x'
        oracle = _oracle(units, observable)
        changed = []
        for setting in ('X', 'Y'):
            word = self.words[units, observable, setting]
            position = word.index(('H', (2,))) + 1
            altered = word[:position] + (('Z', (2,)),) + word[position:]
            output = _propagate(altered, self.zero_branch)
            weights = 8 * (oracle['transformed'][0].real if setting == 'X' else oracle['transformed'][0].imag)
            changed.append(_score_operator(output, weights))
        np.testing.assert_allclose(sum(changed) / 2, -oracle['magnitude'][0] * np.eye(2), atol=ATOL, rtol=0)
        self.assertGreater(abs(oracle['magnitude'][0]), .5)

        phase_word = self.words[units, observable, 'phase']
        position = phase_word.index(('H', (2,))) + 1
        altered = phase_word[:position] + (('Z', (2,)),) + phase_word[position:]
        wrong_phase = _propagate(altered, self.zero_branch)
        leaf = int(np.argmax(abs(oracle['phase'])))
        np.testing.assert_allclose(_score_operator(wrong_phase, 2 * np.eye(4)[leaf]),
                                   -oracle['phase'][leaf] * np.eye(2), atol=ATOL, rtol=0)
        self.assertGreater(abs(oracle['phase'][leaf]), .3)

        transform = WALSH @ oracle['coarse'].conj().T
        bare = np.asarray(oracle['data'].derivatives)
        for wrong_derivatives in (bare, np.exp(1j * oracle['phases']) * bare):
            wrong_f = (transform @ wrong_derivatives.T).T
            decoded = sum(_score_operator(self.outputs[units, observable, setting],
                                         8 * (wrong_f[0].real if setting == 'X' else wrong_f[0].imag))
                          for setting in ('X', 'Y')) / 2
            self.assertGreater(np.linalg.norm(decoded - oracle['magnitude'][0] * np.eye(2), ord=2), .01)
        y_piece = _score_operator(self.outputs[units, observable, 'Y'],
                                  8 * oracle['transformed'][0].imag) / 2
        response = transform @ oracle['operator'] @ oracle['state']
        expected_y_piece = 2 * np.dot(oracle['transformed'][0].imag, response.imag)
        np.testing.assert_allclose(y_piece, expected_y_piece * np.eye(2), atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(y_piece, ord=2), .003)


if __name__ == '__main__':
    unittest.main()
