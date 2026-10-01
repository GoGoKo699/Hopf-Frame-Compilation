"""A bounded native integration audit with all branch/dirty input columns.

Only 64x4 preparation and 64x2 readout columns are propagated. Independent
ideal oracles use analytic Hopf frames and the explicit two-by-two coarse
commutator. This exact finite fixture is not a scalable compiler or an
advantage benchmark; its two reserved compiler flags remain unused.
"""
from __future__ import annotations

from collections import Counter
import unittest

import numpy as np

from compiler_robust_hopf import native_coarse_fixture as native
from compiler_robust_hopf.coarse_frame_decoder import decode_coarse_frame_histograms
from compiler_robust_hopf.conventions import marker_label
from compiler_robust_hopf.frames import direct_real_frame, hopf_ry, real_tree_data
from tests.test_operator_source_compiler import _apply_native_word, _word_matrix


ATOL = 3e-10
H = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2.)
T = np.diag([1., np.exp(1j * np.pi / 4)])
S = np.diag([1., 1j])
X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
WALSH = np.kron(H, H)
CASES = ((1, 1, 1), (0, 1, -1), (1, -1, 0))
ALPHABET = {'H', 'X', 'Z', 'S', 'SDG', 'T', 'TDG', 'CX'}


def _coarse_block():
    """Matrix-order definition, independent of either emitted K3 word."""
    block = T @ H @ T @ H @ T.conj().T @ H @ T.conj().T @ H
    for _ in range(3):
        conjugate = S @ block @ S.conj().T
        block = block @ conjugate @ block.conj().T @ conjugate.conj().T
    return block


def _oracle(units, observable='tilted_x'):
    angles = np.asarray(units) * np.pi / 4
    frame = direct_real_frame(2, angles)
    error = np.eye(4, dtype=complex)
    error[np.ix_([0, 2], [0, 2])] = _coarse_block()
    coarse = frame @ error
    center = T @ X @ T.conj().T if observable == 'tilted_x' else Y
    operator = frame @ np.kron(center, np.eye(2)) @ frame.conj().T
    data = real_tree_data(angles)
    derivatives = np.asarray(data.derivatives)
    transformed = (WALSH @ coarse.conj().T @ derivatives.T).T
    gradient = 2 * np.real(derivatives @ operator @ data.state)
    return angles, frame, coarse, data, transformed, gradient


def _embedding():
    # Columns are branch + 2*dirty. System 0,1 and compiler flags 3,4
    # start at zero; neither the branch nor helper 5 is initialized.
    columns = np.zeros((64, 4), dtype=complex)
    for dirty in (0, 1):
        for branch in (0, 1):
            columns[(branch << 2) | (dirty << 5), branch + 2 * dirty] = 1
    return columns


def _propagate(word, columns):
    # Reuse the separately tested column simulator, not the emitter's own
    # simulator. Its legacy alphabet evaluates a Z as the identical S S.
    flat = []
    for name, qubits in word:
        flat.extend([('S', qubits[0])] * 2 if name == 'Z' else [(name, *qubits)])
    return _apply_native_word(6, flat, columns)


def _operator_error(columns):
    return np.sqrt(max(0., np.linalg.eigvalsh(columns.conj().T @ columns)[-1]))


def _score_operator(columns, coefficients):
    labels = np.arange(64)
    signs = 1 - 2 * ((labels >> 2) & 1)
    scores = signs * coefficients[labels & 3]
    return columns.conj().T @ (scores[:, None] * columns)


class NativeCoarseQBPTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.initial = _embedding()
        cls.zero_branch = cls.initial[:, [0, 2]]
        cls.readouts = {}
        for units in CASES:
            for observable in ('tilted_x', 'y'):
                for quadrature in ('X', 'Y'):
                    word = native.corrected_protocol_word(units, observable, quadrature)
                    cls.readouts[units, observable, quadrature] = _propagate(word, cls.zero_branch)

    def test_full_input_preparation_dirty_echo_and_literal_gate_ledger(self):
        block = _coarse_block()
        self.assertLess(np.linalg.norm(block - np.eye(2), ord=2), 1 / 64)
        self.assertGreater(np.linalg.norm(block.imag), 1e-3)
        small_word = [(name, *qubits) for name, qubits in native.commutator_word(target=0)]
        np.testing.assert_allclose(_word_matrix(1, small_word), block, atol=ATOL, rtol=0)
        toffoli = native.toffoli_word(0, 1, 2)
        self.assertEqual(len(toffoli), 15)
        self.assertEqual(sum(name in ('T', 'TDG') for name, _ in toffoli), 7)
        np.testing.assert_allclose(_word_matrix(3, [(name, *qubits) for name, qubits in toffoli]),
                                   _word_matrix(3, [('CCX', 0, 1, 2)]), atol=ATOL, rtol=0)

        # Test the entire selected E, including low=1 identity, on all
        # sixteen system/branch/helper inputs. Propagation stays in batches
        # of four; only the resulting error columns are concatenated.
        error = np.eye(4, dtype=complex)
        error[np.ix_([0, 2], [0, 2])] = block
        selected = native.selected_error_word(0)
        error_batches = []
        for dirty in (0, 1):
            for branch in (0, 1):
                rows = np.arange(4) | (branch << 2) | (dirty << 5)
                inputs = np.zeros((64, 4), dtype=complex)
                target = np.zeros_like(inputs)
                inputs[rows] = np.eye(4)
                target[rows] = error if branch == 0 else np.eye(4)
                error_batches.append(_propagate(selected, inputs) - target)
        self.assertLess(_operator_error(np.concatenate(error_batches, axis=1)), ATOL)

        for units in CASES:
            _, frame, coarse, _, _, _ = _oracle(units)
            word = native.preparation_word(units)
            output = _propagate(word, self.initial)
            target = np.zeros_like(output)
            for dirty in (0, 1):
                for branch in (0, 1):
                    rows = np.arange(4) | (branch << 2) | (dirty << 5)
                    target[rows, branch + 2 * dirty] = (coarse if branch == 0 else frame)[:, 0]
            self.assertLess(_operator_error(output - target), ATOL)
            np.testing.assert_allclose(output.conj().T @ output, np.eye(4), atol=ATOL, rtol=0)
            np.testing.assert_allclose(native.simulate_columns(word, self.initial), output, atol=ATOL, rtol=0)
            # The initial query flips both possible helper values. Only the
            # completed reflection echoes return it; it is not a spectator.
            first_query = next(i for i, gate in enumerate(word) if gate == ('CX', (0, 5)))
            intermediate = _propagate(word[:first_query + 1], self.initial)
            moved = sum(np.linalg.norm(intermediate[((np.arange(64) >> 5) & 1) != (column // 2), column]) ** 2
                        for column in range(4))
            self.assertGreater(moved, 3.9)
            self.assertFalse(any(qubit in (3, 4) for _, qubits in word for qubit in qubits))

        words = [native.corrected_protocol_word(), native.corrected_protocol_word(readout='Y'),
                 native.original_protocol_word(), native.preparation_word()]
        for word in words:
            self.assertIsInstance(word, tuple)
            count = Counter(name for name, _ in word)
            for name, qubits in word:
                self.assertIn(name, ALPHABET)
                self.assertEqual(len(qubits), 2 if name == 'CX' else 1)
                self.assertEqual(len(set(qubits)), len(qubits))
                self.assertTrue(all(0 <= qubit < 6 for qubit in qubits))
            ledger = native.gate_counts(word)
            self.assertEqual(ledger['total'], len(word))
            self.assertEqual(ledger['T_count'], count['T'] + count['TDG'])
            self.assertEqual(ledger['Clifford_count'], len(word) - count['T'] - count['TDG'])
            for name in ALPHABET:
                self.assertEqual(ledger[name], count[name])
        # Literal unsimplified benchmark, not an advantage or optimality claim.
        self.assertEqual((native.gate_counts(words[0])['total'], native.gate_counts(words[0])['T_count']),
                         (22471, 5914))
        self.assertEqual((native.gate_counts(words[2])['total'], native.gate_counts(words[2])['T_count']),
                         (231, 26))
        self.assertEqual(len(words[1]), len(words[0]) + 1)

    def test_native_score_operators_cover_every_dirty_input_and_original_comparison(self):
        for units in CASES:
            for observable in ('tilted_x', 'y'):
                _, _, _, data, transformed, gradient = _oracle(units, observable)
                np.testing.assert_allclose(gradient, [np.sqrt(2), 0, 0] if observable == 'tilted_x' else 0,
                                           atol=ATOL, rtol=0)
                original = _propagate(native.original_protocol_word(units, observable), self.zero_branch)
                for node in range(1, 4):
                    corrected = sum(_score_operator(self.readouts[units, observable, quadrature],
                                                   8 * (transformed[node - 1].real if quadrature == 'X'
                                                        else transformed[node - 1].imag))
                                    for quadrature in ('X', 'Y')) / 2
                    # Equality of these 2x2 score operators covers coherent
                    # dirty inputs and their arbitrary external references.
                    np.testing.assert_allclose(corrected, gradient[node - 1] * np.eye(2), atol=ATOL, rtol=0)
                    marker = marker_label(node, 2)
                    character = np.array([(-1) ** ((x & marker).bit_count()) for x in range(4)])
                    original_score = _score_operator(original, 2 * data.incoming_amplitude[node - 1] * character)
                    np.testing.assert_allclose(original_score, corrected, atol=ATOL, rtol=0)

    def test_histogram_decoder_matches_independent_native_coarse_scores(self):
        records = [('X', 0, 1), ('Y', 0, -1), ('X', 1, -1), ('Y', 1, 1),
                   ('X', 2, 1), ('Y', 3, 1), ('X', 3, 1), ('Y', 1, -1), ('X', 0, 1)]
        hist_x, hist_y = np.zeros(4, dtype=int), np.zeros(4, dtype=int)
        for quadrature, leaf, sign in records:
            (hist_x if quadrature == 'X' else hist_y)[leaf] += sign
        scale = 1.7
        for units in CASES:
            angles, _, _, _, transformed, _ = _oracle(units)
            expected_blocks = np.array([hopf_ry(angles[0]) @ _coarse_block(),
                                        hopf_ry(angles[1]), hopf_ry(angles[2])])
            blocks = native.logical_blocks(units)['coarse_blocks']
            np.testing.assert_allclose(blocks, expected_blocks, atol=ATOL, rtol=0)
            expected = sum(8 * scale * sign * (transformed[:, leaf].real if quadrature == 'X'
                                             else transformed[:, leaf].imag)
                           for quadrature, leaf, sign in records) / len(records)
            actual = decode_coarse_frame_histograms(angles, blocks, hist_x, hist_y, len(records),
                                                    coefficient_scale=scale)
            np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)

    def test_relative_phase_and_missing_y_are_detected(self):
        units = CASES[0]
        _, _, _, _, transformed, gradient = _oracle(units)
        wrong = []
        for quadrature in ('X', 'Y'):
            word = native.corrected_protocol_word(units, 'tilted_x', quadrature)
            self.assertEqual(word[0], ('H', (2,)))
            # A relative minus after preparing the branch is not a harmless
            # common global phase; it reverses all interference means.
            altered = word[:1] + (('Z', (2,)),) + word[1:]
            output = _propagate(altered, self.zero_branch)
            weights = 8 * (transformed[0].real if quadrature == 'X' else transformed[0].imag)
            wrong.append(_score_operator(output, weights))
        reversed_mean = sum(wrong) / 2
        np.testing.assert_allclose(reversed_mean, -gradient[0] * np.eye(2), atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(reversed_mean - gradient[0] * np.eye(2), ord=2), 2)

        x_piece = _score_operator(self.readouts[units, 'y', 'X'], 8 * transformed[0].real) / 2
        y_piece = _score_operator(self.readouts[units, 'y', 'Y'], 8 * transformed[0].imag) / 2
        self.assertGreater(np.linalg.norm(x_piece, ord=2), 1e-3)
        np.testing.assert_allclose(x_piece + y_piece, 0, atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
