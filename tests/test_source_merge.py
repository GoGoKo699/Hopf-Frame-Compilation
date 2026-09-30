"""Bounded diagnostics for scalar-source flag merging and routed leakage.

These are scalar coefficient filters, not a whole-frame half-unitary
construction. Every routed stage retains its own charged source calls.
The native m=3 fixture checks all dirty input columns on three noncommuting
Hopf edges. Exhausting two-flag Pauli labels verifies only that routing
class; it is not a lower bound for arbitrary two-clean compilation.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

try:
    from .test_operator_source_compiler import (
        X, Y, Z, _adjoint, _expand_toffolis, _operator_source,
        _scalar_block, _source_word, _word_matrix,
    )
except ImportError:
    from test_operator_source_compiler import (
        X, Y, Z, _adjoint, _expand_toffolis, _operator_source,
        _scalar_block, _source_word, _word_matrix,
    )


ATOL = 4e-11


def _native_scalar_word(core_width, signs):
    """Literal S_f word; core bits are low and the scalar flag is high."""
    flag = core_width
    source = _source_word(core_width)
    controlled = _adjoint(source) + [('CX', flag, 0)] + source
    negative = [('X', flag)] + controlled + [('X', flag)]
    uncontrolled = _adjoint(source) + [('X', 0)] + source
    mask = [('Z', j) for j, sign in enumerate(signs) if sign]
    # Matrix convention: S=H C_0(M) N C_1(M) H.
    return ([('H', flag)] + controlled + mask + uncontrolled + mask
            + negative + [('H', flag)])


def _native_hopf_quarter_turn(low, high):
    """Literal R_y(pi/2) on one three-qubit addressed Hopf edge."""
    difference = low ^ high
    assert difference and difference & (difference - 1) == 0
    target = difference.bit_length() - 1
    assert (low >> target) & 1 == 0
    controls = [j for j in range(3) if j != target]
    negative = [('X', j) for j in controls if not ((low >> j) & 1)]
    controlled_x = [('CCX', *controls, target)]
    controlled_z = [('H', target)] + controlled_x + [('H', target)]
    word = negative + controlled_z + controlled_x + negative[::-1]
    return _word_matrix(3, _expand_toffolis(word))


def _routed_stage(pauli, logical_rotation, projector, coefficient, odd):
    """Two flags, logical register, dirty core; full unitary input space."""
    core_dimension = odd.shape[0]
    logical_dimension = projector.shape[0]
    accepted = logical_rotation @ (
        np.eye(logical_dimension) + (coefficient - 1) * projector)
    rejected = logical_rotation @ projector
    actual = (np.kron(np.eye(4), np.kron(accepted, np.eye(core_dimension)))
              + np.kron(pauli, np.kron(rejected, odd)))
    return actual, accepted, rejected


class SourceMergeTests(unittest.TestCase):
    def test_direct_merged_flag_native_word_all_dirty_inputs(self):
        """One flag cannot directly supply both coefficient selection roles."""
        width, target, flag = 3, 3, 4
        source, _, amplitudes = _operator_source(width)
        source_word = _source_word(width)
        controlled = _adjoint(source_word) + [('CX', flag, 0)] + source_word
        negative = [('X', flag)] + controlled + [('X', flag)]
        phase_y = X @ Z
        identity_core = np.eye(8)
        identity_data = np.eye(16)
        # Chronological CZ then CX gives the literal controlled XZ.
        controlled_phase_y = [('H', target), ('CX', flag, target),
                              ('H', target), ('CX', flag, target)]
        identity_case = None
        for cosine_bits, sine_bits in itertools.product(
                itertools.product((0, 1), repeat=width), repeat=2):
            with self.subTest(cosine_bits=cosine_bits, sine_bits=sine_bits):
                masks = [[('Z', j) for j, bit in enumerate(bits) if bit]
                         for bits in (cosine_bits, sine_bits)]
                # Flag 0: N_c then M. Flag 1: M then N_s, then target XZ.
                # Each controlled source is an actual Clifford+T word;
                # surrounding each one by its mask controls the masked source.
                word = ([('H', flag)]
                        + masks[0] + negative + masks[0] + negative
                        + controlled + masks[1] + controlled + masks[1]
                        + controlled_phase_y + [('H', flag)])
                circuit = _word_matrix(5, word)
                np.testing.assert_allclose(circuit.conj().T @ circuit,
                                           np.eye(32), atol=ATOL, rtol=0)
                accepted = circuit[:16, :16]
                coefficients = [1 - 2 * np.dot(amplitudes ** 2, bits)
                                for bits in (cosine_bits, sine_bits)]
                masked = [_word_matrix(width, mask) @ source
                          @ _word_matrix(width, mask) for mask in masks]
                odd_c, odd_s = [(source @ item - item @ source) / 2
                                for item in masked]
                cosine, sine = coefficients
                intended = np.kron((cosine * np.eye(2) + sine * phase_y) / 2,
                                   identity_core)
                unwanted = (np.kron(np.eye(2), odd_c)
                            - np.kron(phase_y, odd_s)) / 2
                np.testing.assert_allclose(accepted, intended + unwanted,
                                           atol=ATOL, rtol=0)
                error = accepted - intended
                squared = error.conj().T @ error
                target_trace = np.trace(squared.reshape(2, 8, 2, 8),
                                        axis1=0, axis2=2) / 2
                bound_squared = max(0, (2 - cosine ** 2 - sine ** 2) / 4)
                np.testing.assert_allclose(target_trace,
                                           bound_squared * identity_core,
                                           atol=ATOL, rtol=0)
                self.assertGreaterEqual(np.linalg.norm(error, 2) + ATOL,
                                        np.sqrt(bound_squared))
                if abs(cosine ** 2 + sine ** 2 - 1) < ATOL:
                    self.assertGreaterEqual(np.linalg.norm(error, 2) + ATOL, .5)
                if cosine_bits == (0, 0, 0) and sine_bits == (1, 0, 0):
                    identity_case = circuit, accepted

        self.assertIsNotNone(identity_case)
        circuit, accepted = identity_case
        np.testing.assert_allclose(accepted.conj().T, accepted,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(accepted @ accepted, accepted,
                                   atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.trace(accepted).real, 8, delta=ATOL)
        reflection = np.diag([-1] * 16 + [1] * 16)
        amplified = -circuit @ reflection @ circuit.conj().T @ reflection @ circuit
        np.testing.assert_allclose(amplified.conj().T @ amplified,
                                   np.eye(32), atol=ATOL, rtol=0)
        np.testing.assert_allclose(amplified[:16, :16], -accepted,
                                   atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(amplified[:16, :16] - identity_data, 2),
                               2, delta=ATOL)

    def test_native_scalar_word_and_source_parity(self):
        source, _, _ = _operator_source(3)
        for signs in itertools.product((0, 1), repeat=3):
            scalar, masked = _scalar_block(source, signs)
            native = _word_matrix(4, _native_scalar_word(3, signs))
            np.testing.assert_allclose(native, scalar, atol=ATOL, rtol=0)
            coefficient = np.trace(source @ masked).real / 8
            odd = (source @ masked - masked @ source) / 2
            expected = (coefficient * np.eye(16)
                        + np.kron(X, odd))
            np.testing.assert_allclose(native, expected, atol=ATOL, rtol=0)
            np.testing.assert_allclose(odd.conj().T, -odd, atol=ATOL, rtol=0)
            np.testing.assert_allclose(
                odd.conj().T @ odd,
                (1 - coefficient ** 2) * np.eye(8), atol=ATOL, rtol=0)
            np.testing.assert_allclose(
                source @ odd @ source, -odd, atol=ATOL, rtol=0)

    def test_all_two_flag_pauli_routes_with_zero_pair_returns(self):
        one_qubit = {'I': np.eye(2), 'X': X, 'Y': Y, 'Z': Z}
        candidates = []
        for names in itertools.product(one_qubit, repeat=2):
            matrix = np.kron(one_qubit[names[0]], one_qubit[names[1]])
            if abs(matrix[0, 0]) < ATOL:
                flip = tuple(int(name in ('X', 'Y')) for name in names)
                candidates.append((names, matrix, flip))
        self.assertEqual(len(candidates), 12)
        valid = 0
        for triple in itertools.product(candidates, repeat=3):
            _, p1, v1 = triple[0]
            _, p2, v2 = triple[1]
            _, p3, v3 = triple[2]
            pair_returns = [(p2 @ p1)[0, 0], (p3 @ p1)[0, 0],
                            (p3 @ p2)[0, 0]]
            zero_pairs = all(abs(value) < ATOL for value in pair_returns)
            self.assertEqual(zero_pairs, len({v1, v2, v3}) == 3)
            if not zero_pairs:
                continue
            valid += 1
            self.assertEqual(tuple(a ^ b ^ c for a, b, c in zip(v1, v2, v3)),
                             (0, 0))
            product = p3 @ p2 @ p1
            np.testing.assert_allclose(
                product, np.diag(np.diag(product)), atol=ATOL, rtol=0)
            self.assertAlmostEqual(abs(product[0, 0]), 1, delta=ATOL)
            # Independent Hermitian Pauli signs cannot alter these facts.
            for signs in itertools.product((-1, 1), repeat=3):
                phase = np.prod(signs) * product[0, 0]
                self.assertAlmostEqual(abs(phase), 1, delta=ATOL)
        self.assertEqual(valid, 6 * 4 ** 3)

    def test_three_noncommuting_hopf_edges_all_dirty_inputs(self):
        source, _, _ = _operator_source(3)
        signs = (0, 1, 0)  # Native geometric source gives c=1/2.
        native = _word_matrix(4, _native_scalar_word(3, signs))
        coefficient = .5
        odd = native[:8, 8:]
        np.testing.assert_allclose(native[:8, :8], .5 * np.eye(8),
                                   atol=ATOL, rtol=0)
        paulis = [np.kron(X, np.eye(2)), np.kron(np.eye(2), X),
                  np.kron(X, X)]
        actual, accepted, rejected, rotations = [], [], [], []
        for pauli, (low, high) in zip(paulis, ((0, 4), (4, 6), (6, 7))):
            rotation = _native_hopf_quarter_turn(low, high)
            expected = np.eye(8, dtype=complex)
            expected[low, low] = expected[high, high] = 0
            expected[high, low] = 1
            expected[low, high] = -1
            np.testing.assert_allclose(rotation, expected, atol=ATOL, rtol=0)
            projector = np.zeros((8, 8), dtype=complex)
            projector[low, low] = projector[high, high] = 1
            circuit, good, bad = _routed_stage(
                pauli, rotation, projector, coefficient, odd)
            np.testing.assert_allclose(circuit.conj().T @ circuit,
                                       np.eye(256), atol=ATOL, rtol=0)
            actual.append(circuit)
            accepted.append(good)
            rejected.append(bad)
            rotations.append(rotation)
        for left, right in ((0, 1), (1, 2)):
            self.assertGreater(np.linalg.norm(
                rotations[left] @ rotations[right]
                - rotations[right] @ rotations[left]), .5)
        embedding = np.eye(256, dtype=complex)[:, :64]
        two_columns = actual[1] @ actual[0] @ embedding
        two_block = embedding.conj().T @ two_columns
        ideal_two = np.kron(accepted[1] @ accepted[0], np.eye(8))
        np.testing.assert_allclose(two_block, ideal_two, atol=ATOL, rtol=0)
        three_block = embedding.conj().T @ actual[2] @ two_columns
        ideal_three = np.kron(
            accepted[2] @ accepted[1] @ accepted[0], np.eye(8))
        logical_mixed = rejected[2] @ rejected[1] @ rejected[0]
        cubic = odd @ odd @ odd
        mixed = np.kron(logical_mixed, cubic)
        np.testing.assert_allclose(three_block, ideal_three + mixed,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(logical_mixed[:, 0], np.eye(8)[:, 7],
                                   atol=ATOL, rtol=0)
        expected_norm = 3 * np.sqrt(3) / 8
        # All columns of the arbitrary dirty input, not a selected core state.
        witness = (three_block - ideal_three)[:, :8]
        np.testing.assert_allclose(
            witness.conj().T @ witness, expected_norm ** 2 * np.eye(8),
            atol=ATOL, rtol=0)
        conjugator = np.kron(np.eye(8), source)
        np.testing.assert_allclose(conjugator @ mixed @ conjugator, -mixed,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(
            (three_block - conjugator @ three_block @ conjugator) / 2,
            mixed, atol=ATOL, rtol=0)

    def test_different_masks_keep_odd_cubic_norm(self):
        source, _, _ = _operator_source(3)
        odds, coefficients = [], []
        for signs in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            scalar, masked = _scalar_block(source, signs)
            odds.append(scalar[:8, 8:])
            coefficients.append(np.trace(source @ masked).real / 8)
        product = odds[2] @ odds[1] @ odds[0]
        norm_squared = np.prod([1 - c * c for c in coefficients])
        np.testing.assert_allclose(product.conj().T @ product,
                                   norm_squared * np.eye(8), atol=ATOL, rtol=0)
        np.testing.assert_allclose(source @ product @ source, -product,
                                   atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.trace(product), 0, delta=ATOL)


if __name__ == '__main__':
    unittest.main()
