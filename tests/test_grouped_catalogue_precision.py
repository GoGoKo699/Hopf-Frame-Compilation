"""Bounded operator boundaries for grouped catalogue precision correction.

The two-bit atom uses an actual Clifford+T word on all nine-wire inputs.
The K=64, J=4 kernel fixture supplies a finite digit table and contracts
explicit SELECT permutations with uniform labels; it is not a table
emitter. The two-group check starts at the accepted-kernel boundary.
Nonwrapping history, padding and recovery remain covered by the shared
binary terminal tests. These checks certify no asymptotic resource bound
or production precision schedule.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.frames import real_frame_matrix

try:
    from .test_grouped_terminal_depth import _mcx
    from .test_operator_source_compiler import (
        _adjoint, _apply_native_word, _word_matrix,
    )
except ImportError:
    from test_grouped_terminal_depth import _mcx
    from test_operator_source_compiler import (
        _adjoint, _apply_native_word, _word_matrix,
    )


ATOL = 2e-10
D, K, J, SOURCE_SIZE = 4, 64, 4, 8
PI = 0.5
PHASES = (1, -1, 1j, -1j)
# Prefix rows differ genuinely. Only digits k>=2 occur, so the coefficients
# epsilon*2^-k obey the residual cap epsilon/4 in this bounded fixture.
DIGITS = (
    {(3, 0, 0): 2, (1, 2, 2): 3, (2, 3, 1): 2},
    {(0, 3, 1): 2, (2, 1, 3): 3, (1, 0, 0): 3},
    {(1, 3, 2): 2, (3, 2, 1): 3, (0, 1, 0): 2},
    {(2, 0, 3): 2, (0, 2, 0): 3, (3, 1, 1): 2},
)


def _atom_word():
    """X=0,1; u=2,3; v=4,5; rejection=6; H=7; dirty borrow=8."""
    xor_v = [('CX', 4, 0), ('CX', 5, 1)]
    negative = [('X', 0), ('X', 1)]
    mismatch = ([('CX', 7, 6)] + negative
                + _mcx([7, 0, 1], 6, 8) + _adjoint(negative))
    move = []
    for target, u, v in ((0, 2, 4), (1, 3, 5)):
        move += _mcx([7, u], target, 8)
        move += _mcx([7, v], target, 8)
    return xor_v + mismatch + _adjoint(xor_v) + move


def _atom_image(value):
    if not ((value >> 7) & 1):
        return value
    x, u, v = value & 3, (value >> 2) & 3, (value >> 4) & 3
    return value ^ (int(x != v) << 6) ^ u ^ v


def _select_image(source, target, under, atom, zero_digit, mode,
                  label, digit, prefix, *, active=True, enabled=True):
    """Complete permutation/phase action, retaining three separate failures."""
    if not active or mode:
        return source, target, under, atom, zero_digit, 1
    u, v, phase = label & 3, (label >> 2) & 3, label >> 4
    bit = enabled and DIGITS[prefix].get((u, v, phase)) == digit
    if not bit:
        return source, target, under, atom, zero_digit ^ 1, 1
    displacement = 2 * digit
    return ((source - displacement) % SOURCE_SIZE, target ^ u ^ v,
            under ^ int(source < displacement), atom ^ int(target != v),
            zero_digit, PHASES[phase])


def _select_accepted_kernel(prefix, *, enabled=True, merge_failures=False):
    """Zero-failure, zero-label block of B-dagger SELECT B.

    B places mode/lambda/k in a uniform superposition. Their unchanged
    labels contribute a second identical amplitude at B-dagger; this
    enumerates all branches of that finite block, including rejections.
    """
    size = SOURCE_SIZE * D
    result = np.zeros((size, size), complex)
    weight = 1 / (2 * K * J)
    for source in range(SOURCE_SIZE):
        for target in range(D):
            column = D * source + target
            for mode in range(2):
                for label in range(K):
                    for digit in range(J):
                        s, x, under, atom, zero, phase = _select_image(
                            source, target, 0, 0, 0, mode, label, digit,
                            prefix, enabled=enabled)
                        rejected = ((under ^ atom) or zero) if merge_failures else (
                            under or atom or zero)
                        if not rejected:
                            result[D * s + x, column] += weight * phase
    return result


def _shift(displacement):
    matrix = np.zeros((SOURCE_SIZE, SOURCE_SIZE), complex)
    for source in range(displacement, SOURCE_SIZE):
        matrix[source - displacement, source] = 1
    return matrix


def _dictionary_kernel(prefix, *, enabled=True):
    """Independent source-shift/matrix-unit formula for the accepted block."""
    result = (1 - PI) * np.eye(SOURCE_SIZE * D, dtype=complex)
    if enabled:
        for (u, v, phase), digit in DIGITS[prefix].items():
            atom = np.zeros((D, D), complex)
            atom[u, v] = PHASES[phase]
            result += PI / (K * J) * np.kron(_shift(2 * digit), atom)
    return result


def _commutator(a, b):
    return _word_matrix(2, _adjoint(b) + _adjoint(a) + b + a)


COARSE = (
    _commutator([('H', 0), ('CX', 0, 1)], [('T', 0)]),
    _commutator([('H', 1)], [('T', 1), ('CX', 1, 0)]),
    _commutator([('H', 0), ('H', 1)], [('T', 0), ('TDG', 1)]),
    _commutator([('CX', 1, 0), ('H', 0)], [('TDG', 0)]),
)


def _apply_group(columns, group, kernels, *, frozen_prefix=None):
    """Apply an accepted local kernel then its actual coarse group.

    Logical bits 0,1 are the first group; bits 2,3 are the second.
    The first group has the required lower-zero predicate. The later
    prefix is read afresh from the current logical state.
    """
    result = np.zeros_like(columns)
    for other in range(D):
        if group == 0:
            prefix, enabled = 0, other == 0
            indices = [16 * source + (other << 2) + target
                       for source in range(SOURCE_SIZE) for target in range(D)]
        else:
            prefix = other if frozen_prefix is None else frozen_prefix
            enabled = True
            indices = [16 * source + (target << 2) + other
                       for source in range(SOURCE_SIZE) for target in range(D)]
        kernel = kernels[prefix] if enabled else (1 - PI) * np.eye(SOURCE_SIZE * D)
        coarse = COARSE[prefix] if enabled else np.eye(D)
        result[indices] = np.kron(np.eye(SOURCE_SIZE), coarse) @ kernel @ columns[indices]
    return result


def _global_dictionary_stage(group):
    """Construct the global accepted stage directly from matrix units."""
    size = SOURCE_SIZE * 16
    contraction = (1 - PI) * np.eye(size, dtype=complex)
    coarse = np.eye(16, dtype=complex)
    for prefix in range(D) if group else (0,):
        for (u, v, phase), digit in DIGITS[prefix].items():
            logical_atom = np.zeros((16, 16), complex)
            row, column = ((u << 2) | prefix, (v << 2) | prefix) if group else (u, v)
            logical_atom[row, column] = PHASES[phase]
            contraction += PI / (K * J) * np.kron(_shift(2 * digit), logical_atom)
        indices = [(target << 2) | prefix for target in range(D)] if group else list(range(D))
        coarse[np.ix_(indices, indices)] = COARSE[prefix]
    return np.kron(np.eye(SOURCE_SIZE), coarse) @ contraction


class GroupedCataloguePrecisionTests(unittest.TestCase):
    def test_native_two_bit_atom_all_inputs_literal_phase_return_and_inverse(self):
        word = _atom_word()
        # Batches cover the entire operator, including arbitrary atom flag,
        # both labels, activity, and both possible dirty-helper inputs.
        for start in range(0, 512, 64):
            inputs = np.arange(start, start + 64)
            columns = np.zeros((512, 64), complex)
            columns[inputs, np.arange(64)] = 1
            expected = np.zeros_like(columns)
            expected[[_atom_image(int(value)) for value in inputs], np.arange(64)] = 1
            actual = _apply_native_word(9, word, columns)
            np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
            np.testing.assert_allclose(_apply_native_word(9, _adjoint(word), actual),
                                       columns, atol=ATOL, rtol=0)

        # Explicitly extract each zero-atom block: row u, column v, including
        # labels with u != v. No root-column-only comparison is used.
        for u in range(D):
            for v in range(D):
                base = (u << 2) | (v << 4) | (1 << 7)
                columns = np.zeros((512, D), complex)
                columns[[base | x for x in range(D)], np.arange(D)] = 1
                actual = _apply_native_word(9, word, columns)
                atom = actual[[base | x for x in range(D)]]
                expected = np.zeros((D, D), complex)
                expected[u, v] = 1
                np.testing.assert_allclose(atom, expected, atol=ATOL, rtol=0)
        witness = (1 << 2) | (1 << 7)  # x=v=0,u=1,H=1.
        column = np.eye(512, dtype=complex)[:, [witness]]
        wrong = _apply_native_word(9, word + word, column)
        self.assertGreater(np.linalg.norm(wrong - column), 1)

    def test_complex_actual_coarse_residual_has_64_literal_types(self):
        target = real_frame_matrix([.23, -.41, .67]).astype(complex)
        residual = COARSE[0].conj().T @ target - np.eye(D)
        self.assertGreater(np.linalg.norm(residual.imag), .1)
        coefficients = np.zeros((D, D, 4))
        coefficients[:, :, 0] = np.maximum(residual.real, 0)
        coefficients[:, :, 1] = np.maximum(-residual.real, 0)
        coefficients[:, :, 2] = np.maximum(residual.imag, 0)
        coefficients[:, :, 3] = np.maximum(-residual.imag, 0)
        self.assertEqual(coefficients.size, K)
        reconstructed = np.einsum('uvp,p->uv', coefficients, np.asarray(PHASES))
        np.testing.assert_allclose(reconstructed, residual, atol=ATOL, rtol=0)
        self.assertLessEqual(np.max(coefficients), np.linalg.norm(residual, 2) + ATOL)
        np.testing.assert_allclose(COARSE[0] @ (np.eye(D) + reconstructed),
                                   target, atol=ATOL, rtol=0)

    def test_64_type_selected_kernel_normalization_source_and_inactive_group(self):
        g = np.array([2 ** (-(j + 1) / 2) for j in range(SOURCE_SIZE - 1)]
                     + [2 ** (-(SOURCE_SIZE - 1) / 2)])
        embedding = np.kron(g[:, None], np.eye(D))
        source_error = max(np.linalg.norm(_shift(2 * k) @ g - 2 ** (-k) * g)
                           for k in range(J))
        self.assertLess(source_error, 2 * 2 ** (-J))
        epsilon = PI / ((1 - PI) * K * J)
        for prefix in range(D):
            actual = _select_accepted_kernel(prefix)
            expected = _dictionary_kernel(prefix)
            np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
            self.assertLessEqual(np.linalg.norm(actual, 2), 1 + ATOL)
            residual = np.zeros((D, D), complex)
            for (u, v, phase), digit in DIGITS[prefix].items():
                residual[u, v] += epsilon * 2 ** (-digit) * PHASES[phase]
            desired = embedding @ ((1 - PI) * (np.eye(D) + residual))
            self.assertLessEqual(np.linalg.norm(actual @ embedding - desired, 2),
                                 PI * source_error + ATOL)
        inactive = _select_accepted_kernel(0, enabled=False)
        np.testing.assert_allclose(inactive, (1 - PI) * np.eye(SOURCE_SIZE * D),
                                   atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(inactive - np.eye(SOURCE_SIZE * D)), 1)
        # Sharing the source/atom rejection bit would let two failures cancel.
        merged = _select_accepted_kernel(0, merge_failures=True)
        self.assertGreater(np.linalg.norm(merged - _dictionary_kernel(0)), .001)
        # H=0 is stronger than G=0: even occupied flags and labels are fixed.
        for label in range(K):
            for digit in range(J):
                got = _select_image(7, 3, 1, 1, 1, 0, label, digit, 2, active=False)
                self.assertEqual(got, (7, 3, 1, 1, 1, 1))

    def test_two_groups_read_changed_prefix_and_keep_the_shared_source(self):
        kernels = [_select_accepted_kernel(prefix) for prefix in range(D)]
        first, second = [_global_dictionary_stage(group) for group in range(2)]
        identity = np.eye(SOURCE_SIZE * 16, dtype=complex)
        actual = _apply_group(_apply_group(identity, 0, kernels), 1, kernels)
        np.testing.assert_allclose(actual, second @ first, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(first @ second - second @ first), .1)

        # Freeze the first prefix separately for each original input column.
        # Earlier coarse gates change it coherently, so this reuse is invalid.
        after_first = _apply_group(identity, 0, kernels)
        wrong = np.zeros_like(identity)
        for prefix in range(D):
            columns = [column for column in range(len(identity)) if column & 3 == prefix]
            wrong[:, columns] = _apply_group(after_first[:, columns], 1, kernels,
                                             frozen_prefix=prefix)
        self.assertGreater(np.linalg.norm(wrong - actual), .1)

        g = np.array([2 ** (-(j + 1) / 2) for j in range(SOURCE_SIZE - 1)]
                     + [2 ** (-(SOURCE_SIZE - 1) / 2)])
        embedding = np.kron(g[:, None], np.eye(16))
        final = embedding.conj().T @ actual @ embedding
        premature = ((embedding.conj().T @ second @ embedding)
                     @ (embedding.conj().T @ first @ embedding))
        self.assertGreater(np.linalg.norm(final - premature), 1e-9)


if __name__ == '__main__':
    unittest.main()
