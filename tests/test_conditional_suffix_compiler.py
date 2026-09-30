"""Finite audits of the conditional-suffix grouping interface.

The scalar source below is the explicit native word from the operator-source
tests. Vectorized application keeps every dirty input column while avoiding
large dense circuit matrices. These fixtures check block composition,
failure flags, coherent inactive sectors, and sparse support; they do not
prove asymptotic gate counts or a complete synthesis theorem.
"""
from __future__ import annotations

import unittest
import math

import numpy as np

try:
    from .test_operator_source_compiler import (
        H, _operator_source, _scalar_block, _sign_word, _word_matrix,
    )
except ImportError:
    from test_operator_source_compiler import (
        H, _operator_source, _scalar_block, _sign_word, _word_matrix,
    )


ATOL = 4e-11


def _apply_axes(state, matrix, axes):
    """Apply an ordinary matrix on named tensor axes, preserving batch axes."""
    axes = tuple(axes)
    order = axes + tuple(axis for axis in range(state.ndim) if axis not in axes)
    arranged = np.transpose(state, order)
    shape = arranged.shape
    result = matrix @ arranged.reshape(matrix.shape[1], -1)
    return np.transpose(result.reshape(shape), np.argsort(order))


def _matrix_unit_dilation(u, v):
    """Flag,target ordering; predicate first, then the endpoint XOR."""
    matrix = np.zeros((4, 4), dtype=complex)
    for flag in range(2):
        for target in range(2):
            new_flag = flag ^ (target != v)
            new_target = target ^ (u ^ v)
            matrix[2 * new_flag + new_target, 2 * flag + target] = 1
    return matrix


class _SmallResidualBlock:
    """Q encodes (I+E)/2 using a mode, four labels, and separate flags.

    Axes are mode, label, scalar flag, matrix-unit flag, local target, dirty
    core, and an arbitrary column batch. The five first scratch bits can
    model a conditionally initialized logical suffix. No state is prepared
    on the three-qubit dirty core.
    """

    def __init__(self, atoms):
        assert len(atoms) == 4
        self.atoms = atoms
        self.core_size = 8
        self.data_size = 2 * self.core_size
        self.scratch_size = 2 * 4 * 2 * 2
        self.size = self.scratch_size * self.data_size
        self.shape = (2, 4, 2, 2, 2, self.core_size)
        source = _operator_source(3)[0]
        self.scalars = []
        for coefficient, _, _, _ in atoms:
            k = int(round(2 * (1 - coefficient)))
            assert abs(coefficient - (1 - k / 2)) < ATOL
            self.scalars.append(_scalar_block(source, _sign_word(k, 2))[0])
        self.label_h = np.kron(H, H)

    def apply(self, columns, inverse=False):
        state = columns.reshape(self.shape + (-1,)).copy()
        state = _apply_axes(state, H, (0,))
        selected = _apply_axes(state[1], self.label_h, (0,))
        for label, (_, phase, u, v) in enumerate(self.atoms):
            term = selected[label]
            scalar = self.scalars[label]
            matrix_unit = _matrix_unit_dilation(u, v)
            if inverse:
                term = np.conj(phase) * _apply_axes(
                    term, matrix_unit.conj().T, (1, 2))
                term = _apply_axes(term, scalar.conj().T, (0, 3))
            else:
                term = _apply_axes(term, scalar, (0, 3))
                term = phase * _apply_axes(term, matrix_unit, (1, 2))
            selected[label] = term
        state[1] = _apply_axes(selected, self.label_h, (0,))
        state = _apply_axes(state, H, (0,))
        return state.reshape(self.size, -1)

    def embedding(self):
        result = np.zeros((self.size, self.data_size), dtype=complex)
        result[:self.data_size] = np.eye(self.data_size)
        return result

    def amplify(self, columns):
        # -Q R Q^dagger R Q with R = I - 2 J J^dagger.
        result = self.apply(columns)
        result[:self.data_size] *= -1
        result = self.apply(result, inverse=True)
        result[:self.data_size] *= -1
        return -self.apply(result)

    def expected_residual(self):
        residual = np.zeros((2, 2), dtype=complex)
        for coefficient, phase, u, v in self.atoms:
            residual[u, v] += coefficient * phase / 4
        return np.kron(residual, np.eye(self.core_size))


def _small_frame(height, word_factory):
    """Complete binary frame in addressed shallow-to-deep gate order."""
    size = 1 << height
    frame = np.eye(size, dtype=complex)
    for depth in range(height):
        stride = 1 << (height - depth - 1)
        for prefix in range(1 << depth):
            rows = [prefix * (stride << 1), prefix * (stride << 1) + stride]
            frame[rows] = word_factory(depth, prefix) @ frame[rows]
    return frame


def _column_support(marker, height):
    if marker == 0:
        return set(range(1 << height))
    lowbit = marker & -marker
    subtree_size = 2 * lowbit
    start = marker - marker % subtree_size
    return set(range(start, start + subtree_size))


def _group_ledger(n, ratio=64, base=256):
    """An illustrative fixed partition, not a synthesis-workspace certificate."""
    r = min(n, base)
    groups = []
    while r < n:
        s = min(n - r, r // ratio)
        # Exact ceil(log2(4*((2*s-1)*2**s+2))), without a large power of two.
        label_bits = 4 if s == 1 else s + (2 * s - 2).bit_length() + 2
        groups.append((r, s, label_bits))
        r += s
    return groups


class ConditionalSuffixCompilerTests(unittest.TestCase):
    def test_controlled_t_dirty_phase_echo_and_helper_return(self):
        # Qubit 0 stores the preserved predicate F; 1 and 2 are arbitrary
        # dirty helpers z,w. CCZ is expanded as H-CCX-H, so this is the
        # chronological X/CX/CCX/H/T/S word, including its literal phase.
        word = [
            ('CX', 0, 1), ('T', 1), ('CX', 0, 1), ('TDG', 1),
            ('CCX', 0, 1, 2), ('S', 2), ('CCX', 0, 1, 2), ('SDG', 2),
            ('H', 2), ('CCX', 0, 1, 2), ('H', 2),
        ]
        actual = _word_matrix(3, word)
        omega = np.exp(1j * np.pi / 4)
        expected = np.diag([omega ** (basis & 1) for basis in range(8)])
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        for predicate in (0, 1):
            for z in (0, 1):
                for w in (0, 1):
                    phase = (omega ** (predicate - 2 * z * predicate)
                             * 1j ** (z * predicate - 2 * w * z * predicate)
                             * (-1) ** (w * z * predicate))
                    self.assertAlmostEqual(phase, omega ** predicate, delta=ATOL)
        # The full diagonal equality checks every helper input, and hence
        # coherent helper superpositions and external entanglement as well.

    def test_matrix_unit_predicate_and_actual_inverse(self):
        for u in range(2):
            for v in range(2):
                dilation = _matrix_unit_dilation(u, v)
                expected = np.zeros((2, 2))
                expected[u, v] = 1
                np.testing.assert_allclose(dilation[:2, :2], expected, atol=0)
                np.testing.assert_allclose(
                    dilation.conj().T @ dilation, np.eye(4), atol=0)

    def test_native_scalar_and_separate_failure_flags(self):
        source = _operator_source(3)[0]
        scalar = _scalar_block(source, _sign_word(1, 2))[0]
        core_size = source.shape[0]
        data_size = 2 * core_size
        embedding = np.zeros((4 * data_size, data_size), dtype=complex)
        embedding[:data_size] = np.eye(data_size)
        state = embedding.reshape(2, 2, 2, core_size, data_size)
        state = _apply_axes(state, scalar, (0, 3))
        state = _apply_axes(state, _matrix_unit_dilation(1, 0), (1, 2))
        atom = np.array([[0, 0], [1, 0]])
        expected = 0.5 * np.kron(atom, np.eye(core_size))
        np.testing.assert_allclose(state[0, 0].reshape(data_size, data_size),
                                   expected, atol=ATOL, rtol=0)

        # Reusing one failure flag allows a rejected scalar path to return.
        shared = np.zeros((2 * data_size, data_size), dtype=complex)
        shared[:data_size] = np.eye(data_size)
        shared = shared.reshape(2, 2, core_size, data_size)
        shared = _apply_axes(shared, scalar, (0, 2))
        shared = _apply_axes(shared, _matrix_unit_dilation(1, 0), (0, 1))
        discrepancy = shared[0].reshape(data_size, data_size) - expected
        self.assertGreater(np.linalg.norm(discrepancy, ord=2), 0.5)

    def test_padded_lcu_preserves_literal_complex_phases(self):
        atoms = [(0.5, 1, 0, 0), (1, -1, 0, 1),
                 (0.5, 1j, 1, 0), (0, -1j, 1, 1)]
        block = _SmallResidualBlock(atoms)
        embedding = block.embedding()
        actual = block.apply(embedding)
        expected = (np.eye(block.data_size) + block.expected_residual()) / 2
        np.testing.assert_allclose(actual[:block.data_size], expected,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(block.apply(actual, inverse=True), embedding,
                                   atol=ATOL, rtol=0)

    def test_oblivious_amplification_includes_every_dirty_input(self):
        atoms = [(0, 1, 0, 0), (0.5, -1, 0, 1),
                 (0.5, 1, 1, 0), (0, 1, 1, 1)]
        block = _SmallResidualBlock(atoms)
        embedding = block.embedding()
        residual = block.expected_residual()
        approximate = np.eye(block.data_size) + residual
        target = approximate / np.sqrt(1 + (1 / 8) ** 2)
        np.testing.assert_allclose(target.conj().T @ target,
                                   np.eye(block.data_size), atol=ATOL, rtol=0)
        b = block.apply(embedding)[:block.data_size]
        amplified = block.amplify(embedding)
        np.testing.assert_allclose(amplified[:block.data_size],
                                   3 * b - 4 * b @ b.conj().T @ b,
                                   atol=ATOL, rtol=0)
        zeta = np.linalg.norm(2 * b - target, ord=2)
        full_error = np.linalg.norm(amplified - embedding @ target, ord=2)
        self.assertLessEqual(full_error, 4 * zeta + ATOL)
        self.assertGreater(np.linalg.norm(amplified[block.data_size:]), 1e-3)

    def test_conditional_suffix_uncompute_keeps_leakage_in_full_error(self):
        atoms = [(0, 1, 0, 0), (0.5, -1, 0, 1),
                 (0.5, 1, 1, 0), (0, 1, 1, 1)]
        block = _SmallResidualBlock(atoms)
        data = block.data_size
        # The external flag h starts zero; every suffix/core input column
        # is included. Compute h=[suffix=0], control A on h, then uncompute.
        columns = np.zeros((2, block.size, block.size), dtype=complex)
        columns[0] = np.eye(block.size)
        columns[:, :data] = columns[::-1, :data]
        columns[1] = block.amplify(columns[1])
        columns[:, :data] = columns[::-1, :data]
        np.testing.assert_allclose(columns[0, :, data:],
                                   np.eye(block.size)[:, data:], atol=0, rtol=0)
        np.testing.assert_allclose(columns[1, :, data:], 0, atol=0, rtol=0)

        target = (np.eye(data) + block.expected_residual()) / np.sqrt(65 / 64)
        ideal = np.zeros_like(columns)
        ideal[0] = np.eye(block.size)
        ideal[0, :data, :data] = target
        difference = (columns - ideal).reshape(2 * block.size, block.size)
        # All nonzero error columns are active; their spectral norm covers
        # arbitrary superpositions and reference-entangled dirty inputs.
        full_error = np.linalg.norm(difference[:, :data], ord=2)
        local_error = np.linalg.norm(
            block.amplify(block.embedding()) - block.embedding() @ target, ord=2)
        self.assertAlmostEqual(full_error, local_error, delta=ATOL)
        self.assertGreater(np.linalg.norm(columns[1]), 1e-3)

    def test_sparse_group_support_with_exact_complex_coarse_words(self):
        native_words = [
            [('H', 0), ('T', 0), ('H', 0)],
            [('T', 0), ('H', 0), ('S', 0)],
            [('H', 0), ('TDG', 0), ('H', 0), ('T', 0)],
        ]
        def target_word(depth, prefix):
            angle = 0.19 + 0.13 * depth + 0.07 * prefix
            return np.array([[np.cos(angle), -np.sin(angle)],
                             [np.sin(angle), np.cos(angle)]])

        for height in (1, 2, 3):
            size = 1 << height
            target = _small_frame(height, target_word)
            coarse = _small_frame(
                height, lambda depth, prefix: _word_matrix(
                    1, native_words[(depth + prefix) % len(native_words)]))
            support = [_column_support(marker, height) for marker in range(size)]
            permitted = np.array([[bool(left & right) for right in support]
                                  for left in support])
            for frame in (target, coarse):
                np.testing.assert_allclose(frame.conj().T @ frame, np.eye(size),
                                           atol=ATOL, rtol=0)
                for marker in range(size):
                    outside = sorted(set(range(size)) - support[marker])
                    np.testing.assert_allclose(frame[outside, marker], 0,
                                               atol=ATOL, rtol=0)
            residual = coarse.conj().T @ target - np.eye(size)
            np.testing.assert_allclose(residual[~permitted], 0, atol=ATOL, rtol=0)
            self.assertEqual(int(permitted.sum()), (2 * height - 1) * size + 2)
            self.assertGreater(np.linalg.norm(coarse.imag), 0.1)

    def test_geometric_group_integer_dirty_and_error_ledgers(self):
        # The proof chooses C sufficiently large for the exact one-qubit
        # synthesis primitive. This finite test isolates the integer dirty,
        # precision, and lookup ledger from that unspecified fixed constant.
        ratio, base = 64, 256
        for n in (1, 2, 64, 255, 256, 257, 511, 1024, 10_000, 1_000_000):
            groups = _group_ledger(n, ratio, base)
            self.assertEqual(sum(s for _, s, _ in groups) + min(n, base), n)
            self.assertLessEqual(len(groups), 2 * ratio * n.bit_length())
            self.assertLessEqual(sum(r for r, _, _ in groups), 2 * ratio * n)
            self.assertLessEqual(sum(r // 4 for r, _, _ in groups), ratio * n)
            normalized_lookup = 0.0
            normalized_error = 0.0
            previous_end = min(n, base)
            for r, s, label_bits in groups:
                self.assertEqual(r, previous_end)
                self.assertLessEqual(s, r // ratio)
                previous_end = r + s
                if s <= 100:
                    padded_entries = 4 * ((2 * s - 1) * (1 << s) + 2)
                    self.assertEqual(label_bits, (padded_entries - 1).bit_length())
                address_bits = n - r - s + label_bits + 2
                for precision in (6, 31, 1000):
                    source_width = precision + r // 4 + 8
                    self.assertLessEqual(source_width + address_bits + 8,
                                         precision + n + 7)
                # Q/N = 2**(label_bits-s-r), with the unchanged prefix used
                # directly as an address. There is no extra copied prefix.
                normalized_lookup += math.ldexp(1.0, label_bits - s - r)
                normalized_error += 10 * math.ldexp(1.0, -(r // 4) - 8)
            self.assertLessEqual(normalized_lookup, 32 * (base + 1) * 2.0 ** -base)
            self.assertLessEqual(normalized_error, 80 * 2.0 ** (-8 - base // 4))
            baseline_error = 5 * np.sqrt(2) / 8 * (1 - 2.0 ** -min(n, base))
            self.assertLess(baseline_error + normalized_error, 1)


if __name__ == '__main__':
    unittest.main()
