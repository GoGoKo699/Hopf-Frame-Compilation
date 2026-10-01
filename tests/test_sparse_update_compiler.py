"""Bounded algebra audits for packing ancestor-closed Hopf updates.

Tree matrices have dimension at most 16. The dense residual fixture applies
small ideal scalar dilations to batched columns, without constructing its
256-dimensional circuit matrix. It checks coefficient placement, distinct
failure flags, literal phases, inverse order, and the complete amplified
isometry. It is not a native high-precision emitter or a gate-count proof.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_conditional_suffix_compiler import H, _apply_axes
from tests.test_operator_source_compiler import _word_matrix
from tests.test_tree_residual_structure import _addressed_frame, _marker, _su2


ATOL = 5e-11
PHASES = (1, -1, 1j, -1j)
X = np.array([[0, 1], [1, 0]], dtype=complex)


def _ancestor_closure(nodes):
    result = set()
    for node in nodes:
        while node:
            result.add(node)
            node //= 2
    return result


def _masked_frame(height, words, nodes):
    return _addressed_frame(height, {
        node: word if node in nodes else np.eye(2)
        for node, word in words.items()
    })


def _native_words(height):
    programs = ([('H', 0), ('T', 0), ('H', 0)],
                [('T', 0), ('H', 0), ('S', 0)],
                [('H', 0), ('TDG', 0), ('H', 0), ('T', 0)])
    result = {}
    for node in range(1, 1 << height):
        program = _word_matrix(1, programs[node % len(programs)])
        echo = X @ program.conj().T @ X @ program
        result[node] = echo @ echo
    return result


def _affine_transposition(height, first, second):
    """X/CNOT affine conjugation of one all-zero-controlled target flip.

    The returned ledger lists literal affine X/CNOT gates and the one
    multi-controlled X. Its native decomposition is deliberately not priced
    by this finite permutation test.
    """
    assert first != second
    difference = first ^ second
    pivot = (difference & -difference).bit_length() - 1
    targets = [bit for bit in range(height)
               if bit != pivot and (difference >> bit) & 1]
    images = []
    for original in range(1 << height):
        value = original ^ first
        for bit in targets:
            value ^= ((value >> pivot) & 1) << bit
        if value & ~(1 << pivot) == 0:
            value ^= 1 << pivot
        for bit in reversed(targets):
            value ^= ((value >> pivot) & 1) << bit
        images.append(value ^ first)
    ledger = (2 * first.bit_count(), 2 * len(targets), height - 1)
    return np.array(images), ledger


def _pack_modes(height, support):
    images = np.arange(1 << height)
    transpositions = []
    for destination, original in enumerate(support):
        current = int(images[original])
        if current != destination:
            swap, ledger = _affine_transposition(height, current, destination)
            images = swap[images]
            transpositions.append((current, destination, ledger))
    return images, transpositions


def _column_atom(size, column, reverse_order=False):
    """Accepted block |+><column|: predicate, XOR, then Hadamards."""
    assert size == 2
    predicate = np.zeros((2 * size, 2 * size))
    for flag in range(2):
        for data in range(size):
            output_flag = flag ^ (data != column)
            predicate[output_flag * size + data, flag * size + data] = 1
    displacement = np.eye(size)[:, np.arange(size) ^ column]
    spread = np.kron(np.eye(2), H @ displacement)
    return predicate @ spread if reverse_order else spread @ predicate


def _scalar_dilation(coefficient, perturbation=0):
    """Ideal accepted c I + epsilon X on an arbitrary dirty spectator."""
    eigenvalues = np.array([coefficient + perturbation,
                            coefficient - perturbation])
    assert np.max(np.abs(eigenvalues)) < 1
    roots = np.sqrt(1 - eigenvalues ** 2)
    accepted = coefficient * np.eye(2) + perturbation * X
    rejected = roots.mean() * np.eye(2) + (roots[0] - roots[1]) / 2 * X
    return np.block([[accepted, -rejected], [rejected, accepted]])


class _DenseResidualBlock:
    """Mode, (column, phase), scalar flag, atom flag, data, dirty core.

    The label register has K=4M values. The only dense circuit components
    are 2x2, 4x4 and 8x8; every dirty input column is retained in the final
    isometry. Scalar perturbations can rotate the dirty input coherently.
    """

    def __init__(self, target, perturbation=0):
        self.target = target
        self.logical_size = 2
        self.count = 8
        self.data_size = 4
        self.shape = (2, self.count, 2, 2, 2, 2)
        self.size = int(np.prod(self.shape))
        self.label_h = np.kron(np.kron(H, H), H)
        residual = target - np.eye(2)
        components = np.stack((np.maximum(residual.real, 0),
                               np.maximum(-residual.real, 0),
                               np.maximum(residual.imag, 0),
                               np.maximum(-residual.imag, 0)), axis=-1)
        self.coefficients = self.count * np.sqrt(2) * components
        self.scalars = {}
        for column in range(2):
            for phase in range(4):
                for row in range(2):
                    epsilon = perturbation * (row + 1) if column == phase == 0 else 0
                    self.scalars[row, column, phase] = _scalar_dilation(
                        self.coefficients[row, column, phase], epsilon)

    def apply(self, columns, inverse=False):
        state = columns.reshape(self.shape + (-1,)).copy()
        state = _apply_axes(state, H, (0,))
        selected = _apply_axes(state[1], self.label_h, (0,))
        for column in range(2):
            atom = _column_atom(2, column)
            for phase, omega in enumerate(PHASES):
                index = 4 * column + phase
                term = selected[index]
                if not inverse:
                    term = _apply_axes(term, atom, (1, 2))
                for row in range(2):
                    scalar = self.scalars[row, column, phase]
                    term[:, :, row] = _apply_axes(
                        term[:, :, row], scalar.conj().T if inverse else scalar,
                        (0, 2))
                if inverse:
                    term = _apply_axes(term, atom.conj().T, (1, 2))
                selected[index] = (np.conj(omega) if inverse else omega) * term
        state[1] = _apply_axes(selected, self.label_h, (0,))
        return _apply_axes(state, H, (0,)).reshape(self.size, -1)

    def embedding(self):
        result = np.zeros((self.size, self.data_size), dtype=complex)
        result[:self.data_size] = np.eye(self.data_size)
        return result

    def amplify(self, columns):
        result = self.apply(columns)
        result[:self.data_size] *= -1
        result = self.apply(result, inverse=True)
        result[:self.data_size] *= -1
        return -self.apply(result)

    def controlled_amplify(self, columns):
        """Both reflections may be shared; the final minus is controlled."""
        result = columns.copy()
        result[1] = self.apply(result[1])
        result[:, :self.data_size] *= -1
        result[1] = self.apply(result[1], inverse=True)
        result[:, :self.data_size] *= -1
        result[1] = -self.apply(result[1])
        return result


class SparseUpdateCompilerTests(unittest.TestCase):
    def test_ancestor_closure_factorization_and_exact_complement(self):
        cases = ((1, (1,)), (2, (2,)), (3, (3, 4)),
                 (4, (5, 12)), (4, (9, 14)), (4, (1,)), (4, ()))
        for height, changed in cases:
            with self.subTest(height=height, changed=changed):
                size = 1 << height
                coarse = _native_words(height)
                target = {node: word.copy() for node, word in coarse.items()}
                for node in changed:
                    target[node] = _su2(0.07 * node, -0.03 * node, 0.05)
                closure = _ancestor_closure(changed)
                outside = set(coarse) - closure
                forest = _masked_frame(height, coarse, outside)
                q_target = _masked_frame(height, target, closure)
                q_coarse = _masked_frame(height, coarse, closure)
                w = _addressed_frame(height, target)
                c = _addressed_frame(height, coarse)
                np.testing.assert_allclose(w, forest @ q_target, atol=ATOL, rtol=0)
                np.testing.assert_allclose(c, forest @ q_coarse, atol=ATOL, rtol=0)
                residual = c.conj().T @ w
                np.testing.assert_allclose(residual, q_coarse.conj().T @ q_target,
                                           atol=ATOL, rtol=0)
                support = sorted({0} | {_marker(node, height) for node in closure})
                self.assertEqual(len(support), len(closure) + 1)
                complement = sorted(set(range(size)) - set(support))
                for unitary in (q_target, q_coarse, residual):
                    np.testing.assert_allclose(unitary[:, complement],
                                               np.eye(size)[:, complement],
                                               atol=ATOL, rtol=0)
                images, swaps = _pack_modes(height, support)
                self.assertLessEqual(len(swaps), len(support))
                np.testing.assert_array_equal(images[support], np.arange(len(support)))
                p = np.eye(size)[:, images]
                packed = p @ residual @ p.conj().T
                width = (len(support) - 1).bit_length()
                self.assertLessEqual(len(support), 1 << width)
                np.testing.assert_allclose(packed[:, len(support):],
                                           np.eye(size)[:, len(support):],
                                           atol=ATOL, rtol=0)

    def test_affine_basis_transpositions_on_every_small_input(self):
        for height in range(1, 5):
            size = 1 << height
            for first in range(size):
                for second in range(first + 1, size):
                    images, (x_count, cx_count, controls) = _affine_transposition(
                        height, first, second)
                    expected = np.arange(size)
                    expected[first], expected[second] = second, first
                    np.testing.assert_array_equal(images, expected)
                    self.assertLessEqual(x_count, 2 * height)
                    self.assertLessEqual(cx_count, 2 * (height - 1))
                    self.assertEqual(controls, height - 1)

    def test_dense_labels_reconstruct_complex_rows_after_column_atom(self):
        target = _su2(0.002, 0.001, -0.0003)
        block = _DenseResidualBlock(target)
        self.assertLess(np.max(block.coefficients), 1 / 4)
        reconstructed = np.zeros((2, 2), dtype=complex)
        for column in range(2):
            atom = _column_atom(2, column)
            expected_atom = np.zeros((2, 2))
            expected_atom[:, column] = 1 / np.sqrt(2)
            np.testing.assert_allclose(atom[:2, :2], expected_atom, atol=ATOL, rtol=0)
            for phase, omega in enumerate(PHASES):
                reconstructed += omega * np.diag(block.coefficients[:, column, phase]) \
                    @ atom[:2, :2] / block.count
        np.testing.assert_allclose(reconstructed, target - np.eye(2), atol=ATOL, rtol=0)
        embedded = block.embedding()
        actual = block.apply(embedded)
        expected = np.kron(target, np.eye(2)) / 2
        np.testing.assert_allclose(actual[:block.data_size], expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(block.apply(actual, inverse=True), embedded,
                                   atol=ATOL, rtol=0)

    def test_atom_order_and_separate_scalar_flag_are_necessary(self):
        coefficient = np.array([0.04, 0.11])
        atom = _column_atom(2, 1)
        expected = np.diag(coefficient) @ atom[:2, :2]
        reversed_atom = _column_atom(2, 1, reverse_order=True)
        self.assertGreater(np.linalg.norm(np.diag(coefficient)
                                          @ reversed_atom[:2, :2] - expected), 0.05)
        # A single shared flag permits atom rejection to re-enter acceptance.
        shared_scalar = np.zeros((4, 4))
        for data, value in enumerate(coefficient):
            rows = [data, 2 + data]
            root = np.sqrt(1 - value ** 2)
            shared_scalar[np.ix_(rows, rows)] = [[value, -root], [root, value]]
        shared = shared_scalar @ atom
        self.assertGreater(np.linalg.norm(shared[:2, :2] - expected, ord=2), 0.9)
        # Filtering before the atom reads the input column, not its output row.
        self.assertGreater(np.linalg.norm(atom[:2, :2] @ np.diag(coefficient)
                                          - expected, ord=2), 0.04)

    def test_actual_inverse_full_amplification_and_coherent_inactive_branch(self):
        block = _DenseResidualBlock(_su2(0.002, 0.001, -0.0003))
        embedding = block.embedding()
        target = np.kron(block.target, np.eye(2))
        np.testing.assert_allclose(block.amplify(embedding), embedding @ target,
                                   atol=ATOL, rtol=0)
        rng = np.random.default_rng(171)
        columns = rng.normal(size=(block.size, 3)) + 1j * rng.normal(size=(block.size, 3))
        np.testing.assert_allclose(block.apply(block.apply(columns), inverse=True),
                                   columns, atol=ATOL, rtol=0)
        # Check every inactive scratch/data/core basis input in batches.
        for start in range(0, block.size, 8):
            inactive = np.zeros((2, block.size, 8), dtype=complex)
            inactive[0, start + np.arange(8), np.arange(8)] = 1
            np.testing.assert_array_equal(block.controlled_amplify(inactive), inactive)
        coherent = np.zeros((2, block.size, block.data_size), dtype=complex)
        coherent[0] = coherent[1] = embedding / np.sqrt(2)
        expected = coherent.copy()
        expected[1] = embedding @ target / np.sqrt(2)
        np.testing.assert_allclose(block.controlled_amplify(coherent), expected,
                                   atol=ATOL, rtol=0)
        # An unconditional final minus would put a relative phase on h=0.
        wrong = expected.copy()
        wrong[0] *= -1
        self.assertGreater(np.linalg.norm(wrong - expected), 1)

    def test_scalar_dirty_perturbation_is_counted_with_all_rejection(self):
        block = _DenseResidualBlock(_su2(0.002, 0.001, -0.0003), perturbation=0.001)
        embedding = block.embedding()
        target = np.kron(block.target, np.eye(2))
        accepted = block.apply(embedding)[:block.data_size]
        amplified = block.amplify(embedding)
        np.testing.assert_allclose(amplified[:block.data_size],
                                   3 * accepted - 4 * accepted @ accepted.conj().T @ accepted,
                                   atol=ATOL, rtol=0)
        zeta = np.linalg.norm(2 * accepted - target, ord=2)
        self.assertGreater(zeta, 1e-5)
        full_error = np.linalg.norm(amplified - embedding @ target, ord=2)
        self.assertLessEqual(full_error, 3 * zeta + zeta ** 2 + ATOL)
        self.assertGreater(np.linalg.norm(amplified[block.data_size:]), 1e-5)
        self.assertGreater(np.linalg.norm((2 * accepted - target)[0::2, 1::2]), 1e-5)


if __name__ == '__main__':
    unittest.main()
