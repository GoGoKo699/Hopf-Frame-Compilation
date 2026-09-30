"""Finite checks of the nilpotent-source obstruction and its boundaries.

These deterministic, small matrix fixtures check the stated formulas and
assumption boundaries. They do not prove the dimension-independent lemma or
a lower bound for unrestricted Hopf-frame compilation.
"""
from __future__ import annotations

import unittest

import numpy as np


ATOL = 2e-12


def _bound(alpha, dimension_ratio):
    return alpha ** dimension_ratio / (1 + alpha) ** (dimension_ratio - 1)


def _shift(dimension):
    """S|j> = |j-1> for j > 0, and S|0> = 0."""
    return np.diag(np.ones(dimension - 1), 1).astype(complex)


def _geometric_embedding(dimension, dirty_dimension, alpha):
    normalization = np.sqrt((1 - alpha ** 2) / (1 - alpha ** (2 * dimension)))
    vector = normalization * alpha ** np.arange(dimension)
    return np.kron(vector[:, None], np.eye(dirty_dimension)), normalization


def _random_unitary(rng, dimension):
    matrix = rng.normal(size=(dimension, dimension))
    matrix = matrix + 1j * rng.normal(size=matrix.shape)
    return np.linalg.qr(matrix)[0]


def _pauli(width, factors):
    """Qubit zero is the least significant tensor factor."""
    result = np.ones((1, 1), dtype=complex)
    for qubit in reversed(range(width)):
        result = np.kron(result, factors.get(qubit, np.eye(2)))
    return result


class SourceReuseLimitsTests(unittest.TestCase):
    def test_smallest_singular_subspaces_of_nilpotent_contractions(self):
        rng = np.random.default_rng(20260930)
        for ratio in (1, 2, 4):
            for input_dimension in (1, 2, 3):
                dimension = ratio * input_dimension
                matrix = rng.normal(size=(dimension, dimension))
                matrix = matrix + 1j * rng.normal(size=matrix.shape)
                contraction = np.triu(matrix, 1)
                norm = np.linalg.norm(contraction, ord=2)
                if norm:
                    contraction /= norm
                np.testing.assert_allclose(
                    np.linalg.matrix_power(contraction, dimension), 0,
                    atol=ATOL, rtol=0)
                self.assertLessEqual(np.linalg.norm(contraction, ord=2), 1 + ATOL)

                for alpha in (0.25, 0.5, 0.75):
                    with self.subTest(ratio=ratio, input_dimension=input_dimension,
                                      alpha=alpha):
                        difference = alpha * np.eye(dimension) - contraction
                        _, singular_values, vh = np.linalg.svd(difference)
                        # This D-dimensional subspace minimizes the residual
                        # operator norm over all encoding isometries of size D.
                        embedding = vh.conj().T[:, -input_dimension:]
                        np.testing.assert_allclose(
                            embedding.conj().T @ embedding,
                            np.eye(input_dimension), atol=ATOL, rtol=0)
                        residual = np.linalg.norm(difference @ embedding, ord=2)
                        self.assertAlmostEqual(
                            residual, singular_values[-input_dimension], delta=ATOL)
                        self.assertGreaterEqual(residual + ATOL, _bound(alpha, ratio))
                        self.assertAlmostEqual(
                            np.linalg.slogdet(difference)[1],
                            dimension * np.log(alpha), delta=ATOL)

    def test_geometric_source_and_dirty_tensor_factor(self):
        for ratio in (1, 2, 4, 8):
            for alpha in (0.25, 0.5, 1 / np.sqrt(2)):
                expected_errors = []
                for dirty_dimension in (1, 2, 3):
                    with self.subTest(ratio=ratio, alpha=alpha,
                                      dirty_dimension=dirty_dimension):
                        contraction = np.kron(_shift(ratio), np.eye(dirty_dimension))
                        embedding, normalization = _geometric_embedding(
                            ratio, dirty_dimension, alpha)
                        np.testing.assert_allclose(
                            embedding.conj().T @ embedding,
                            np.eye(dirty_dimension), atol=ATOL, rtol=0)
                        residual_matrix = contraction @ embedding - alpha * embedding
                        expected_matrix = np.zeros_like(embedding)
                        expected_matrix[-dirty_dimension:] = (
                            -normalization * alpha ** ratio * np.eye(dirty_dimension))
                        np.testing.assert_allclose(
                            residual_matrix, expected_matrix, atol=ATOL, rtol=0)
                        residual = np.linalg.norm(residual_matrix, ord=2)
                        self.assertAlmostEqual(
                            residual, normalization * alpha ** ratio, delta=ATOL)
                        self.assertGreaterEqual(residual + ATOL, _bound(alpha, ratio))
                        expected_errors.append(residual)
                np.testing.assert_allclose(
                    expected_errors, expected_errors[0], atol=ATOL, rtol=0)

    def test_general_encoding_and_basis_change_preserve_residual(self):
        rng = np.random.default_rng(630729)
        ratio, input_dimension, alpha = 4, 3, 0.5
        dimension = ratio * input_dimension
        contraction = np.kron(_shift(ratio), np.eye(input_dimension))
        embedding, _ = _geometric_embedding(ratio, input_dimension, alpha)
        physical_basis = _random_unitary(rng, dimension)
        input_basis = _random_unitary(rng, input_dimension)
        changed_contraction = physical_basis @ contraction @ physical_basis.conj().T
        changed_embedding = physical_basis @ embedding @ input_basis
        residual = contraction @ embedding - alpha * embedding
        changed_residual = changed_contraction @ changed_embedding - alpha * changed_embedding
        np.testing.assert_allclose(
            changed_residual, physical_basis @ residual @ input_basis, atol=ATOL, rtol=0)
        np.testing.assert_allclose(
            changed_embedding.conj().T @ changed_embedding,
            np.eye(input_dimension), atol=ATOL, rtol=0)
        np.testing.assert_allclose(
            np.linalg.matrix_power(changed_contraction, ratio), 0, atol=ATOL, rtol=0)
        self.assertAlmostEqual(
            np.linalg.norm(changed_residual, ord=2),
            np.linalg.norm(residual, ord=2), delta=ATOL)

    def test_clean_width_corollary_and_geometric_dimension_upper_witness(self):
        self.assertAlmostEqual(_bound(0.5, 4), 1 / 54, delta=ATOL)
        for clean_qubits in range(5):
            ratio = 2 ** clean_qubits
            self.assertAlmostEqual(
                _bound(0.5, ratio), 1.5 * 3.0 ** (-ratio), delta=ATOL)
        for precision in (6, 10, 16, 32):
            clean_qubits = (precision - 1).bit_length()
            ratio = 2 ** clean_qubits
            embedding, normalization = _geometric_embedding(ratio, 1, 0.5)
            residual = np.linalg.norm(_shift(ratio) @ embedding - 0.5 * embedding)
            expected = normalization * 0.5 ** ratio
            self.assertAlmostEqual(residual / expected, 1, delta=ATOL)
            self.assertLessEqual(residual, 2.0 ** (-precision))
            required_ratio = (precision + np.log2(1.5)) / np.log2(3)
            self.assertGreaterEqual(ratio, required_ratio)

    def test_nonnilpotent_contraction_invalidates_the_obstruction(self):
        alpha, ratio, input_dimension = 0.5, 4, 2
        dimension = ratio * input_dimension
        contraction = alpha * np.eye(dimension)
        embedding = np.eye(dimension)[:, :input_dimension]
        residual = np.linalg.norm(contraction @ embedding - alpha * embedding, ord=2)
        self.assertEqual(residual, 0)
        self.assertGreater(_bound(alpha, ratio), residual)
        self.assertGreater(np.linalg.norm(np.linalg.matrix_power(contraction, dimension)), 0)

    def test_compression_of_a_nilpotent_is_not_necessarily_nilpotent(self):
        contraction = _shift(2)
        embedding = np.ones((2, 1)) / np.sqrt(2)
        compressed = embedding.conj().T @ contraction @ embedding
        np.testing.assert_allclose(contraction @ contraction, 0, atol=ATOL, rtol=0)
        np.testing.assert_allclose(compressed, [[0.5]], atol=ATOL, rtol=0)
        np.testing.assert_allclose(compressed @ compressed, [[0.25]], atol=ATOL, rtol=0)
        # Its compression has an exact scalar block, while its complete output
        # still has a nonzero residual. These are different source interfaces.
        self.assertAlmostEqual(
            np.linalg.norm(contraction @ embedding - 0.5 * embedding), 0.5, delta=ATOL)

    def test_dressed_mask_transfer_witness_independent_source_construction(self):
        x = np.array([[0, 1], [1, 0]], dtype=complex)
        y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        z = np.diag([1, -1]).astype(complex)
        for width in range(2, 7):
            with self.subTest(width=width):
                dimension = 2 ** width
                identity = np.eye(dimension, dtype=complex)
                source_basis = identity.copy()
                rotations = []
                # Build U_m=R_(m-2)...R_0 directly from its Pauli exponentials,
                # independently of the existing native gate-word fixtures.
                for j in range(width - 1):
                    generator = _pauli(width, {j: y, j + 1: x})
                    rotation = (np.cos(np.pi / 8) * identity
                                + 1j * np.sin(np.pi / 8) * generator)
                    rotations.append(rotation)
                    source_basis = rotation @ source_basis
                np.testing.assert_allclose(
                    source_basis.conj().T @ source_basis, identity, atol=ATOL, rtol=0)
                for j in range(width - 1):
                    with self.subTest(width=width, mask_bit=j):
                        mask = _pauli(width, {j: z})
                        dressed = source_basis.conj().T @ mask @ source_basis
                        witness = _pauli(width, {0: z})
                        np.testing.assert_allclose(
                            dressed @ dressed, identity, atol=ATOL, rtol=0)
                        transfer = np.trace(
                            witness @ dressed @ witness @ dressed.conj().T) / dimension
                        self.assertAlmostEqual(
                            transfer.real, 1 - 2.0 ** (-j), delta=ATOL)
                        self.assertAlmostEqual(transfer.imag, 0, delta=ATOL)

                        # Check the Clifford central factor and the shortened
                        # word separately, including its sign and order.
                        central = (mask + _pauli(width, {j: x, j + 1: x})) / np.sqrt(2)
                        cnot = identity.copy()
                        for basis in range(dimension):
                            target = basis ^ (1 << (j + 1)) if (basis >> j) & 1 else basis
                            cnot[:, basis] = identity[:, target]
                        hadamard = _pauli(width, {j: (x + z) / np.sqrt(2)})
                        np.testing.assert_allclose(
                            central, cnot @ hadamard @ cnot, atol=ATOL, rtol=0)
                        prefix = identity.copy()
                        for rotation in rotations[:j]:
                            prefix = rotation @ prefix
                        np.testing.assert_allclose(
                            dressed, prefix.conj().T @ central @ prefix,
                            atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
