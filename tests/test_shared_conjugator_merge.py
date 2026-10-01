"""A full-word audit of the common-parity source-conjugator proposal.

The one-clean scalar conjugator can be shared by every physical Hopf pair.
Its literal cancellation is valid, but rejected amplitudes return at the
first fork.  Local inactive modes remain identity, so even the product of
the local accepted blocks has no uniform global beta normalization.  This
fixture rules out this unamplified word, not an arbitrary residual merger.
"""
from __future__ import annotations

import unittest
import itertools

import numpy as np

from tests.test_one_clean_compiler import (
    _coefficients, _encoded_signs, _literal_mask, _paired_source_data,
    _paired_weights, _scalar_word,
)
from tests.test_operator_source_compiler import (
    X, _adjoint, _expand_toffolis, _word_matrix,
)


ATOL = 8e-11
EDGES = ((0, 2), (0, 1), (2, 3))
ANGLES = (.23, -.51, .87)
PARITY = np.diag([1, -1, -1, 1])


def _edge_matrices(low, high):
    swap = np.zeros((4, 4))
    swap[low, high] = swap[high, low] = 1
    return swap @ swap, swap, swap @ PARITY


def _global_a_word(q):
    """One fixed scalar source conjugated by flag-controlled logical parity."""
    core, flag = q + 1, q + 3
    fixed = _paired_weights(q)[1]
    routing = []
    for target in (core, core + 1):
        routing += [('H', target), ('CX', flag, target), ('H', target)]
    return routing + _scalar_word(q, _literal_mask(fixed), flag) + _adjoint(routing)


def _edge_b_word(q, edge, signs):
    """Conditional scalar word, inactive identity, with literal X routing."""
    core, flag = q + 1, q + 3
    low, high = edge
    target_bit = (low ^ high).bit_length() - 1
    address_bit = 1 - target_bit
    target, address = core + target_bit, core + address_bit
    negative = [('X', address)] if not ((low >> address_bit) & 1) else []
    scalar = _scalar_word(q, _literal_mask(signs), flag, predicate=(address,))
    routing = [('CX', flag, target)]
    return (negative + routing + scalar + routing + negative)


def _native_complex_coarse_word():
    """Three literal addressed determinant-one complex diagonal local words.

    Each active pair receives diag(exp(-i*pi/4), exp(i*pi/4)).  This is a
    native Hopf coarse *word*, without asserting a small coarse-error bound.
    """
    word = []
    for low, high in EDGES:
        target = (low ^ high).bit_length() - 1
        address = 1 - target
        negative = [('X', address)] if not ((low >> address) & 1) else []
        word += negative + [('T', target), ('CX', address, target),
                            ('TDG', target), ('CX', address, target)] + negative
    return word


class SharedConjugatorMergeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.q = 3
        cls.core = cls.q + 1
        cls.core_dimension = 1 << cls.core
        cls.total = cls.core + 3  # Core, two logical bits, one used flag.
        cls.source, _, cls.weights, cls.fixed = _paired_source_data(cls.q)
        fixed_mask = _word_matrix(cls.core, _literal_mask(cls.fixed))
        fixed_source = fixed_mask @ cls.source @ fixed_mask.conj().T
        cls.df = (cls.source @ fixed_source - fixed_source @ cls.source) / 2
        cls.a_word = _global_a_word(cls.q)
        cls.a = _word_matrix(cls.total, _expand_toffolis(cls.a_word))
        cls.b_words, cls.bs, cls.qs, cls.es, cls.ks, cls.kappas = [], [], [], [], [], []
        for edge, angle in zip(EDGES, ANGLES):
            # X_e Z_parity has the opposite orientation on the (2,3) pair.
            oriented = PARITY[edge[0], edge[0]] * angle
            signs = _encoded_signs(cls.q, oriented)
            s, r = _coefficients(cls.weights, cls.fixed, signs)
            mask = _word_matrix(cls.core, _literal_mask(signs))
            programmed = mask @ cls.source @ mask.conj().T
            dg = (cls.source @ programmed - programmed @ cls.source) / 2
            k = dg + 2 * r * cls.df
            projector, swap, generator = _edge_matrices(*edge)
            e = np.eye(4) + (s - 1) * projector + r * generator
            b_word = _edge_b_word(cls.q, edge, signs)
            b = _word_matrix(cls.total, _expand_toffolis(b_word))
            q_word = cls.a_word + b_word + _adjoint(cls.a_word)
            actual = _word_matrix(cls.total, _expand_toffolis(q_word))
            cls.b_words.append(b_word)
            cls.bs.append(b)
            cls.qs.append(actual)
            cls.es.append(e)
            cls.ks.append(k)
            cls.kappas.append(np.sqrt(1 - s * s - r * r))

    def test_global_conjugator_and_each_native_boundary_all_ports(self):
        dim = 1 << self.total
        expected_a = (.5 * np.eye(dim)
                      + np.kron(X, np.kron(PARITY, self.df)))
        np.testing.assert_allclose(self.a, expected_a, atol=ATOL, rtol=0)
        for edge, actual, e, k, kappa in zip(
                EDGES, self.qs, self.es, self.ks, self.kappas):
            with self.subTest(edge=edge):
                swap = _edge_matrices(*edge)[1]
                expected = (np.kron(np.eye(2), np.kron(e, np.eye(self.core_dimension)))
                            + np.kron(X, np.kron(swap, k)))
                # A matrix equality checks occupied flag ports and all dirty columns.
                np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                np.testing.assert_allclose(actual.conj().T @ actual,
                                           np.eye(dim), atol=ATOL, rtol=0)
                np.testing.assert_allclose(k.conj().T, -k, atol=ATOL, rtol=0)
                np.testing.assert_allclose(k.conj().T @ k,
                                           kappa ** 2 * np.eye(self.core_dimension),
                                           atol=ATOL, rtol=0)

    def test_exact_full_word_telescope_with_unused_second_flag(self):
        fork = self.qs[2] @ self.qs[1] @ self.qs[0]
        fused_word = (self.a_word + self.b_words[0] + self.b_words[1]
                      + self.b_words[2] + _adjoint(self.a_word))
        fused = _word_matrix(self.total, _expand_toffolis(fused_word))
        np.testing.assert_allclose(fused, fork, atol=ATOL, rtol=0)
        np.testing.assert_allclose(
            fork, self.a.conj().T @ self.bs[2] @ self.bs[1] @ self.bs[0] @ self.a,
            atol=ATOL, rtol=0)
        # A second available flag is untouched on both its accepted and occupied ports.
        np.testing.assert_allclose(np.kron(np.eye(2), fused),
                                   np.kron(np.eye(2), fork), atol=ATOL, rtol=0)

    def test_fork_rejected_returns_have_constant_full_core_defect(self):
        block_size = 4 * self.core_dimension
        fork = self.qs[2] @ self.qs[1] @ self.qs[0]
        accepted = fork[:block_size, :block_size]
        local_product = np.kron(self.es[2] @ self.es[1] @ self.es[0],
                                np.eye(self.core_dimension))
        left = np.zeros((4, 4))
        right = np.zeros((4, 4))
        left[1, 2] = right[3, 0] = 1
        defect = (np.kron(left, self.ks[1] @ self.ks[0])
                  + np.kron(right, self.ks[2] @ self.ks[0]))
        np.testing.assert_allclose(accepted - local_product, defect,
                                   atol=ATOL, rtol=0)
        expected = max(self.kappas[0] * self.kappas[1],
                       self.kappas[0] * self.kappas[2])
        self.assertAlmostEqual(np.linalg.norm(defect, 2), expected, delta=ATOL)
        self.assertGreater(expected, .8)
        # Each dirty input on a supported logical column attains its branch norm.
        for logical, child in ((2, 1), (0, 2)):
            columns = defect[:, logical * self.core_dimension:(logical + 1) * self.core_dimension]
            np.testing.assert_allclose(
                columns.conj().T @ columns,
                (self.kappas[0] * self.kappas[child]) ** 2 * np.eye(self.core_dimension),
                atol=ATOL, rtol=0)

    def test_literal_complex_native_coarse_wrapper_preserves_defect(self):
        coarse = _word_matrix(2, _native_complex_coarse_word())
        np.testing.assert_allclose(coarse.conj().T @ coarse, np.eye(4),
                                   atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.det(coarse), 1, delta=ATOL)
        self.assertGreater(np.linalg.norm(coarse.imag), .5)
        block_size = 4 * self.core_dimension
        fork = self.qs[2] @ self.qs[1] @ self.qs[0]
        actual = fork[:block_size, :block_size]
        local_product = np.kron(self.es[2] @ self.es[1] @ self.es[0],
                                np.eye(self.core_dimension))
        wrapper = np.kron(coarse.conj().T, np.eye(self.core_dimension))
        self.assertAlmostEqual(np.linalg.norm(wrapper @ (actual - local_product), 2),
                               np.linalg.norm(actual - local_product, 2), delta=ATOL)

    def test_candidate_fails_the_half_unitary_target_on_a_full_dirty_column(self):
        target = np.eye(4)
        for edge, angle in zip(EDGES, ANGLES):
            projector, _, generator = _edge_matrices(*edge)
            oriented = PARITY[edge[0], edge[0]] * angle
            rotation = (np.eye(4) + (np.cos(oriented) - 1) * projector
                        + np.sin(oriented) * generator)
            target = rotation @ target
        block_size = 4 * self.core_dimension
        accepted = (self.qs[2] @ self.qs[1] @ self.qs[0])[:block_size, :block_size]
        local_product = self.es[2] @ self.es[1] @ self.es[0]
        scalar_difference = local_product[1, 2] - target[1, 2] / 2
        lower = self.kappas[0] * self.kappas[1] - abs(scalar_difference)
        self.assertGreater(lower, .8)
        error = accepted - np.kron(target / 2, np.eye(self.core_dimension))
        witness = error[self.core_dimension:2 * self.core_dimension,
                        2 * self.core_dimension:3 * self.core_dimension]
        np.testing.assert_allclose(
            witness, scalar_difference * np.eye(self.core_dimension) + self.ks[1] @ self.ks[0],
            atol=ATOL, rtol=0)
        self.assertGreaterEqual(np.linalg.svd(witness, compute_uv=False)[-1] + ATOL, lower)

    def test_arbitrary_program_masks_obey_the_feasible_ellipse(self):
        weights, fixed = _paired_weights(2)
        for bits in itertools.product((0, 1), repeat=len(weights)):
            signs = np.array(bits)
            s, r = _coefficients(weights, fixed, signs)
            self.assertLessEqual(s * s + 4 * r * r / 3, 1 + ATOL)

    def test_untouched_rows_and_retuning_bound_for_all_small_mask_pairs(self):
        block_size = 4 * self.core_dimension
        accepted = (self.qs[2] @ self.qs[1] @ self.qs[0])[:block_size, :block_size]
        s0, r0, sl = self.es[0][0, 0], self.es[0][2, 0], self.es[1][0, 0]
        d = self.core_dimension
        np.testing.assert_allclose(accepted[d:2 * d, d:2 * d],
                                   sl * np.eye(d), atol=ATOL, rtol=0)
        root_columns = np.concatenate((accepted[:d, :d], accepted[:d, 2 * d:3 * d]), axis=1)
        np.testing.assert_allclose(root_columns,
                                   np.kron([[sl * s0, -sl * r0]], np.eye(d)),
                                   atol=ATOL, rtol=0)
        # The row and spectator constraints alone force a positive error for
        # every source width; this finite sweep is only an indexing check.
        weights, fixed = _paired_weights(2)
        pairs = [_coefficients(weights, fixed, np.array(bits))
                 for bits in itertools.product((0, 1), repeat=len(weights))]
        theta_root = theta_left = np.pi / 4
        target_scale = .5 * np.cos(theta_left)
        direction = np.array([np.cos(theta_root), np.sin(theta_root)])
        lower = abs(target_scale) * (1 - np.sqrt(1 - np.sin(theta_root) ** 2 / 4)) / 2
        self.assertGreater(lower, .011)
        for (s0, r0), (sl, _) in itertools.product(pairs, repeat=2):
            row_error = np.linalg.norm(sl * np.array([s0, r0]) - target_scale * direction)
            spectator_error = abs(sl - target_scale)
            self.assertGreaterEqual(max(row_error, spectator_error) + ATOL, lower)

    def test_precision_limit_and_nonuniform_local_normalization(self):
        beta = np.sin(np.pi / 10)
        ideal_blocks = []
        for edge, angle in zip(EDGES, ANGLES):
            projector, _, generator = _edge_matrices(*edge)
            oriented = PARITY[edge[0], edge[0]] * angle
            ideal_blocks.append(np.eye(4) + (beta * np.cos(oriented) - 1) * projector
                                + beta * np.sin(oriented) * generator)
        singular = np.linalg.svd(ideal_blocks[2] @ ideal_blocks[1] @ ideal_blocks[0],
                                compute_uv=False)
        np.testing.assert_allclose(singular, [beta, beta, beta ** 2, beta ** 2],
                                   atol=ATOL, rtol=0)
        # Large-q scalar arithmetic, without exponentially large core matrices.
        weights, fixed = _paired_weights(24)
        kappas = []
        for edge, angle in zip(EDGES, ANGLES):
            signs = _encoded_signs(24, PARITY[edge[0], edge[0]] * angle)
            s, r = _coefficients(weights, fixed, signs)
            kappas.append(np.sqrt(1 - s * s - r * r))
        limit = 1 - beta ** 2
        self.assertAlmostEqual(max(kappas[0] * kappas[1], kappas[0] * kappas[2]),
                               limit, delta=1e-6)

    def test_literal_telescope_retains_linear_precision_cost_per_program(self):
        counts = []
        for q in (3, 4):
            a_word = _global_a_word(q)
            bs = [_edge_b_word(q, edge, _encoded_signs(
                q, PARITY[edge[0], edge[0]] * angle)) for edge, angle in zip(EDGES, ANGLES)]
            word = _expand_toffolis(a_word + sum(bs, []) + _adjoint(a_word))
            counts.append(sum(gate[0] in ('T', 'TDG') for gate in word))
        # Three M calls in each of the three B programs and two outside A factors;
        # each M contains U and U^dagger, each with 2q T gates.  Controls add a
        # q-independent count.  This is the charged literal word, not a lower bound.
        self.assertEqual(counts[1] - counts[0], 12 * (len(EDGES) + 2))


if __name__ == '__main__':
    unittest.main()
