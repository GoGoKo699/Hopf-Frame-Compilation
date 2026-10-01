"""Small audits of affine-tree fusion and specific source-reuse diagnostics.

Full occupied-input products are compared with outer-product fusion, direct
path maps, and independent source words. Matrices have dimension at most 64.
Support and depth-band witnesses concern only the stated restricted routes;
these tests do not establish a general compiler or gate-count lower bound.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

import numpy as np

from tests.test_one_clean_compiler import (
    _encoded_signs, _literal_mask, _paired_source_data, _paired_source_word,
    _paired_weights, _scalar_word,
)
from tests.test_operator_source_compiler import _adjoint, _pauli, _word_matrix, X
from tests.test_residual_assembly import (
    EPSILON, SCALE, _affine_data, _affine_network, _close_fixtures,
    _fixed_weighted_data,
)
from tests.test_tree_residual_structure import (
    ATOL, _addressed_frame, _edge_product, _marker, _rotation,
)
from tests.test_tree_transport import (
    _embed_rows, _path_matrix, _shift_and_injection,
)


def _local_data(height, coarse, target):
    data = _affine_data(height, coarse, target)
    gates = _affine_network(height, target, data)[2]
    return data, {node: gates[node] for node in range(1, 1 << height)}, gates[0]


def _fuse_subtree(node, size, gates):
    if node >= size // 2:
        labels = [('c', node), ('m', node)]
        return gates[node], labels, labels
    left, left_inputs, left_outputs = _fuse_subtree(2 * node, size, gates)
    right, right_inputs, right_outputs = _fuse_subtree(2 * node + 1, size, gates)
    dl, dr = len(left), len(right)
    result = np.zeros((2 + dl + dr, 2 + dl + dr), dtype=complex)
    parent = gates[node]
    result[:2, :4] = parent[:2]
    result[2:2 + dl, :4] = np.outer(left[:, 0], parent[2])
    result[2 + dl:, :4] = np.outer(right[:, 0], parent[3])
    result[2:2 + dl, 4:3 + dl] = left[:, 1:]
    result[2 + dl:, 3 + dl:] = right[:, 1:]
    inputs = [('c', node), ('m', node), ('c', 2 * node), ('c', 2 * node + 1)]
    inputs += left_inputs[1:] + right_inputs[1:]
    outputs = [('c', node), ('m', node)] + left_outputs + right_outputs
    return result, inputs, outputs


def _literal_subtree_product(node, size, gates, inputs, outputs):
    descendants = [candidate for candidate in range(node, size)
                   if candidate >> (candidate.bit_length() - node.bit_length()) == node]
    product = np.eye(len(inputs), dtype=complex)
    instructions = []
    for current in descendants:
        labels = [('c', current), ('m', current)]
        if current < size // 2:
            labels += [('c', 2 * current), ('c', 2 * current + 1)]
        rows = [inputs.index(label) for label in labels]
        product[rows] = gates[current] @ product[rows]
        instructions.append((rows, gates[current]))
    inverse = np.eye(len(inputs), dtype=complex)
    for rows, gate in reversed(instructions):
        inverse[rows] = gate.conj().T @ inverse[rows]
    order = [inputs.index(label) for label in outputs]
    return product[order], inverse[:, order], product


def _positive_fixture(height):
    size = 1 << height
    parameter = EPSILON * 2.0 ** (-height - 6)
    cosine = (1 - parameter * parameter) / (1 + parameter * parameter)
    sine = 2 * parameter / (1 + parameter * parameter)
    rotation = np.array([[cosine, -sine], [sine, cosine]], dtype=complex)
    coarse = {node: np.eye(2, dtype=complex) for node in range(1, size)}
    target = {node: rotation for node in range(1, size)}
    return coarse, target


class AffineTreeFusionTests(unittest.TestCase):
    def test_ten_mode_outer_product_merge_retains_all_occupied_inputs_and_literal_phase(self):
        for fixture, (coarse, target) in enumerate(_close_fixtures(3)):
            with self.subTest(fixture=fixture):
                _, gates, _ = _local_data(3, coarse, target)
                parent, left, right = gates[1], gates[2], gates[3]
                merged = np.block([
                    [parent[:2], np.zeros((2, 3)), np.zeros((2, 3))],
                    [np.outer(left[:, 0], parent[2]), left[:, 1:], np.zeros((4, 3))],
                    [np.outer(right[:, 0], parent[3]), np.zeros((4, 3)), right[:, 1:]],
                ])
                product = np.eye(10, dtype=complex)
                instructions = [([0, 1, 2, 3], parent), ([2, 4, 5, 6], left),
                                ([3, 7, 8, 9], right)]
                for rows, gate in instructions:
                    product[rows] = gate @ product[rows]
                output_order = [0, 1, 2, 4, 5, 6, 3, 7, 8, 9]
                np.testing.assert_allclose(merged, product[output_order], atol=ATOL, rtol=0)
                np.testing.assert_allclose(merged.conj().T @ merged, np.eye(10),
                                           atol=ATOL, rtol=0)
                # Different displayed row/column orders introduce an odd
                # permutation. The common-physical-order circuit has det +1.
                np.testing.assert_allclose(np.linalg.det(merged), -1, atol=ATOL, rtol=0)
                np.testing.assert_allclose(np.linalg.det(product), 1, atol=ATOL, rtol=0)
                inverse = np.eye(10, dtype=complex)
                for rows, gate in reversed(instructions):
                    inverse[rows] = gate.conj().T @ inverse[rows]
                np.testing.assert_allclose(inverse @ product, np.eye(10), atol=ATOL, rtol=0)
                selected = np.eye(20, dtype=complex)
                for rows, gate in instructions:
                    active = [10 + row for row in rows]
                    selected[active] = gate @ selected[active]
                expected = np.block([[np.eye(10), np.zeros((10, 10))],
                                     [np.zeros((10, 10)), product]])
                np.testing.assert_allclose(np.kron(selected, np.eye(2)),
                                           np.kron(expected, np.eye(2)), atol=ATOL, rtol=0)

    def test_recursive_two_m_mode_fusion_matches_local_factors_and_independent_paths(self):
        for height in (2, 3, 4):
            size = 1 << height
            for fixture, (coarse, target) in enumerate(_close_fixtures(height)):
                with self.subTest(height=height, fixture=fixture):
                    (blocks, rho, _, _), gates, _ = _local_data(height, coarse, target)
                    fused, inputs, outputs = _fuse_subtree(1, size, gates)
                    literal, inverse, physical = _literal_subtree_product(
                        1, size, gates, inputs, outputs)
                    self.assertEqual(len(fused), 2 * (size - 1))
                    np.testing.assert_allclose(fused, literal, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(inverse @ fused, np.eye(len(fused)),
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(np.linalg.det(physical), 1, atol=ATOL, rtol=0)
                    stop_rows = [outputs.index(('c', node)) for node in range(1, size)]
                    marker_columns = [inputs.index(('m', node)) for node in range(1, size)]
                    expected = np.zeros((size - 1, size - 1), dtype=complex)
                    continuation = np.zeros(size - 1, dtype=complex)
                    for row, node in enumerate(range(1, size)):
                        continuation[row] = blocks[node][1, 0] * _edge_product(target, 1, node)
                        for column, ancestor in enumerate(range(1, size)):
                            if ancestor == node:
                                expected[row, column] = blocks[node][1, 1]
                            else:
                                path = _edge_product(target, ancestor, node, True)
                                if path is not None:
                                    expected[row, column] = blocks[node][1, 0] * path
                    np.testing.assert_allclose(fused[np.ix_(stop_rows, marker_columns)],
                                               expected / SCALE, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(fused[stop_rows, 0] * np.sqrt(rho[1]),
                                               continuation, atol=ATOL, rtol=0)

    def test_positive_tree_has_full_frontier_rank_and_growing_exact_column_support(self):
        for height in (2, 3, 4):
            size = 1 << height
            coarse, target = _positive_fixture(height)
            (blocks, rho, _, _), gates, root_gate = _local_data(height, coarse, target)
            paths = _path_matrix(height, target)
            partial = np.eye(2 * size, dtype=complex)
            partial[[0, size, size + 1]] = root_gate @ partial[[0, size, size + 1]]
            heap_markers = [0] + [_marker(node, height) for node in range(1, size)]
            physical_columns = np.argsort(heap_markers)
            band = np.eye(2 * size, dtype=complex)
            for depth in range(height - 1):
                for node in range(1 << depth, 1 << (depth + 1)):
                    rows = [size + node, node, size + 2 * node, size + 2 * node + 1]
                    partial[rows] = gates[node] @ partial[rows]
                    band[rows] = gates[node] @ band[rows]
                next_nodes = list(range(1 << (depth + 1), 1 << (depth + 2)))
                frontier = partial[[size + node for node in next_nodes]]
                np.testing.assert_allclose(frontier @ frontier.conj().T,
                                           np.eye(len(next_nodes)), atol=ATOL, rtol=0)
                accepted = frontier[:, physical_columns]
                expected = np.array([np.sqrt(rho[node]) / SCALE * paths[node - 1]
                                     for node in next_nodes])
                np.testing.assert_allclose(accepted, expected, atol=ATOL, rtol=0)
                self.assertEqual(np.linalg.matrix_rank(accepted), len(next_nodes))
                # No threshold truncation is used: the exact support concerns
                # this positive family, not an approximation lower bound.
                self.assertEqual(np.count_nonzero(band[:, size + 1]),
                                 3 * (1 << (depth + 1)) - 2)
            h = np.array([blocks[node][1, 0] for node in range(1, size)])
            forward = _embed_rows(height, h[:, None] * paths)
            diagonal = np.empty(size, dtype=complex)
            diagonal[0] = blocks[1][0, 0]
            for node in range(1, size):
                diagonal[_marker(node, height)] = blocks[node][1, 1]
            affine = np.diag(diagonal) + forward
            self.assertEqual(np.count_nonzero(affine), height * size + 1)
            self.assertEqual(np.count_nonzero(forward.conj().T), (height - 1) * size + 1)

    def test_balanced_terminal_witness_has_undamped_transport_but_a_one_angle_frame_formula(self):
        parameter = EPSILON / 16
        cosine = (1 - parameter ** 2) / (1 + parameter ** 2)
        sine = 2 * parameter / (1 + parameter ** 2)
        theta = 2 * np.arctan(parameter)
        for height in (2, 3, 4):
            size = 1 << height
            coarse = {node: _rotation(np.pi / 4) for node in range(1, size)}
            target = {node: _rotation(np.pi / 4 + theta) if node >= size // 2
                      else coarse[node] for node in range(1, size)}
            (blocks, rho, _, _), gates, _ = _local_data(height, coarse, target)
            expected_rho = sine ** 2 * SCALE ** 2 / (SCALE ** 2 - cosine ** 2)
            np.testing.assert_allclose([rho[node] for node in range(1, size)],
                                       expected_rho, atol=ATOL, rtol=0)
            reverse_data = _fixed_weighted_data(height, target, coarse)
            np.testing.assert_allclose([reverse_data[4][node] for node in range(1, size)],
                                       sine ** 2, atol=ATOL, rtol=0)
            for node in range(1, size // 2):
                np.testing.assert_allclose(gates[node][:, 0],
                                           [0, 0, 1 / np.sqrt(2), 1 / np.sqrt(2)],
                                           atol=ATOL, rtol=0)
            shift, injection = _shift_and_injection(height, target)
            root = np.eye(size - 1)[:, 0]
            self.assertAlmostEqual(np.linalg.norm(np.linalg.matrix_power(shift, height - 1) @ root),
                                   1, delta=ATOL)
            paths = _path_matrix(height, target)
            stops = np.zeros((size - 1, size - 1))
            stops[size // 2 - 1:, size // 2 - 1:] = np.eye(size // 2)
            transport = _embed_rows(height, stops @ paths)
            even_projector = np.diag([int(index % 2 == 0) for index in range(size)])
            np.testing.assert_allclose(transport.conj().T @ transport, even_projector,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(transport @ transport.conj().T, np.eye(size) - even_projector,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(transport @ transport, 0, atol=ATOL, rtol=0)
            weights = np.array([blocks[node][1, 0] for node in range(1, size)])
            forward = _embed_rows(height, weights[:, None] * paths)
            truncated = np.zeros_like(paths)
            power = np.eye(size - 1, dtype=complex)
            for _ in range(height - 1):
                truncated += power @ injection
                power = power @ shift
            missed = _embed_rows(height, weights[:, None] * (paths - truncated))[:, 0]
            self.assertAlmostEqual(np.linalg.norm(missed), sine, delta=ATOL)
            np.testing.assert_allclose(forward, sine * transport, atol=ATOL, rtol=0)
            c_frame = _addressed_frame(height, coarse)
            w_frame = _addressed_frame(height, target)
            np.testing.assert_allclose(w_frame,
                                       np.kron(np.eye(size // 2), _rotation(theta)) @ c_frame,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(c_frame.conj().T @ w_frame,
                                       cosine * np.eye(size) + sine * (transport - transport.conj().T),
                                       atol=ATOL, rtol=0)

    def test_actual_paired_loader_correlation_matches_programmed_signed_moments(self):
        root = (np.sqrt(5) - 1) / 8
        self.assertAlmostEqual(16 * root * root + 4 * root - 1, 0, delta=ATOL)
        for q in (2, 3, 4, 5):
            core = q + 1
            signs = _encoded_signs(q, np.pi / 3)
            loader = _paired_source_word(q)
            transformed_mask = _word_matrix(
                core, loader + _literal_mask(signs) + _adjoint(loader))
            x0 = _pauli(core, {0: X})
            correlation = np.trace(x0 @ transformed_mask.conj().T @ x0 @ transformed_mask) / (1 << core)
            weights, _ = _paired_weights(q)
            signed_moment = float(weights @ (1 - 2 * signs))
            np.testing.assert_allclose(correlation, signed_moment, atol=ATOL, rtol=0)
            self.assertLessEqual(abs(signed_moment - root), (15 / 8) * 2.0 ** (-q) + ATOL)
            exact_moment = sum(Fraction(float(weight)) * int(1 - 2 * sign)
                               for weight, sign in zip(weights, signs))
            polynomial = 16 * exact_moment ** 2 + 4 * exact_moment - 1
            self.assertNotEqual(polynomial, 0)

    def test_changed_address_and_reused_fixed_scalar_keep_their_unwanted_actions(self):
        # Native Ry(pi/4)=HZ on address x=0; z=1 is arbitrary dirty work.
        phase_mask = [('H', 1), ('CX', 0, 1), ('H', 1)]
        addressed_rotation = [('Z', 0), ('H', 0)]
        routed = _word_matrix(2, phase_mask + addressed_rotation + phase_mask)
        expected = np.block([[_rotation(np.pi / 4), np.zeros((2, 2))],
                             [np.zeros((2, 2)), _rotation(-np.pi / 4)]])
        np.testing.assert_allclose(routed, expected, atol=ATOL, rtol=0)
        unchanged_dirty = _word_matrix(2, addressed_rotation)
        self.assertAlmostEqual(np.linalg.norm(routed - unchanged_dirty, ord=2),
                               2 * np.sin(np.pi / 4), delta=ATOL)
        for q in (2, 3):
            core = q + 1
            fixed = _paired_source_data(q)[3]
            scalar_word = _scalar_word(q, _literal_mask(fixed), core)
            scalar = _word_matrix(core + 1, scalar_word)
            twice = _word_matrix(core + 1, scalar_word + scalar_word)
            accepted = scalar[:1 << core, :1 << core]
            np.testing.assert_allclose(accepted, np.eye(1 << core) / 2, atol=ATOL, rtol=0)
            np.testing.assert_allclose(twice[:1 << core, :1 << core],
                                       -np.eye(1 << core) / 2, atol=ATOL, rtol=0)
            self.assertAlmostEqual(np.linalg.norm(
                twice[:1 << core, :1 << core] - accepted @ accepted, ord=2), 3 / 4,
                delta=ATOL)
            inverse = _word_matrix(core + 1, _adjoint(scalar_word))
            np.testing.assert_allclose(inverse @ scalar, np.eye(1 << (core + 1)),
                                       atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
