"""Small checks of weighted transport and candidate block completions.

Weighted operators have dimension at most 16 and complete dilations at most
32. These checks distinguish a bounded operator norm from orthogonal stopping
columns, and verify a logical one-flag dilation. They do not provide an
elementary Clifford+T implementation or an asymptotic gate bound.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_tree_residual_structure import (
    ATOL, _addressed_frame, _edge_product, _marker, _rotation,
)
from tests.test_tree_transport import (
    _embed_rows, _overlap_generators, _path_matrix, _word_fixtures,
)


def _riccati_data(height, coarse, target):
    """Evaluate the scalar recurrence; no matrix inverses enter this helper."""
    size = 1 << height
    blocks = _overlap_generators(height, coarse, target)
    weights = {node: blocks[node][1, 0] for node in range(1, size)}
    defects = {node: 1 - abs(blocks[node][0, 0]) ** 2
               for node in range(1, size)}
    # The floor gives a positive normalization for coincident frames, including
    # exact zero defects. It is a fixture choice, not a new endpoint parameter.
    gamma = max(0.0, max(defects.values()))
    alpha = 4 * np.sqrt(max(gamma, 1e-8))
    rho = {leaf: 0.0 for leaf in range(size, 2 * size)}
    beta, cross = {}, {}
    for node in range(size - 1, 0, -1):
        word = target[node]
        children = np.array([rho[2 * node], rho[2 * node + 1]])
        continuation = np.vdot(word[:, 0], children * word[:, 0]).real
        injection = np.vdot(word[:, 1], children * word[:, 1]).real
        cross[node] = np.vdot(word[:, 0], children * word[:, 1])
        beta[node] = alpha ** 2 - injection
        rho[node] = (abs(weights[node]) ** 2 + continuation
                     + abs(cross[node]) ** 2 / beta[node])
    return weights, defects, gamma, alpha, rho, beta, cross


def _complete_special_unitary(first, second):
    """Complete two prescribed columns using canonical vectors, with det one."""
    columns = [first.copy(), second.copy()]
    for seed in np.eye(len(first), dtype=complex):
        candidate = seed.copy()
        for _ in range(2):
            for column in columns:
                candidate -= np.vdot(column, candidate) * column
        norm = np.linalg.norm(candidate)
        if norm > 1e-10:
            columns.append(candidate / norm)
        if len(columns) == len(first):
            break
    result = np.column_stack(columns)
    determinant = np.linalg.det(result)
    result[:, -1] *= determinant.conjugate() / abs(determinant)
    return result


def _riccati_network(height, target, data):
    """Build full gates and their reversed adjoints, including both flag sectors."""
    weights, _, _, alpha, rho, beta, cross = data
    size = 1 << height
    gates = []

    def carrier(node):
        if node == 1:
            return 0
        return node // 2 if node % 2 == 0 else size + node - 2

    accepted = np.sqrt(rho[1]) / alpha
    rejected = np.sqrt(1 - accepted ** 2)
    gates.append(([0, size], np.array([
        [accepted, -rejected], [rejected, accepted],
    ], dtype=complex)))
    for node in range(1, size):
        if node < size // 2:
            word = target[node]
            children = np.sqrt([rho[2 * node], rho[2 * node + 1]])
            if rho[node] == 0:
                first = np.array([1, 0, 0, 0], dtype=complex)
            else:
                first = np.r_[weights[node], children * word[:, 0],
                              -cross[node].conjugate() / np.sqrt(beta[node])]
                first /= np.sqrt(rho[node])
            second = np.r_[0, children * word[:, 1], np.sqrt(beta[node])] / alpha
            gate = _complete_special_unitary(first, second)
            rows = [carrier(node), node, size + 2 * node - 1, size + 2 * node]
        else:
            phase = weights[node] / abs(weights[node]) if weights[node] != 0 else 1
            gate = np.diag([phase, np.conj(phase)])
            rows = [carrier(node), node]
        gates.append((rows, gate))

    network = np.eye(2 * size, dtype=complex)
    for rows, gate in gates:
        network[rows] = gate @ network[rows]
    inverse = np.eye(2 * size, dtype=complex)
    for rows, gate in reversed(gates):
        inverse[rows] = gate.conj().T @ inverse[rows]

    accepted_rows = [2 * size - 1] + [carrier(node) for node in range(1, size)]
    used = set(accepted_rows)
    row_order = accepted_rows + [row for row in range(2 * size) if row not in used]
    network = network[row_order]
    inverse = inverse[:, row_order]
    # Both signal sectors use the repository's prescribed physical marker
    # order. This is an explicit permutation, not a free omitted operation.
    heap_to_physical = [0] + [_marker(node, height) for node in range(1, size)]
    permutation = np.argsort(heap_to_physical)
    full_order = np.r_[permutation, size + permutation]
    network = network[np.ix_(full_order, full_order)]
    inverse = inverse[np.ix_(full_order, full_order)]
    return network, inverse, gates


def _dilation_fixtures(height):
    native, generic, singular = _word_fixtures(height)
    identity = {node: np.eye(2, dtype=complex) for node in range(1, 1 << height)}
    return ((native, generic), (native, singular), (singular, native),
            (generic, generic), (identity, identity))


class WeightedTransportBlockTests(unittest.TestCase):
    def test_weighted_gram_from_child_defects_and_ancestor_paths(self):
        for height in (2, 3, 4):
            size = 1 << height
            native, generic, singular = _word_fixtures(height)
            for fixture, (coarse, target) in enumerate((
                    (native, generic), (generic, singular), (singular, native),
                    (singular, singular))):
                with self.subTest(height=height, fixture=fixture):
                    blocks = _overlap_generators(height, coarse, target)
                    defects = {node: 1 - abs(blocks[node][0, 0]) ** 2
                               for node in range(1, size)}
                    defects.update({leaf: 0.0 for leaf in range(size, 2 * size)})
                    h = np.array([blocks[node][1, 0] for node in range(1, size)])
                    weighted = h[:, None] * _path_matrix(height, target)

                    # Independently predict every Gram entry from the two
                    # child defects, rather than summing weighted descendants.
                    expected = np.zeros((size, size), dtype=complex)
                    expected[0, 0] = defects[1]
                    for node in range(1, size):
                        marker = _marker(node, height)
                        child_defects = np.array([defects[2 * node],
                                                  defects[2 * node + 1]])
                        word = target[node]
                        zeta = np.sum(abs(word[:, 1]) ** 2 * child_defects)
                        tau = np.sum(word[:, 0].conj() * word[:, 1]
                                     * child_defects)
                        block = blocks[node]
                        np.testing.assert_allclose(
                            tau, -(block[0, 0].conj() * block[0, 1]
                                   + block[1, 0].conj() * block[1, 1]),
                            atol=ATOL, rtol=0)
                        expected[marker, marker] = zeta
                        expected[0, marker] = np.conj(
                            _edge_product(target, 1, node)) * tau
                        expected[marker, 0] = expected[0, marker].conj()
                        for ancestor in range(1, size):
                            if ancestor == node:
                                continue
                            path = _edge_product(target, ancestor, node, True)
                            if path is None:
                                continue
                            ancestor_marker = _marker(ancestor, height)
                            expected[ancestor_marker, marker] = np.conj(path) * tau
                            expected[marker, ancestor_marker] = expected[
                                ancestor_marker, marker].conj()

                    np.testing.assert_allclose(
                        weighted.conj().T @ weighted, expected, atol=ATOL, rtol=0)

    def test_close_native_identity_witness_has_nonorthogonal_stopping_columns(self):
        height = 2
        size = 1 << height
        coarse = {node: np.eye(2, dtype=complex) for node in range(1, size)}
        coarse_allowance = 1 / 64
        # Finite members of the analytic family t -> 0. Its exact formulas,
        # rather than extrapolation from these cases, give the limiting overlap.
        for angle in (1 / 128, 1 / 512, 1 / 2048):
            with self.subTest(angle=angle):
                target = {1: _rotation(angle), 2: _rotation(angle),
                          3: np.eye(2, dtype=complex)}
                cosine, sine = np.cos(angle), np.sin(angle)
                blocks = _overlap_generators(height, coarse, target)
                h = np.array([blocks[node][1, 0] for node in range(1, size)])
                weighted = h[:, None] * _path_matrix(height, target)
                logical_order = [0] + [_marker(node, height)
                                       for node in range(1, size)]
                weighted = weighted[:, logical_order]
                expected = np.array([
                    [sine, 0, 0, 0],
                    [cosine * sine, -sine ** 2, 0, 0],
                    [0, 0, 0, 0],
                ], dtype=complex)
                np.testing.assert_allclose(weighted, expected, atol=0, rtol=ATOL)

                active_gram = weighted[:, :2].conj().T @ weighted[:, :2]
                expected_gram = np.array([
                    [sine ** 2 * (1 + cosine ** 2), -cosine * sine ** 3],
                    [-cosine * sine ** 3, sine ** 4],
                ])
                np.testing.assert_allclose(
                    active_gram, expected_gram, atol=0, rtol=ATOL)
                expected_eigenvalues = sine ** 2 * np.array([1 - cosine,
                                                            1 + cosine])
                np.testing.assert_allclose(
                    np.linalg.eigvalsh(active_gram), expected_eigenvalues,
                    atol=1e-20, rtol=1e-8)
                normalized = weighted[:, :2] / np.linalg.norm(weighted[:, :2], axis=0)
                overlap = np.vdot(normalized[:, 0], normalized[:, 1])
                np.testing.assert_allclose(
                    overlap, -cosine / np.sqrt(1 + cosine ** 2), atol=ATOL, rtol=0)
                self.assertGreater(abs(overlap), 0.7)

                # The root column already saturates gamma. Its nonzero cross
                # term forces a norm strictly larger than sqrt(gamma).
                gamma = max(1 - abs(blocks[node][0, 0]) ** 2
                            for node in range(1, size))
                stable_gamma = sine ** 2 * (1 + cosine ** 2)
                self.assertAlmostEqual(gamma, stable_gamma, delta=ATOL)
                norm = np.linalg.norm(weighted, ord=2)
                np.testing.assert_allclose(
                    norm, sine * np.sqrt(1 + cosine), atol=ATOL, rtol=0)
                self.assertGreater(norm, np.sqrt(gamma))
                self.assertGreater(norm, np.sqrt(stable_gamma))

                # This counterexample remains inside the retained constant
                # subtree-error contract, including both prescribed columns.
                frame_error = np.linalg.norm(
                    _addressed_frame(height, target) - np.eye(size), ord=2)
                self.assertLessEqual(frame_error, 2 * angle + ATOL)
                self.assertLessEqual(frame_error, coarse_allowance)
                for node in (2, 3):
                    self.assertLessEqual(
                        np.linalg.norm(target[node] - coarse[node], ord=2),
                        coarse_allowance)

    def test_scalar_riccati_values_match_independent_subtree_schur_complements(self):
        for height in (1, 2, 3, 4):
            size = 1 << height
            for fixture, (coarse, target) in enumerate(_dilation_fixtures(height)):
                data = _riccati_data(height, coarse, target)
                weights, defects, gamma, alpha, rho, beta, _ = data
                with self.subTest(height=height, fixture=fixture):
                    for root in range(1, size):
                        descendants = [node for node in range(root, size)
                                       if _edge_product(target, root, node) is not None]
                        root_column = np.array([
                            weights[node] * _edge_product(target, root, node)
                            for node in descendants])
                        markers = np.zeros((len(descendants), len(descendants)),
                                           dtype=complex)
                        for row, node in enumerate(descendants):
                            for column, ancestor in enumerate(descendants):
                                if ancestor == node:
                                    continue
                                path = _edge_product(target, ancestor, node, True)
                                if path is not None:
                                    markers[row, column] = weights[node] * path
                        # This dense Schur solve uses direct path products and
                        # is independent of the bottom-up scalar recurrence.
                        gram = markers.conj().T @ markers
                        schur_matrix = alpha ** 2 * np.eye(len(descendants)) - gram
                        self.assertGreaterEqual(
                            np.linalg.eigvalsh(schur_matrix)[0],
                            alpha ** 2 - 4 * gamma - ATOL)
                        coupling = markers.conj().T @ root_column
                        independent = (np.vdot(root_column, root_column)
                                       + np.vdot(coupling,
                                                 np.linalg.solve(schur_matrix, coupling)))
                        np.testing.assert_allclose(
                            rho[root], independent, atol=ATOL, rtol=0)
                        self.assertLessEqual(rho[root], 4 * defects[root] / 3 + ATOL)
                        self.assertGreaterEqual(beta[root], 11 * alpha ** 2 / 12 - ATOL)
                        self.assertGreaterEqual(rho[root], 0)

                    if fixture == 4:
                        self.assertEqual(gamma, 0)
                        self.assertTrue(all(value == 0 for value in rho.values()))

    def test_complete_one_flag_network_and_literal_reverse_on_all_inputs(self):
        for height in (1, 2, 3, 4):
            size = 1 << height
            for fixture, (coarse, target) in enumerate(_dilation_fixtures(height)):
                with self.subTest(height=height, fixture=fixture):
                    data = _riccati_data(height, coarse, target)
                    weights, _, _, alpha, _, _, _ = data
                    network, inverse, gates = _riccati_network(height, target, data)
                    weighted = np.array([weights[node] for node in range(1, size)])
                    weighted = weighted[:, None] * _path_matrix(height, target)
                    expected = _embed_rows(height, weighted) / alpha
                    np.testing.assert_allclose(
                        network[:size, :size], expected, atol=ATOL, rtol=0)

                    # Verify every input column, including the occupied/rejected
                    # signal sector, instead of testing initialized inputs only.
                    np.testing.assert_allclose(
                        network.conj().T @ network, np.eye(2 * size),
                        atol=ATOL, rtol=0)
                    np.testing.assert_allclose(
                        inverse @ network, np.eye(2 * size), atol=ATOL, rtol=0)
                    np.testing.assert_allclose(
                        network @ inverse, np.eye(2 * size), atol=ATOL, rtol=0)
                    np.testing.assert_allclose(
                        inverse, network.conj().T, atol=ATOL, rtol=0)
                    for _, gate in gates:
                        np.testing.assert_allclose(
                            gate.conj().T @ gate, np.eye(len(gate)), atol=ATOL, rtol=0)
                        np.testing.assert_allclose(np.linalg.det(gate), 1,
                                                   atol=ATOL, rtol=0)

    def test_weighted_map_stability_without_dividing_by_small_defects(self):
        for height in (1, 2, 3, 4):
            size = 1 << height
            native, generic, singular = _word_fixtures(height)
            identity = {node: np.eye(2, dtype=complex) for node in range(1, size)}
            for fixture, (coarse, target) in enumerate((
                    (native, generic), (native, singular), (generic, generic),
                    (identity, identity))):
                for perturbation in (1e-3, 1e-8):
                    with self.subTest(height=height, fixture=fixture,
                                      perturbation=perturbation):
                        changed = {node: target[node] @ _rotation(
                            perturbation * (0.5 + node / (2 * size)))
                            for node in range(1, size)}
                        delta = max(np.linalg.norm(changed[node] - target[node], ord=2)
                                    for node in range(1, size))
                        maps = []
                        for words in (target, changed):
                            blocks = _overlap_generators(height, coarse, words)
                            h = np.array([blocks[node][1, 0] for node in range(1, size)])
                            maps.append(h[:, None] * _path_matrix(height, words))
                        actual = np.linalg.norm(maps[0] - maps[1], ord=2)
                        self.assertGreater(delta, 0)
                        self.assertLessEqual(actual, 2 * height ** 1.5 * delta + ATOL)
                        if fixture == 3:
                            np.testing.assert_array_equal(maps[0], 0)
                            self.assertGreater(actual, 0)


if __name__ == "__main__":
    unittest.main()
