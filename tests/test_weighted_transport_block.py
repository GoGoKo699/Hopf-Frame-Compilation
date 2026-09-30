"""Small checks of weighted transport and candidate block completions.

Weighted operators have dimension at most 16 and complete dilations at most
32; permutation checks enumerate at most 512 labels. These checks distinguish
a bounded operator norm from orthogonal stopping columns, and verify logical
one-flag dilations and their level packing. They do not provide an elementary
Clifford+T implementation or an asymptotic gate bound.
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


def _swap_bits(value, first, second):
    if ((value >> first) ^ (value >> second)) & 1:
        value ^= (1 << first) | (1 << second)
    return value


def _heap_marker_program(height):
    """Literal guarded suffix reversals, followed by global index reversal."""
    size = 1 << height
    result = []
    for original in range(2 * size):
        value = original
        for leading_zeroes in range(height):
            leading_one = height - leading_zeroes - 1
            # The signal bit is excluded from the guard. Each guarded suffix
            # reversal preserves every leading-one prefix, including its own.
            if (value & (size - 1)) >> leading_one == 1:
                for low in range(leading_one // 2):
                    value = _swap_bits(value, low, leading_one - low - 1)
        for low in range(height // 2):
            value = _swap_bits(value, low, height - low - 1)
        result.append(value)
    return np.array(result)


def _level_packing_program(height, depth):
    """Execute the controlled cyclic shift and a fixed three-bit permutation."""
    size = 1 << height
    result = []
    # Continuation, marker, left child, right child codes become 00, 01, 10, 11.
    # The unused codes complete the same literal permutation for every level.
    three_bit_permutation = (4, 1, 5, 6, 7, 0, 2, 3)
    for original in range(2 * size):
        value = original
        if depth == 0:
            if value & 2:
                value ^= size
        else:
            # A=1 selects the child slots. Move the last child bit b to the
            # front of x,b, giving b,x and leaving the parent slots untouched.
            if value & (1 << (depth + 1)):
                for low in range(depth):
                    value = _swap_bits(value, low, low + 1)
            positions = (height, depth + 1, depth)
            code = sum(((value >> position) & 1) << (2 - index)
                       for index, position in enumerate(positions))
            replacement = three_bit_permutation[code]
            for index, position in enumerate(positions):
                value &= ~(1 << position)
                value |= ((replacement >> (2 - index)) & 1) << position
        result.append(value)
    return np.array(result)


def _heap_gather_program(height):
    """Gather stop modes using a predicate flag flip and one three-cycle."""
    size = 1 << height
    result = []
    for original in range(2 * size):
        value = original
        if (value & (size - 1)) >= 2:
            value ^= size
        if (value & (size - 1)) >> 1 == 0:
            if value & 1:
                value ^= size
            if not value & size:
                value ^= 1
        result.append(value)
    return np.array(result)


def _batched_riccati_network(height, target, data):
    """Use packed disjoint level blocks and explicit physical permutations."""
    weights, _, _, alpha, rho, beta, cross = data
    size = 1 << height
    steps, local_gates = [], []

    def append_permutation(mapping):
        steps.append(np.eye(2 * size, dtype=complex)[np.argsort(mapping)])

    def append_gate(rows, gate):
        step = np.eye(2 * size, dtype=complex)
        step[np.ix_(rows, rows)] = gate
        steps.append(step)
        local_gates.append(gate)

    physical_map = _heap_marker_program(height)
    append_permutation(np.argsort(physical_map))
    accepted = np.sqrt(rho[1]) / alpha
    rejected = np.sqrt(1 - accepted ** 2)
    append_gate([0, size], np.array([
        [accepted, -rejected], [rejected, accepted],
    ], dtype=complex))
    for depth in range(height - 1):
        packing = _level_packing_program(height, depth)
        append_permutation(packing)
        for node in range(1 << depth, 1 << (depth + 1)):
            word = target[node]
            children = np.sqrt([rho[2 * node], rho[2 * node + 1]])
            if rho[node] == 0:
                first = np.array([1, 0, 0, 0], dtype=complex)
            else:
                first = np.r_[weights[node], children * word[:, 0],
                              -cross[node].conjugate() / np.sqrt(beta[node])]
                first /= np.sqrt(rho[node])
            second = np.r_[0, children * word[:, 1], np.sqrt(beta[node])] / alpha
            gate = _complete_special_unitary(first, second)[[0, 3, 1, 2]]
            carrier = 0 if node == 1 else size + node
            rows = packing[[carrier, node, size + 2 * node, size + 2 * node + 1]]
            append_gate(rows, gate)
        append_permutation(np.argsort(packing))
    for node in range(size // 2, size):
        carrier = 0 if node == 1 else size + node
        phase = weights[node] / abs(weights[node]) if weights[node] != 0 else 1
        append_gate([carrier, node], np.diag([phase, phase.conjugate()]))
    append_permutation(_heap_gather_program(height))
    append_permutation(physical_map)

    network = np.eye(2 * size, dtype=complex)
    for step in steps:
        network = step @ network
    inverse = np.eye(2 * size, dtype=complex)
    for step in reversed(steps):
        inverse = step.conj().T @ inverse
    return network, inverse, local_gates


def _six_givens_elimination(unitary):
    """Fixed QR schedule with determinant-one two-level factors, including zeros."""
    reduced = unitary.copy()
    factors = []
    for column, first, second in ((0, 2, 3), (0, 1, 2), (0, 0, 1),
                                  (1, 2, 3), (1, 1, 2), (2, 2, 3)):
        a, b = reduced[first, column], reduced[second, column]
        norm = np.hypot(abs(a), abs(b))
        if norm == 0:
            gate = np.eye(2, dtype=complex)
        else:
            gate = np.array([[a.conjugate(), b.conjugate()], [-b, a]]) / norm
        rows = [first, second]
        reduced[rows] = gate @ reduced[rows]
        factor = np.eye(4, dtype=complex)
        factor[np.ix_(rows, rows)] = gate
        factors.append(factor)
    return reduced, factors


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

    def test_literal_packing_gather_and_physical_marker_permutations(self):
        for height in range(1, 9):
            size = 1 << height
            with self.subTest(height=height):
                physical = _heap_marker_program(height)
                expected = np.array([0] + [_marker(node, height)
                                           for node in range(1, size)])
                np.testing.assert_array_equal(physical, np.r_[expected, size + expected])
                np.testing.assert_array_equal(np.sort(physical), np.arange(2 * size))
                gather = _heap_gather_program(height)
                np.testing.assert_array_equal(np.sort(gather), np.arange(2 * size))
                self.assertEqual(gather[size + 1], 0)
                self.assertEqual(gather[0], 1)
                for node in range(2, size):
                    self.assertEqual(gather[size + node], node)
                    self.assertEqual(gather[node], size + node)
                for depth in range(height - 1):
                    packing = _level_packing_program(height, depth)
                    np.testing.assert_array_equal(np.sort(packing), np.arange(2 * size))
                    packed_rows = []
                    for address in range(1 << depth):
                        node = (1 << depth) + address
                        carrier = 0 if node == 1 else size + node
                        rows = [carrier, node, size + 2 * node, size + 2 * node + 1]
                        expected_rows = ([0, 1, 2, 3] if depth == 0 else
                                         [slot * (1 << depth) + address
                                          for slot in range(4)])
                        np.testing.assert_array_equal(packing[rows], expected_rows)
                        packed_rows.extend(packing[rows])
                    self.assertEqual(len(packed_rows), len(set(packed_rows)))

    def test_batched_one_flag_network_preserves_the_accepted_block_and_actual_inverse(self):
        for height in (1, 2, 3, 4):
            size = 1 << height
            for fixture, (coarse, target) in enumerate(_dilation_fixtures(height)):
                with self.subTest(height=height, fixture=fixture):
                    data = _riccati_data(height, coarse, target)
                    weights, _, _, alpha, _, _, _ = data
                    network, inverse, gates = _batched_riccati_network(height, target, data)
                    weighted = np.array([weights[node] for node in range(1, size)])
                    weighted = weighted[:, None] * _path_matrix(height, target)
                    np.testing.assert_allclose(
                        network[:size, :size], _embed_rows(height, weighted) / alpha,
                        atol=ATOL, rtol=0)
                    for product in (network.conj().T @ network, inverse @ network,
                                    network @ inverse):
                        np.testing.assert_allclose(product, np.eye(2 * size),
                                                   atol=ATOL, rtol=0)
                    np.testing.assert_allclose(inverse, network.conj().T,
                                               atol=ATOL, rtol=0)
                    for gate in gates:
                        np.testing.assert_allclose(gate.conj().T @ gate,
                                                   np.eye(len(gate)), atol=ATOL, rtol=0)
                        np.testing.assert_allclose(np.linalg.det(gate), 1,
                                                   atol=ATOL, rtol=0)

    def test_six_fixed_givens_factors_reconstruct_special_unitary_level_blocks(self):
        fixtures = [np.eye(4, dtype=complex), np.eye(4, dtype=complex)[[1, 2, 0, 3]],
                    np.diag([1j, -1j, -1, -1])]
        for height in (2, 3, 4):
            for coarse, target in _dilation_fixtures(height):
                data = _riccati_data(height, coarse, target)
                _, _, gates = _batched_riccati_network(height, target, data)
                fixtures.extend(gate for gate in gates if len(gate) == 4)
        for fixture, unitary in enumerate(fixtures):
            with self.subTest(fixture=fixture):
                reduced, factors = _six_givens_elimination(unitary)
                np.testing.assert_allclose(reduced, np.eye(4), atol=ATOL, rtol=0)
                reconstructed = np.eye(4, dtype=complex)
                for factor in factors:
                    np.testing.assert_allclose(np.linalg.det(factor), 1, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(factor.conj().T @ factor, np.eye(4),
                                               atol=ATOL, rtol=0)
                    reconstructed = reconstructed @ factor.conj().T
                np.testing.assert_allclose(reconstructed, unitary, atol=ATOL, rtol=0)

    def test_selected_completion_can_stay_far_while_accepted_block_tends_to_zero(self):
        # This concerns these particular completion conventions. It is not a
        # lower bound on compilers free to choose a different rejected block.
        alpha = 1 / 16
        coarse = {1: np.eye(2, dtype=complex)}

        def fixed_alpha_data(target):
            data = _riccati_data(1, coarse, target)
            weight = data[0][1]
            return (data[0], data[1], data[2], alpha,
                    {1: abs(weight) ** 2, 2: 0.0, 3: 0.0},
                    {1: alpha ** 2}, {1: 0j})

        for builder in (_riccati_network, _batched_riccati_network):
            zero, _, _ = builder(1, coarse, fixed_alpha_data(coarse))
            for angle in (-1e-3, -1e-5):
                with self.subTest(builder=builder.__name__, angle=angle):
                    target = {1: _rotation(angle)}
                    network, _, _ = builder(1, target, fixed_alpha_data(target))
                    difference = network - zero
                    np.testing.assert_allclose(np.linalg.norm(difference, ord=2), 2,
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(
                        np.linalg.norm(difference[:2, :2], ord=2),
                        abs(np.sin(angle)) / alpha, atol=ATOL, rtol=0)


if __name__ == "__main__":
    unittest.main()
