"""Finite whole-space checks of the two-flag residual assembly.

Independent frame matrices and subtree Schur solves audit the affine forward
block, the swapped-pair actual inverse, and the final amplified isometry.
Circuit matrices have dimension at most 64. These are logical unitary checks,
not native emitters, asymptotic cost measurements, or large simulations.
Separate integer/rational checks cover the independent forest support and
finite packing inequalities; their lower bound concerns the enlarged family.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

import numpy as np

from tests.test_tree_residual_structure import (
    ATOL, _addressed_frame, _edge_product, _marker, _su2,
)
from tests.test_tree_transport import (
    _embed_rows, _overlap_generators, _path_matrix, _word_fixtures,
)
from tests.test_weighted_transport_block import (
    _batched_riccati_network, _complete_special_unitary, _heap_marker_program,
    _level_packing_program,
)


EPSILON = 1 / 64
ALPHA = 4 * EPSILON
SCALE = 2 - ALPHA


def _close_fixtures(height):
    size = 1 << height
    native, generic, singular = _word_fixtures(height)
    identity = {node: np.eye(2, dtype=complex) for node in range(1, size)}

    def changed(words, odd_only=False):
        result = {}
        for node in range(1, size):
            delta = EPSILON * (0.5 + node / size) / (8 * height)
            perturbation = (_su2(delta, 0.4 * delta, -0.2 * delta)
                            if not odd_only or node % 2 else np.eye(2))
            result[node] = words[node] @ perturbation
        return result

    phases = {node: _su2(0, EPSILON / (8 * height), 0)
              if node % 2 else np.eye(2, dtype=complex) for node in range(1, size)}
    return ((native, changed(native)), (singular, changed(singular, True)),
            (generic, generic), (identity, identity), (identity, phases))


def _affine_data(height, coarse, target):
    size = 1 << height
    blocks = _overlap_generators(height, coarse, target)
    rho = {leaf: 0.0 for leaf in range(size, 2 * size)}
    beta, cross = {}, {}
    for node in range(size - 1, 0, -1):
        h, diagonal = blocks[node][1]
        word = target[node]
        children = np.array([rho[2 * node], rho[2 * node + 1]])
        continuation = np.vdot(word[:, 0], children * word[:, 0]).real
        injection = np.vdot(word[:, 1], children * word[:, 1]).real
        cross[node] = np.conj(h) * diagonal + np.vdot(word[:, 0], children * word[:, 1])
        beta[node] = SCALE ** 2 - abs(diagonal) ** 2 - injection
        rho[node] = abs(h) ** 2 + continuation + abs(cross[node]) ** 2 / beta[node]
    return blocks, rho, beta, cross


def _fixed_weighted_data(height, coarse, target):
    """Use the assembly's fixed alpha, including exactly coincident frames."""
    size = 1 << height
    blocks = _overlap_generators(height, coarse, target)
    weights = {node: blocks[node][1, 0] for node in range(1, size)}
    defects = {node: 1 - abs(blocks[node][0, 0]) ** 2 for node in range(1, size)}
    rho = {leaf: 0.0 for leaf in range(size, 2 * size)}
    beta, cross = {}, {}
    for node in range(size - 1, 0, -1):
        word = target[node]
        children = np.array([rho[2 * node], rho[2 * node + 1]])
        continuation = np.vdot(word[:, 0], children * word[:, 0]).real
        injection = np.vdot(word[:, 1], children * word[:, 1]).real
        cross[node] = np.vdot(word[:, 0], children * word[:, 1])
        beta[node] = ALPHA ** 2 - injection
        rho[node] = (abs(weights[node]) ** 2 + continuation
                     + abs(cross[node]) ** 2 / beta[node])
    return weights, defects, max(0, max(defects.values())), ALPHA, rho, beta, cross


def _orthogonal_column(column):
    seed = np.eye(len(column), dtype=complex)[:, 0]
    result = seed - column * np.vdot(column, seed)
    return result / np.linalg.norm(result)


def _affine_packing(height, depth):
    if depth:
        return _level_packing_program(height, depth)
    # Affine root continuation is (1,1), so it uses the same three-bit
    # permutation as other levels, rather than the weighted root's CNOT.
    size = 1 << height
    permutation = (4, 1, 5, 6, 7, 0, 2, 3)
    result = []
    for original in range(2 * size):
        code = ((original >> height) << 2) | (original & 3)
        replacement = permutation[code]
        result.append((original & ~(size | 3)) | ((replacement >> 2) << height)
                      | (replacement & 3))
    return np.array(result)


def _affine_network(height, target, data):
    blocks, rho, beta, cross = data
    size = 1 << height
    steps, local_gates = [], []

    def append_gate(rows, gate):
        step = np.eye(2 * size, dtype=complex)
        step[np.ix_(rows, rows)] = gate
        steps.append(step)
        local_gates.append(gate)

    def append_permutation(mapping):
        steps.append(np.eye(2 * size, dtype=complex)[np.argsort(mapping)])

    physical_map = _heap_marker_program(height)
    append_permutation(np.argsort(physical_map))
    root_overlap = blocks[1][0, 0]
    first = np.array([root_overlap / SCALE,
                      np.sqrt(1 - (abs(root_overlap) ** 2 + rho[1]) / SCALE ** 2),
                      np.sqrt(rho[1]) / SCALE], dtype=complex)
    root = _complete_special_unitary(first, _orthogonal_column(first))
    append_gate([0, size, size + 1], root)
    for depth in range(height - 1):
        packing = _affine_packing(height, depth)
        append_permutation(packing)
        for node in range(1 << depth, 1 << (depth + 1)):
            word = target[node]
            h, diagonal = blocks[node][1]
            children = np.sqrt([rho[2 * node], rho[2 * node + 1]])
            second = np.r_[diagonal, children * word[:, 1], np.sqrt(beta[node])] / SCALE
            if rho[node] == 0:
                first = _orthogonal_column(second)
            else:
                first = np.r_[h, children * word[:, 0],
                              -cross[node].conjugate() / np.sqrt(beta[node])]
                first /= np.sqrt(rho[node])
            gate = _complete_special_unitary(first, second)[[0, 3, 1, 2]]
            rows = packing[[size + node, node, size + 2 * node, size + 2 * node + 1]]
            append_gate(rows, gate)
        append_permutation(np.argsort(packing))
    for node in range(size // 2, size):
        h, diagonal = blocks[node][1]
        phase = h / abs(h) if h != 0 else 1
        root_beta = np.sqrt(beta[node])
        gate = np.array([[phase * root_beta, diagonal],
                         [-diagonal.conjugate(), np.conj(phase) * root_beta]]) / SCALE
        append_gate([size + node, node], gate)
    gather = [value ^ size if value % size else value for value in range(2 * size)]
    append_permutation(gather)
    append_permutation(physical_map)
    network = np.eye(2 * size, dtype=complex)
    for step in steps:
        network = step @ network
    inverse = np.eye(2 * size, dtype=complex)
    for step in reversed(steps):
        inverse = step.conj().T @ inverse
    return network, inverse, local_gates


def _select_assembly(forward, forward_inverse, reverse, reverse_inverse):
    dimension = forward.shape[0]
    coefficients = np.array([[np.sqrt(SCALE / 2), -np.sqrt(ALPHA / 2)],
                             [np.sqrt(ALPHA / 2), np.sqrt(SCALE / 2)]])
    prep = np.kron(coefficients, np.eye(dimension))
    zero = np.zeros_like(forward)
    select = np.block([[forward, zero], [zero, reverse]])
    select_inverse = np.block([[forward_inverse, zero], [zero, reverse_inverse]])
    return prep.conj().T @ select @ prep, prep.conj().T @ select_inverse @ prep


def _assembly(height, coarse, target):
    forward, forward_inverse, gates = _affine_network(
        height, target, _affine_data(height, coarse, target))
    swapped, swapped_inverse, _ = _batched_riccati_network(
        height, coarse, _fixed_weighted_data(height, target, coarse))
    unitary, inverse = _select_assembly(forward, forward_inverse, swapped_inverse, swapped)
    size = 1 << height
    reflection = np.eye(4 * size, dtype=complex)
    reflection[:size, :size] *= -1
    amplified = -unitary @ reflection @ inverse @ reflection @ unitary
    return unitary, inverse, amplified, (forward, forward_inverse, swapped, swapped_inverse), gates


class ResidualAssemblyTests(unittest.TestCase):
    def test_independent_forest_support_and_positive_selected_walsh_columns(self):
        # Integer supports are compared with a separate marker-depth rule.
        # Raw Walsh signs suffice: each column's square-root scale cancels
        # exactly against its diagonal filter in the analytic parametrization.
        for height in range(3, 9):
            size = 1 << height
            columns, support = {0}, {(row, 0) for row in range(1, size)}
            for depth in range(height):
                span = 1 << (height - depth)
                for prefix in range(1 << depth):
                    anchor = prefix * span
                    marker = anchor + span // 2
                    self.assertNotIn(marker, columns)
                    columns.add(marker)
                    flipped_suffix = (marker % span) ^ (span // 2)
                    for row in range(anchor, anchor + span):
                        walsh_sign = (-1) ** ((row % span) & flipped_suffix).bit_count()
                        self.assertEqual(walsh_sign, 1)
                        if row not in (anchor, marker):
                            support.add((row, marker))
            independent_support = {(row, 0) for row in range(1, size)}
            for marker in range(1, size):
                span = 2 * (marker & -marker)
                anchor = (marker // span) * span
                for row in range(anchor, anchor + span):
                    if row and (row & -row) < span // 2:
                        independent_support.add((row, marker))
            self.assertEqual(columns, set(range(size)))
            self.assertEqual(support, independent_support)
            dimension = (height - 1) * size + 1
            self.assertEqual(len(support), dimension)
            # All cube corners obey the Frobenius and filter bounds, with
            # rational squares so no floating-point square root is needed.
            radius_squared = Fraction(1, 16 * dimension)
            self.assertEqual(dimension * radius_squared, Fraction(1, 16))
            self.assertLessEqual(size * radius_squared, Fraction(1, 16))

    def test_independent_forest_rational_grid_exceeds_native_word_capacity(self):
        # A rational inner cube of radius 1/(4d) is enough. Its grid has at
        # least 2**(N-ceil(log2(12d))) points per entry at eta=2**(-N).
        # Compare exponents only; no exponential-size matrix or grid exists.
        for height in (8, 16, 32, 64, 128):
            size = 1 << height
            dimension = (height - 1) * size + 1
            width = size + 2 * height + 9
            radius = Fraction(1, 4 * dimension)
            self.assertLessEqual(dimension * radius ** 2, Fraction(1, 16))
            self.assertLessEqual(size * radius ** 2, Fraction(1, 16))
            ceiling_log = (12 * dimension - 1).bit_length()
            grid_exponent = size - ceiling_log
            t_count = size * height // 8
            packing_exponent = dimension * grid_exponent
            word_exponent = (2 * width ** 2 + 3 * width + 5
                             + (2 * width + 1) * t_count)
            self.assertGreater(grid_exponent, 0)
            self.assertGreater(packing_exponent, word_exponent)
        eta = Fraction(1, 1 << 16)
        spacing = 3 * eta
        # Differing in a single coordinate is an exact rank-one matrix
        # separation; two eta-balls cannot cover both of those grid targets.
        self.assertGreater(spacing, 2 * eta)

    def test_affine_riccati_matches_independent_dense_subtree_schur_solve(self):
        for height in (1, 2, 3, 4):
            size = 1 << height
            for fixture, (coarse, target) in enumerate(_close_fixtures(height)):
                blocks, rho, beta, _ = _affine_data(height, coarse, target)
                with self.subTest(height=height, fixture=fixture):
                    self.assertLessEqual(np.linalg.norm(
                        _addressed_frame(height, coarse) - _addressed_frame(height, target),
                        ord=2), EPSILON)
                    for root in range(1, size):
                        descendants = [node for node in range(root, size)
                                       if _edge_product(target, root, node) is not None]
                        continuation = np.array([blocks[node][1, 0]
                                                 * _edge_product(target, root, node)
                                                 for node in descendants])
                        affine = np.zeros((len(descendants), len(descendants)), dtype=complex)
                        for row, node in enumerate(descendants):
                            for column, ancestor in enumerate(descendants):
                                if ancestor == node:
                                    affine[row, column] = blocks[node][1, 1]
                                else:
                                    path = _edge_product(target, ancestor, node, True)
                                    if path is not None:
                                        affine[row, column] = blocks[node][1, 0] * path
                        schur = SCALE ** 2 * np.eye(len(descendants)) - affine.conj().T @ affine
                        coupling = affine.conj().T @ continuation
                        independent = (np.vdot(continuation, continuation)
                                       + np.vdot(coupling, np.linalg.solve(schur, coupling)))
                        np.testing.assert_allclose(rho[root], independent, atol=ATOL, rtol=0)
                        self.assertGreater(np.linalg.eigvalsh(schur)[0], 0)
                        self.assertGreater(beta[root], 0)
                        defect = max(0, 1 - abs(blocks[root][0, 0]) ** 2)
                        self.assertLessEqual(rho[root], 2 * defect + ATOL)
                    self.assertLess(abs(blocks[1][0, 0]) ** 2 + rho[1], SCALE ** 2)

    def test_complete_affine_reverse_and_two_flag_amplified_frame(self):
        for height in (1, 2, 3, 4):
            size = 1 << height
            for fixture, (coarse, target) in enumerate(_close_fixtures(height)):
                with self.subTest(height=height, fixture=fixture):
                    unitary, inverse, amplified, components, gates = _assembly(height, coarse, target)
                    forward, forward_inverse, swapped, reverse = components
                    blocks = _overlap_generators(height, coarse, target)
                    diagonal = np.empty(size, dtype=complex)
                    diagonal[0] = blocks[1][0, 0]
                    for node in range(1, size):
                        diagonal[_marker(node, height)] = blocks[node][1, 1]
                    h = np.array([blocks[node][1, 0] for node in range(1, size)])
                    k = np.array([blocks[node][0, 1] for node in range(1, size)])
                    expected_forward = np.diag(diagonal) + _embed_rows(
                        height, h[:, None] * _path_matrix(height, target))
                    expected_reverse = _embed_rows(
                        height, k.conj()[:, None] * _path_matrix(height, coarse)).conj().T
                    np.testing.assert_allclose(forward[:size, :size], expected_forward / SCALE,
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(reverse[:size, :size], expected_reverse / ALPHA,
                                               atol=ATOL, rtol=0)
                    for circuit, actual_inverse in ((forward, forward_inverse), (swapped, reverse),
                                                     (unitary, inverse)):
                        np.testing.assert_allclose(actual_inverse, circuit.conj().T,
                                                   atol=ATOL, rtol=0)
                        np.testing.assert_allclose(actual_inverse @ circuit,
                                                   np.eye(len(circuit)), atol=ATOL, rtol=0)
                    for gate in gates:
                        np.testing.assert_allclose(gate.conj().T @ gate, np.eye(len(gate)),
                                                   atol=ATOL, rtol=0)
                        np.testing.assert_allclose(np.linalg.det(gate), 1, atol=ATOL, rtol=0)
                    c_frame = _addressed_frame(height, coarse)
                    w_frame = _addressed_frame(height, target)
                    relative_frame = c_frame.conj().T @ w_frame
                    np.testing.assert_allclose(unitary[:size, :size], relative_frame / 2,
                                               atol=ATOL, rtol=0)
                    embedding = np.eye(4 * size)[:, :size]
                    np.testing.assert_allclose(amplified @ embedding,
                                               embedding @ relative_frame, atol=ATOL, rtol=0)
                    complete = np.kron(np.eye(4), c_frame) @ amplified
                    np.testing.assert_allclose(complete @ embedding, embedding @ w_frame,
                                               atol=ATOL, rtol=0)

    def test_affine_mode_packing_and_gather_are_literal_permutations(self):
        for height in range(2, 8):
            size = 1 << height
            for depth in range(height - 1):
                packing = _affine_packing(height, depth)
                np.testing.assert_array_equal(np.sort(packing), np.arange(2 * size))
                for address in range(1 << depth):
                    node = (1 << depth) + address
                    rows = [size + node, node, size + 2 * node, size + 2 * node + 1]
                    expected = [slot * (1 << depth) + address for slot in range(4)]
                    np.testing.assert_array_equal(packing[rows], expected)
            gather = np.array([value ^ size if value % size else value
                               for value in range(2 * size)])
            np.testing.assert_array_equal(np.sort(gather), np.arange(2 * size))
            np.testing.assert_array_equal(gather[size + 1:], np.arange(1, size))
            self.assertEqual(gather[0], 0)
            self.assertEqual(gather[size], size)

    def test_actual_inverse_hybrid_controls_all_dirty_columns_and_rejected_leakage(self):
        height, size = 2, 4
        coarse, target = _close_fixtures(height)[0]
        unitary, _, _, components, _ = _assembly(height, coarse, target)
        forward, _, swapped, _ = components
        pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
        pauli_z = np.diag([1, -1])
        forward = np.kron(forward, np.eye(2))
        swapped = np.kron(swapped, np.eye(2))
        generator = np.kron(pauli_x, np.kron(np.eye(size), pauli_z))
        reverse_generator = np.kron(pauli_z, np.kron(np.eye(size), pauli_x))
        for delta in (1e-3, 1e-6):
            with self.subTest(delta=delta):
                noise = np.cos(delta) * np.eye(4 * size) + 1j * np.sin(delta) * generator
                reverse_noise = (np.cos(delta) * np.eye(4 * size)
                                 + 1j * np.sin(delta) * reverse_generator)
                forward_changed = noise @ forward
                swapped_changed = reverse_noise @ swapped
                changed, actual_inverse = _select_assembly(
                    forward_changed, forward.conj().T @ noise.conj().T,
                    swapped.conj().T @ reverse_noise.conj().T, swapped_changed)
                np.testing.assert_allclose(actual_inverse @ changed, np.eye(8 * size),
                                           atol=ATOL, rtol=0)
                reference = np.kron(unitary, np.eye(2))
                query_error = np.linalg.norm(changed - reference, ord=2)
                self.assertLessEqual(query_error, 2 * np.sin(delta / 2) + ATOL)
                reflection = np.eye(8 * size)
                reflection[:2 * size, :2 * size] *= -1
                amplified = -changed @ reflection @ actual_inverse @ reflection @ changed
                embedding = np.eye(8 * size)[:, :2 * size]
                c_frame = _addressed_frame(height, coarse)
                w_frame = _addressed_frame(height, target)
                complete = np.kron(np.eye(4), np.kron(c_frame, np.eye(2))) @ amplified
                expected = embedding @ np.kron(w_frame, np.eye(2))
                isometry_error = np.linalg.norm(complete @ embedding - expected, ord=2)
                self.assertLessEqual(isometry_error, 3 * query_error + ATOL)
                self.assertGreater(np.linalg.norm((complete @ embedding)[2 * size:], ord=2), 0)

    def test_wrong_reverse_or_relative_root_phase_cannot_pass_the_frame_contract(self):
        height, size = 3, 8
        coarse, target = _close_fixtures(height)[0]
        unitary, _, _, components, _ = _assembly(height, coarse, target)
        forward, forward_inverse, swapped, reverse = components
        expected = _addressed_frame(height, coarse).conj().T @ _addressed_frame(height, target)
        wrong_reverse, _ = _select_assembly(forward, forward_inverse, swapped, reverse)
        self.assertGreater(np.linalg.norm(wrong_reverse[:size, :size] - expected / 2, ord=2), 1e-5)
        altered_forward = forward.copy()
        altered_forward[0] *= -1  # A coherent phase on the root's accepted output.
        wrong_phase, _ = _select_assembly(
            altered_forward, altered_forward.conj().T, reverse, swapped)
        self.assertGreater(np.linalg.norm(wrong_phase[:size, :size] - expected / 2, ord=2), 0.5)
        reflection = np.eye(4 * size)
        reflection[:size, :size] *= -1
        wrong_amplified = -unitary @ reflection @ unitary @ reflection @ unitary
        embedding = np.eye(4 * size)[:, :size]
        self.assertGreater(np.linalg.norm(wrong_amplified @ embedding - embedding @ expected,
                                         ord=2), 1e-4)


if __name__ == '__main__':
    unittest.main()
