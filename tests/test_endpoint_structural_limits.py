"""Finite regression evidence for the endpoint structural-limit arguments.

These fixtures check all matrix columns and preserve literal scalar phases.
They are not dimension-independent proofs, an unrestricted T-count lower
bound, or a compiler for arbitrary frames. The Haar native words are bounded
two/three-qubit examples; the source counterexamples are mathematical matrix
checks, not implementations of their arbitrary-angle phase operations.
"""
from __future__ import annotations

from itertools import product
import unittest

import numpy as np

from compiler_robust_hopf.conventions import marker_label
from compiler_robust_hopf.frames import hopf_ry, real_frame_matrix, real_tree_data
from compiler_robust_hopf.native_coarse_fixture import (
    ELEMENTARY_GATES,
    controlled_reflection_word,
    toffoli_word,
)
from tests.test_operator_source_compiler import _apply_native_word


ATOL = 3e-12
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)


def _native_columns(width, word, columns):
    """Use the existing independent simulator; its alphabet spells Z as SS."""
    flat = []
    for name, wires in word:
        flat.extend([("S", wires[0])] * 2 if name == "Z" else [(name, *wires)])
    return _apply_native_word(width, flat, columns)


def _haar_core_word(n):
    """Bounded native C: early angles pi/4, final angles zero, n<=3."""
    word = ()
    for target in range(n - 1, 0, -1):
        suffix_zero = tuple((wire, 0) for wire in range(target))
        # Chronological Z then H gives HZ=Ry(pi/4).
        word += controlled_reflection_word(None, target, suffix_zero)
        word += controlled_reflection_word(0, target, suffix_zero)
    return word


def _nested_projectors(theta, n):
    """Local root columns, computed independently of the global frame."""
    dimension = 1 << n
    result = {}
    for depth in range(n + 1):
        height, size = n - depth, 1 << (n - depth)
        for prefix in range(1 << depth):
            anchor = prefix * size
            local_angles = [theta[(1 << (depth + level))
                                  + (prefix << level) + j - 1]
                            for level in range(height) for j in range(1 << level)]
            local_root = (real_frame_matrix(local_angles)[:, 0]
                          if height else np.ones(1))
            state = np.zeros(dimension)
            state[anchor:anchor + size] = local_root
            coordinate = np.zeros((dimension, dimension))
            coordinate[anchor:anchor + size, anchor:anchor + size] = np.eye(size)
            marker_space = coordinate.copy()
            marker_space[anchor, anchor] = 0
            result[(depth, prefix)] = (
                marker_space, coordinate - np.outer(state, state), coordinate, state)
    return result


class EndpointStructuralLimitsTests(unittest.TestCase):
    def test_nested_projectors_determine_all_columns_in_singular_charts(self):
        for n in (2, 3, 4):
            dimension = 1 << n
            unequal = np.linspace(-1.1, 2.3, dimension - 1)
            singular = unequal.copy()
            singular[0], singular[2] = 0, np.pi / 2
            for theta in (unequal, singular, np.zeros(dimension - 1)):
                with self.subTest(n=n, chart=tuple(theta)):
                    frame = real_frame_matrix(theta)
                    projectors = _nested_projectors(theta, n)
                    for (depth, prefix), (e, p, coordinate, state) in projectors.items():
                        np.testing.assert_allclose(frame @ e @ frame.T, p,
                                                   atol=ATOL, rtol=0)
                        self.assertAlmostEqual(np.linalg.norm(state), 1, delta=ATOL)
                        other = np.eye(dimension) - coordinate
                        for block in (other @ frame @ coordinate,
                                      coordinate @ frame @ other):
                            # Support/coupling claims concern all input columns.
                            self.assertLessEqual(np.linalg.svd(block, compute_uv=False)[1], ATOL)
                        if depth < n:
                            left = projectors[(depth + 1, 2 * prefix)]
                            right = projectors[(depth + 1, 2 * prefix + 1)]
                            marker = (2 * prefix + 1) << (n - depth - 1)
                            column = frame[:, marker]
                            np.testing.assert_allclose(p - left[1] - right[1],
                                                       np.outer(column, column),
                                                       atol=ATOL, rtol=0)
                    root = frame[:, 0]
                    np.testing.assert_allclose(np.eye(dimension) - projectors[(0, 0)][1],
                                               np.outer(root, root), atol=ATOL, rtol=0)
                    # Same prepared state, different completion: a marker swap
                    # must violate at least one individual subtree constraint.
                    changed = frame.copy()
                    changed[:, [1, 2]] = changed[:, [2, 1]]
                    np.testing.assert_allclose(changed[:, 0], root, atol=ATOL, rtol=0)
                    self.assertGreater(max(np.linalg.norm(changed @ e @ changed.T - p, 2)
                                           for e, p, _, _ in projectors.values()), 0.9)

    def test_nested_projector_stability_and_tree_offdiagonal_partition(self):
        for n in (2, 3, 4):
            dimension = 1 << n
            theta = np.linspace(-0.8, 1.9, dimension - 1)
            theta[0] = 0  # Zero global mass does not remove local constraints.
            frame = real_frame_matrix(theta)
            projectors = _nested_projectors(theta, n)
            indices = np.arange(dimension)
            raw = np.exp(1j * (indices[:, None] + 1) * (indices[None, :] + 2))
            hermitian = (raw + raw.conj().T) / (2 * dimension)
            values, vectors = np.linalg.eigh(hermitian)
            u = (vectors * np.exp(0.013j * values)) @ vectors.conj().T
            u = u @ np.diag(np.exp(0.37j * indices))
            candidate = frame @ u
            delta = max(np.linalg.norm(candidate @ e @ candidate.conj().T - p, 2)
                        for e, p, _, _ in projectors.values())
            kappa = (2 * n - 1) * delta
            with self.subTest(n=n):
                self.assertLess(kappa, 1)
                root_e = projectors[(0, 0)][0]
                root_zero = np.eye(dimension) - root_e
                parts = root_zero @ u @ root_e + root_e @ u @ root_zero
                self.assertLessEqual(np.linalg.norm(parts, 2), delta + ATOL)
                for depth in range(n - 1):
                    depth_part = np.zeros_like(u)
                    for prefix in range(1 << depth):
                        e = projectors[(depth, prefix)][0]
                        left = projectors[(depth + 1, 2 * prefix)][0]
                        right = projectors[(depth + 1, 2 * prefix + 1)][0]
                        marker = e - left - right
                        depth_part += (e @ u @ e - left @ u @ left
                                       - right @ u @ right - marker @ u @ marker)
                    self.assertLessEqual(np.linalg.norm(depth_part, 2), 2 * delta + ATOL)
                    parts += depth_part
                np.testing.assert_allclose(parts, u - np.diag(np.diag(u)),
                                           atol=ATOL, rtol=0)
                self.assertLessEqual(np.linalg.norm(parts, 2), kappa + ATOL)
                correction = np.diag(np.conj(np.diag(u)) / np.abs(np.diag(u)))
                self.assertLessEqual(np.linalg.norm(candidate @ correction - frame, 2),
                                     kappa + 1 - np.sqrt(1 - kappa ** 2) + ATOL)
                # Exact phase freedom satisfies the projectors and is removed
                # by an actual diagonal; no scalar quotient is used.
                phases = np.exp(0.21j * indices)
                exact = frame @ np.diag(phases)
                for e, p, _, _ in projectors.values():
                    np.testing.assert_allclose(exact @ e @ exact.conj().T, p,
                                               atol=ATOL, rtol=0)
                np.testing.assert_allclose(exact @ np.diag(phases.conj()), frame,
                                           atol=ATOL, rtol=0)

    def test_all_quarter_turn_edges_form_one_cycle_and_diagonal_transitions(self):
        for n in (2, 3, 4):
            dimension = 1 << n
            frame = real_frame_matrix(np.full(dimension - 1, np.pi / 2))
            permutation = np.argmax(np.abs(frame), axis=0)
            visited, vertex = set(), 0
            while vertex not in visited:
                visited.add(vertex)
                vertex = int(permutation[vertex])
            with self.subTest(n=n):
                self.assertEqual(vertex, 0)
                self.assertEqual(len(visited), dimension)
                for labels in (np.arange(dimension), np.arange(dimension) % 3,
                               np.repeat([0, 1], dimension // 2)):
                    diagonal = np.diag(labels)
                    displacement = diagonal @ frame - frame @ diagonal
                    transitions = np.count_nonzero(labels != labels[permutation])
                    self.assertEqual(np.linalg.matrix_rank(displacement, tol=ATOL), transitions)
                    self.assertGreaterEqual(transitions, len(set(labels)))

    def test_complex_pair_rotation_energy_identity_for_fixed_baselines(self):
        # The fixed-menu proof permits complex baseline entries. This fixture
        # checks the real-angle energy formula without a real-column assumption.
        pair = np.array([0.3 + 0.4j, -0.2 + 0.7j])
        rho = np.vdot(pair, pair).real
        a = (abs(pair[0]) ** 2 - abs(pair[1]) ** 2) / 2
        b = np.real(pair[0] * np.conj(pair[1]))
        self.assertLessEqual(np.hypot(a, b), rho / 2)
        for angle in np.linspace(-2.1, 2.7, 11):
            rotated = hopf_ry(angle).T @ pair
            wave = a * np.cos(2 * angle) + b * np.sin(2 * angle)
            np.testing.assert_allclose(np.abs(rotated) ** 2,
                                       [rho / 2 + wave, rho / 2 - wave],
                                       atol=ATOL, rtol=0)

    def test_full_frame_tangent_gram_in_unequal_and_singular_charts(self):
        for n in (2, 3, 4):
            count = (1 << n) - 1
            unequal = 0.37 * np.arange(1, count + 1) - 0.61
            singular = unequal.copy()
            singular[0] = 0.0
            singular[1] = np.pi / 2
            cases = (("unequal", unequal), ("singular", singular),
                     ("all_zero", np.zeros(count)),
                     ("all_quarter_turn", np.full(count, np.pi / 2)))
            for name, theta in cases:
                with self.subTest(n=n, chart=name):
                    frame = real_frame_matrix(theta)
                    tangents = []
                    for j in range(count):
                        shift = np.zeros(count)
                        shift[j] = np.pi / 2
                        # Exact parameter-shift identity, not a small-step
                        # finite difference. The angle occurs only once.
                        derivative = (real_frame_matrix(theta + shift)
                                      - real_frame_matrix(theta - shift)) / 2
                        tangents.append(frame.T @ derivative)
                    tangents = np.asarray(tangents)
                    np.testing.assert_allclose(
                        tangents + tangents.transpose(0, 2, 1), 0,
                        atol=ATOL, rtol=0)
                    gram = np.einsum("aij,bij->ab", tangents, tangents)
                    np.testing.assert_allclose(
                        gram, 2 * np.eye(count), atol=ATOL, rtol=0)
                    if name == "singular":
                        # The root's right child has a vanishing STATE
                        # derivative while its full-frame tangent survives.
                        state_data = real_tree_data(theta)
                        np.testing.assert_allclose(
                            state_data.derivatives[2], 0, atol=ATOL, rtol=0)
                        self.assertAlmostEqual(
                            np.linalg.norm(tangents[2]), np.sqrt(2), delta=ATOL)

    def test_marker_columns_recover_all_sines_and_cosines_without_chart_division(self):
        for n in (2, 3, 4):
            count = (1 << n) - 1
            cases = (np.linspace(-1.13, 2.71, count), np.zeros(count),
                     np.full(count, np.pi / 2))
            for case, theta in enumerate(cases):
                with self.subTest(n=n, case=case):
                    peeled = real_frame_matrix(theta).copy()
                    for depth in reversed(range(n)):
                        suffix = n - depth - 1
                        for prefix in range(1 << depth):
                            node = (1 << depth) + prefix
                            anchor = (2 * prefix) << suffix
                            marker = marker_label(node, n)
                            cosine = peeled[marker, marker]
                            sine = -peeled[anchor, marker]
                            np.testing.assert_allclose(
                                [cosine, sine],
                                [np.cos(theta[node - 1]), np.sin(theta[node - 1])],
                                atol=ATOL, rtol=0)
                            anchor_row = peeled[anchor].copy()
                            marker_row = peeled[marker].copy()
                            peeled[anchor] = cosine * anchor_row + sine * marker_row
                            peeled[marker] = -sine * anchor_row + cosine * marker_row
                    np.testing.assert_allclose(
                        peeled, np.eye(1 << n), atol=ATOL, rtol=0)

    def test_row_support_is_root_state_plus_ancestor_markers(self):
        for n in (2, 3, 5):
            dimension = 1 << n
            allowed = np.zeros((dimension, dimension), dtype=bool)
            allowed[:, 0] = True
            for leaf in range(dimension):
                for depth in range(n):
                    prefix = leaf >> (n - depth)
                    node = (1 << depth) + prefix
                    allowed[leaf, marker_label(node, n)] = True
            self.assertTrue(np.all(np.sum(allowed, axis=1) == n + 1))
            unequal = np.linspace(-0.83, 1.17, dimension - 1)
            singular = unequal.copy()
            singular[::3] = 0
            for theta in (unequal, singular):
                with self.subTest(n=n, singular=theta is singular):
                    frame = real_frame_matrix(theta)
                    np.testing.assert_allclose(
                        frame[~allowed], 0, atol=ATOL, rtol=0)
                    self.assertLessEqual(
                        np.linalg.norm(np.abs(frame), 2), np.sqrt(n + 1) + ATOL)

    def test_balanced_haar_saturates_the_absolute_norm_bound(self):
        for n in (1, 2, 3, 5):
            frame = real_frame_matrix(np.full((1 << n) - 1, np.pi / 4))
            with self.subTest(n=n):
                self.assertTrue(np.all(
                    np.count_nonzero(np.abs(frame) > ATOL, axis=1) == n + 1))
                self.assertAlmostEqual(
                    np.linalg.norm(np.abs(frame), 2), np.sqrt(n + 1), delta=ATOL)

    def test_residual_block_and_norms_with_unequal_final_angles(self):
        for n in (2, 3, 5):
            dimension = 1 << n
            half = dimension // 2
            coarse_angles = np.full(dimension - 1, np.pi / 4)
            # Unequal final angles cancel from C^T W; equality of those
            # angles is not an assumption of the residual formula.
            coarse_angles[half - 1:] = np.linspace(-0.43, 1.07, half)
            coarse = real_frame_matrix(coarse_angles)
            haar = real_frame_matrix(np.full(half - 1, np.pi / 4))
            permutation = list(range(0, dimension, 2)) + list(range(1, dimension, 2))
            for delta in (0.02, -0.017, 2.0 ** -12):
                with self.subTest(n=n, delta=delta):
                    target_angles = coarse_angles.copy()
                    target_angles[half - 1:] += delta
                    residual = coarse.T @ real_frame_matrix(target_angles) - np.eye(dimension)
                    cosine, sine = np.cos(delta), np.sin(delta)
                    block = np.block([
                        [(cosine - 1) * np.eye(half), -sine * haar.T],
                        [sine * haar, (cosine - 1) * np.eye(half)],
                    ])
                    np.testing.assert_allclose(
                        residual[np.ix_(permutation, permutation)], block,
                        atol=ATOL, rtol=0)
                    signed_norm = 2 * abs(np.sin(delta / 2))
                    absolute_norm = 1 - cosine + abs(sine) * np.sqrt(n)
                    self.assertAlmostEqual(
                        np.linalg.norm(residual, 2), signed_norm, delta=ATOL)
                    self.assertAlmostEqual(
                        np.linalg.norm(np.abs(residual), 2), absolute_norm, delta=ATOL)
                    # The ratio tends to sqrt(n); the finite-delta identity
                    # includes the small diagonal term rather than dropping it.
                    self.assertAlmostEqual(
                        absolute_norm / signed_norm,
                        abs(np.sin(delta / 2)) + abs(np.cos(delta / 2)) * np.sqrt(n),
                        delta=ATOL)

    def test_literal_controlled_h_and_rotation_have_two_t_gates(self):
        word = controlled_reflection_word(0, 0, ((1, 1),))
        controlled_h = _native_columns(2, word, np.eye(4, dtype=complex))
        expected_h = np.block([[I2, np.zeros((2, 2))], [np.zeros((2, 2)), H]])
        np.testing.assert_allclose(controlled_h, expected_h, atol=ATOL, rtol=0)
        self.assertEqual(sum(name in ("T", "TDG") for name, _ in word), 2)
        word = controlled_reflection_word(None, 0, ((1, 1),)) + word
        actual = _native_columns(2, word, np.eye(4, dtype=complex))
        expected = np.block([
            [I2, np.zeros((2, 2))], [np.zeros((2, 2)), hopf_ry(np.pi / 4)],
        ])
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        self.assertEqual(sum(name in ("T", "TDG") for name, _ in word), 2)

    def test_bounded_native_haar_core_matches_every_column(self):
        for n, expected_t, expected_length in ((2, 2, 20), (3, 18, 72)):
            with self.subTest(n=n):
                word = _haar_core_word(n)
                self.assertTrue(all(name in ELEMENTARY_GATES for name, _ in word))
                self.assertEqual(sum(name in ("T", "TDG") for name, _ in word), expected_t)
                self.assertEqual(len(word), expected_length)
                self.assertTrue(all(0 <= q < n for _, wires in word for q in wires))
                angles = np.zeros((1 << n) - 1)
                angles[:(1 << (n - 1)) - 1] = np.pi / 4
                actual = _native_columns(n, word, np.eye(1 << n, dtype=complex))
                np.testing.assert_allclose(
                    actual, real_frame_matrix(angles), atol=ATOL, rtol=0)

    def test_parallel_coin_boolean_parity_witness(self):
        for n in (2, 3, 4):
            signs = np.asarray(list(product((-1, 1), repeat=n)), dtype=int)
            parity = np.prod(signs, axis=1)
            amplitudes = []
            for row in signs:
                angles = np.zeros((1 << n) - 1)
                for depth, sign in enumerate(row):
                    # All-right path node has one-based heap label
                    # 2^(depth+1)-1.
                    angles[(1 << (depth + 1)) - 2] = sign * np.pi / 2
                amplitudes.append(real_frame_matrix(angles)[-1, 0])
            with self.subTest(n=n):
                np.testing.assert_allclose(amplitudes, parity, atol=ATOL, rtol=0)
                # Every Boolean monomial of degree <n is exactly orthogonal
                # to parity. This finite witness does not simulate arbitrary
                # oracle interlayers or prove the query-degree induction.
                for mask in range((1 << n) - 1):
                    columns = [q for q in range(n) if mask & (1 << q)]
                    monomial = np.prod(signs[:, columns], axis=1)
                    self.assertEqual(int(parity @ monomial), 0)

    def test_whole_program_dirty_echo_depends_on_original_dirty_word(self):
        # U(z1,z2)=H^z2 X^z1 and loaded word s=(1,0). The commutator
        # U(z xor s) U(z)^dagger gives X at z=00 but Z at z=01.
        def program(z1, z2):
            return (H if z2 else I2) @ (X if z1 else I2)

        clean_result = program(1, 0) @ program(0, 0).conj().T
        dirty_result = program(1, 1) @ program(0, 1).conj().T
        np.testing.assert_allclose(clean_result, X, atol=ATOL, rtol=0)
        np.testing.assert_allclose(dirty_result, Z, atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(dirty_result - X, 2), np.sqrt(2), delta=ATOL)

    def test_first_one_select_changes_a_coherent_dirty_tuple(self):
        # Tuple length four, followed by two logical data bits. Selecting the
        # first 1 applies Z to data bit 1 or 2, otherwise identity.
        select = np.eye(64, dtype=complex)
        for z in range(16):
            bits = [(z >> (3 - q)) & 1 for q in range(4)]
            first = bits.index(1) if 1 in bits else 4
            for data in range(4):
                if first < 2:
                    select[4 * z + data, 4 * z + data] = (-1) ** ((data >> (1 - first)) & 1)
        plus = np.ones(2) / np.sqrt(2)
        zero, one = np.eye(2)
        a = np.kron(np.kron(np.kron(one, plus), plus), plus)
        b = np.kron(np.kron(np.kron(zero, one), plus), plus)
        data = np.kron(one, zero)
        initial = np.kron((a + b) / np.sqrt(2), data)
        expected = np.kron((-a + b) / np.sqrt(2), data)
        np.testing.assert_allclose(select @ initial, expected, atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(select @ initial - initial), np.sqrt(2), delta=ATOL)
        uniform = np.kron(np.kron(np.kron(plus, plus), plus), plus)
        embedding = np.kron(uniform[:, None], np.eye(4))
        np.testing.assert_allclose(
            embedding.conj().T @ select @ embedding,
            0.5 * np.kron(Z, I2) + 0.25 * np.kron(I2, Z) + 0.25 * np.eye(4),
            atol=ATOL, rtol=0)

    def test_phase_halving_works_for_prepared_catalyst_not_arbitrary_input(self):
        theta = 0.37
        phase2 = np.diag([1, np.exp(2j * theta)])
        phase1 = np.diag([1, np.exp(1j * theta)])
        p00 = np.diag([1, 0, 0, 0])
        podd = np.diag([0, 1, 1, 0])
        p11 = np.diag([0, 0, 0, 1])
        effective = (np.kron(I2, p00) + np.kron(phase2 @ X, podd)
                     + np.exp(2j * theta) * np.kron(I2, p11))
        # Registers in integer order: catalyst c, data x,y, clean helper a.
        # The first two CNOTs are deliberately not uncomputed by this gadget.
        before = (("CX", (2, 3)), ("CX", (1, 3)),
                  ("CX", (3, 2)), ("CX", (3, 1)))
        before += toffoli_word(2, 1, 0) + (("CX", (3, 0)),)
        after = (("CX", (3, 0)),) + toffoli_word(2, 1, 0)
        after += (("CX", (3, 1)), ("CX", (3, 2)))
        embedding = np.eye(16, dtype=complex)[:, ::2]
        output = _native_columns(4, before, embedding)
        output[1::2] *= np.exp(2j * theta)
        output = _native_columns(4, after, output)
        np.testing.assert_allclose(output, embedding @ effective, atol=ATOL, rtol=0)
        ideal = np.kron(I2, np.kron(phase1, phase1))
        zero, one = np.eye(2)
        data = np.kron(zero, one)
        prepared = np.kron(np.array([1, np.exp(1j * theta)]) / np.sqrt(2), data)
        opposite = np.kron(np.array([1, -np.exp(1j * theta)]) / np.sqrt(2), data)
        dirty = np.kron(zero, data)
        np.testing.assert_allclose(effective @ prepared, ideal @ prepared, atol=ATOL, rtol=0)
        np.testing.assert_allclose(effective @ opposite, -ideal @ opposite, atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(effective @ dirty - ideal @ dirty), np.sqrt(2), delta=ATOL)

    def test_reused_majorana_flag_gives_second_chebyshev_moment(self):
        # Two arbitrary source qubits, two data qubits, and a most-significant
        # clean flag. Build the actual controlled-M / N / controlled-M word,
        # not a block matrix defined in terms of the desired first moment.
        y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        gammas = (np.kron(X, I2), np.kron(y, I2), np.kron(Z, X))
        data_involutions = (np.kron(Z, I2), np.kron(I2, Z), np.eye(4))
        probabilities = (0.5, 0.25, 0.25)
        source = sum(np.sqrt(p) * np.kron(gamma, np.eye(4))
                     for p, gamma in zip(probabilities, gammas, strict=True))
        masked_source = sum(np.sqrt(p) * np.kron(gamma, pauli)
                            for p, gamma, pauli in
                            zip(probabilities, gammas, data_involutions, strict=True))
        identity = np.eye(16, dtype=complex)
        zero = np.zeros_like(identity)
        np.testing.assert_allclose(source @ source, identity, atol=ATOL, rtol=0)
        np.testing.assert_allclose(masked_source @ masked_source, identity, atol=ATOL, rtol=0)
        control_zero = np.block([[source, zero], [zero, identity]])
        control_one = np.block([[identity, zero], [zero, source]])
        flag_h = np.kron(H, identity)
        actual = flag_h @ control_zero @ np.kron(I2, masked_source) @ control_one @ flag_h
        np.testing.assert_allclose(
            actual.conj().T @ actual, np.eye(32), atol=ATOL, rtol=0)
        h0 = sum(p * pauli for p, pauli in
                 zip(probabilities, data_involutions, strict=True))
        np.testing.assert_allclose(
            actual[:16, :16], np.kron(np.eye(4), h0), atol=ATOL, rtol=0)
        second_moment = (actual @ actual)[:16, :16]
        np.testing.assert_allclose(
            second_moment, np.kron(np.eye(4), 2 * h0 @ h0 - np.eye(4)),
            atol=ATOL, rtol=0)
        # Four columns span EVERY source input with data |10>. Here H0=0,
        # but two uses return the source with a literal minus sign.
        source_columns = np.kron(np.eye(4), np.eye(4)[:, 2:3])
        np.testing.assert_allclose(
            np.kron(np.eye(4), h0 @ h0) @ source_columns, 0, atol=ATOL, rtol=0)
        np.testing.assert_allclose(
            second_moment @ source_columns, -source_columns, atol=ATOL, rtol=0)
        initialized = np.vstack([source_columns, np.zeros_like(source_columns)])
        np.testing.assert_allclose(
            actual @ actual @ initialized, -initialized, atol=ATOL, rtol=0)


if __name__ == "__main__":
    unittest.main()
