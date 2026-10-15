"""Bounded algebra audits for a bounded-density Haar correction.

The weighted frame is reconstructed independently from normalized subtree
vectors. Row and column source fixtures occupy one leaf register, one depth
tag, and disjoint failure flags. Their unitary completions are verification
fixtures, not supplied native source oracles. The complete small LCU word
retains its branch, factor, damping, and source work through actual inverses
and literal amplitude amplification. These checks establish no asymptotic
resource bound or scalable emitter and change no scientific receipt.
"""
from __future__ import annotations

from fractions import Fraction as F
import unittest

import numpy as np


ATOL = 4e-11


def _nodes(size):
    width = (size - 1).bit_length()
    for depth in range(width):
        length = size >> depth
        for prefix in range(1 << depth):
            marker = (2 * prefix + 1) << (width - depth - 1)
            yield marker, depth + 1, prefix * length, length


def _weighted_frame(weights):
    """Literal normalized complements, without the Haar-gauge identity."""
    weights = np.asarray(weights, dtype=float)
    size = len(weights)
    result = np.zeros((size, size))
    result[:, 0] = np.sqrt(weights / weights.sum())
    for marker, _, start, length in _nodes(size):
        middle, stop = start + length // 2, start + length
        left, right = weights[start:middle], weights[middle:stop]
        cosine = np.sqrt(left.sum() / (left.sum() + right.sum()))
        sine = np.sqrt(right.sum() / (left.sum() + right.sum()))
        result[start:middle, marker] = -sine * np.sqrt(left / left.sum())
        result[middle:stop, marker] = cosine * np.sqrt(right / right.sum())
    return result


def _gauge(weights):
    weights = np.asarray(weights, dtype=float)
    size = len(weights)
    haar = np.zeros((size, size))
    scaling = np.zeros_like(haar)
    c, a = np.ones(size), np.zeros(size)
    haar[:, 0] = 1 / np.sqrt(size)
    for marker, _, start, length in _nodes(size):
        middle, stop = start + length // 2, start + length
        left = weights[start:middle].mean()
        right = weights[middle:stop].mean()
        mean = (left + right) / 2
        haar[start:middle, marker] = -1 / np.sqrt(length)
        haar[middle:stop, marker] = 1 / np.sqrt(length)
        scaling[start:stop, marker] = 1 / np.sqrt(length)
        c[marker] = np.sqrt(mean / (left * right))
        a[marker] = (right - left) / (2 * np.sqrt(mean * left * right))
    return haar, scaling, c, a


def _reference_gauge(reference, density):
    """General positive reference, with its literal prescribed frame."""
    reference = np.asarray(reference, dtype=float)
    density = np.asarray(density, dtype=float)
    size = len(reference)
    frame = _weighted_frame(reference)
    scaling = np.zeros((size, size))
    c, correction = np.ones(size), np.zeros(size)
    for marker, _, start, length in _nodes(size):
        middle, stop = start + length // 2, start + length
        left_mass = reference[start:middle].sum()
        right_mass = reference[middle:stop].sum()
        mass = left_mass + right_mass
        left = np.dot(reference[start:middle], density[start:middle]) / left_mass
        right = np.dot(reference[middle:stop], density[middle:stop]) / right_mass
        mean = (left_mass * left + right_mass * right) / mass
        cosine, sine = np.sqrt(left_mass / mass), np.sqrt(right_mass / mass)
        scaling[start:stop, marker] = np.sqrt(reference[start:stop] / mass)
        c[marker] = np.sqrt(mean / (left * right))
        correction[marker] = cosine * sine * (right - left) / np.sqrt(mean * left * right)
    return frame, scaling, c, correction


def _source_index(size, leaf, tag, failure):
    return (4 * tag + failure) * size + leaf


def _schur_sources(matrix, alpha):
    """Two normalized isometries with a single shared leaf/depth label."""
    size = len(matrix)
    width = (size - 1).bit_length()
    labels = 1 << width.bit_length()
    doubled = np.block([[np.zeros_like(matrix), matrix],
                        [matrix.T, np.zeros_like(matrix)]])
    weights = np.linalg.solve(np.eye(2 * size) - doubled / alpha,
                              np.ones(2 * size))
    u, v = weights[:size], weights[size:]
    row = np.zeros((size * labels * 4, size), dtype=complex)
    column = np.zeros_like(row)
    # Column packing is j -> (anchor(v), t(v)). The state column has its
    # own tag 0 and routes entirely to column failure.
    column[_source_index(size, 0, 0, 1), 0] = 1
    for marker, tag, start, length in _nodes(size):
        for leaf in range(start, start + length):
            column[_source_index(size, leaf, tag, 0), marker] = np.sqrt(
                matrix[leaf, marker] * u[leaf] / (alpha * v[marker]))
        column[_source_index(size, start, tag, 1), marker] = 1 / np.sqrt(v[marker])
    # Independently recover the unique ancestor marker from leaf and depth;
    # no second full-size input index is retained in this common space.
    for leaf in range(size):
        for tag in range(1, width + 1):
            prefix = leaf >> (width - tag + 1)
            marker = (2 * prefix + 1) << (width - tag)
            row[_source_index(size, leaf, tag, 0), leaf] = np.sqrt(
                matrix[leaf, marker] * v[marker] / (alpha * u[leaf]))
        row[_source_index(size, leaf, 0, 2), leaf] = 1 / np.sqrt(u[leaf])
    return row, column, u, v


def _completion(isometry):
    """Actual unitary completion with its first input columns prescribed."""
    dimension, count = isometry.shape
    basis, triangular = np.linalg.qr(np.column_stack((isometry, np.eye(dimension))))
    for column in range(count):
        diagonal = triangular[column, column]
        basis[:, column] *= diagonal / abs(diagonal)
    return basis


def _axis_unitary(matrix, state, axis):
    moved = np.moveaxis(state, axis, 0)
    return np.moveaxis(np.tensordot(matrix, moved, axes=(1, 0)), 0, axis)


def _flag_rotation(state, ratios, axis, *, inverse=False):
    """Real attenuation rotation; the data axis is penultimate."""
    moved = np.moveaxis(state, axis, 0)
    ratio = np.asarray(ratios).reshape((1,) * (moved.ndim - 3) + (-1, 1))
    sine = np.sqrt(np.maximum(0, 1 - ratio * ratio)) * (-1 if inverse else 1)
    result = np.empty_like(moved)
    result[0] = ratio * moved[0] - sine * moved[1]
    result[1] = sine * moved[0] + ratio * moved[1]
    return np.moveaxis(result, 0, axis)


def _haar_contraction(haar, ratios, source_dimension):
    """Actual source-space unitary: attenuation on F0, then Haar on data."""
    size = len(haar)
    result = np.eye(source_dimension, dtype=complex)
    work = source_dimension // size
    for tag in range(work // 4):
        for failure in (0, 2):
            for leaf, ratio in enumerate(ratios):
                low = _source_index(size, leaf, tag, failure)
                high = _source_index(size, leaf, tag, failure | 1)
                pair = np.stack((result[low].copy(), result[high].copy()))
                sine = np.sqrt(max(0, 1 - ratio * ratio))
                result[low], result[high] = np.array([[ratio, -sine],
                                                     [sine, ratio]]) @ pair
    return np.einsum('ab,wbk->wak', haar, result.reshape(work, size, source_dimension)
                     ).reshape(source_dimension, source_dimension)


def _lcu_fixture(weights, *, perturbation=0.0):
    weights = np.asarray(weights, dtype=float)
    size = len(weights)
    haar, scaling, c, a = _gauge(weights)
    lower, upper = weights.min(), weights.max()
    correction_bound = 2 * (upper - lower) / lower ** 1.5
    alpha = max(1, 2 * correction_bound)
    c_bound, d_bound = 1 / np.sqrt(lower), np.sqrt(upper)
    sources = []
    for sign in (1, -1):
        matrix = scaling * np.maximum(sign * a, 0)[None, :]
        row, column, _, _ = _schur_sources(matrix, alpha)
        left, right = _completion(row), _completion(column)
        if sign == 1 and perturbation:
            # Mix a column-failure coordinate into an actual main-overlap
            # coordinate. Use the perturbed word's own inverse below.
            i = _source_index(size, 0, 0, 1)
            j = _source_index(size, 0, 1, 0)
            pair = np.stack((right[i].copy(), right[j].copy()))
            cosine, sine = np.cos(perturbation), np.sin(perturbation)
            right[i], right[j] = np.array([[cosine, -sine],
                                           [sine, cosine]]) @ pair
        sources.append(left.conj().T @ right)
    dimension = len(sources[0])
    branches = [_haar_contraction(haar, c / c_bound, dimension),
                -sources[0], sources[1], np.eye(dimension)]
    coefficients = np.array([c_bound, alpha, alpha, 0.0])
    selector = _completion(np.sqrt(coefficients / coefficients.sum())[:, None])
    normalization = d_bound * coefficients.sum()
    rounds = 0
    while np.sin(np.pi / (4 * rounds + 2)) > 1 / normalization + 1e-15:
        rounds += 1
    amplitude = np.sin(np.pi / (4 * rounds + 2))
    damping = normalization * amplitude
    return branches, selector, np.sqrt(weights) / d_bound, damping, rounds, amplitude


def _select(state, branches, *, inverse=False):
    result = np.empty_like(state)
    for selector, branch in enumerate(branches):
        old = state[:, :, selector]
        matrix = branch.conj().T if inverse else branch
        block = old.reshape(4, len(branch), state.shape[-1])
        result[:, :, selector] = np.einsum('ab,kbc->kac', matrix, block).reshape(old.shape)
    return result


def _lcu_word(state, fixture, *, inverse=False, active=True):
    branches, selector, diagonal, damping, _, _ = fixture
    if inverse and active:
        state = _flag_rotation(state, np.full(len(diagonal), damping), 0, inverse=True)
        state = _flag_rotation(state, diagonal, 1, inverse=True)
    state = _axis_unitary(selector, state, 2)
    if active:
        state = _select(state, branches, inverse=inverse)
    state = _axis_unitary(selector.conj().T, state, 2)
    if not inverse and active:
        # Outer D and damping each have a fresh flag; neither can convert
        # source or LCU failure into the all-zero accepted work subspace.
        state = _flag_rotation(state, diagonal, 1)
        state = _flag_rotation(state, np.full(len(diagonal), damping), 0)
    return state


def _zero_reflection(state, *, active=True):
    result = state.copy()
    if active:
        result[0, 0, 0, 0] *= -1
    return result


def _amplified(state, fixture, *, active=True):
    state = _lcu_word(state, fixture, active=active)
    for _ in range(fixture[4]):
        state = _zero_reflection(state, active=active)
        state = _lcu_word(state, fixture, inverse=True, active=active)
        state = _zero_reflection(state, active=active)
        state = _lcu_word(state, fixture, active=active)
        if active:
            state = -state  # Literal Grover sign, guarded by activity.
    return state


def _initial_columns(size, source_dimension):
    result = np.zeros((2, 2, 4, source_dimension // size, size, size), dtype=complex)
    result[0, 0, 0, 0] = np.eye(size)
    return result


class WeightedHaarDepthTests(unittest.TestCase):
    def test_weighted_frame_gauge_with_zero_correction_columns(self):
        for weights in ([1, 1], [1, 1, 1, 1], [0.9, 1.1, 0.9, 1.1],
                        [0.25, 1.75, 0.5, 1.5],
                        [0.5, 1.5, 0.75, 1.25, 1.75, 0.25, 1.25, 0.75]):
            with self.subTest(weights=weights):
                haar, scaling, c, a = _gauge(weights)
                frame = _weighted_frame(weights)
                actual = np.sqrt(weights)[:, None] * (haar * c - scaling * a)
                np.testing.assert_allclose(actual, frame, atol=ATOL)
                np.testing.assert_allclose(frame.T @ frame, np.eye(len(weights)), atol=ATOL)
                self.assertEqual(a[0], 0)
                self.assertEqual(c[0], 1)
        self.assertTrue(np.all(_gauge([1, 1, 1, 1])[3] == 0))
        self.assertEqual(_gauge([0.9, 1.1, 0.9, 1.1])[3][2], 0)

    def test_exact_subtree_martingale_energy_and_carleson_bound(self):
        weights = [F(1, 4), F(7, 4), F(1, 2), F(3, 2),
                   F(5, 4), F(3, 4), F(1), F(1)]
        lower = min(weights)
        nodes = list(_nodes(len(weights)))
        for _, _, start, length in nodes:
            mean = sum(weights[start:start + length]) / length
            variance = sum((value - mean) ** 2 for value in weights[start:start + length])
            energy, carleson = F(0), F(0)
            for _, _, child_start, child_length in nodes:
                if not (start <= child_start and child_start + child_length <= start + length):
                    continue
                middle = child_start + child_length // 2
                left = sum(weights[child_start:middle]) / (child_length // 2)
                right = sum(weights[middle:child_start + child_length]) / (child_length // 2)
                delta_squared = (right - left) ** 2 / 4
                energy += child_length * delta_squared
                carleson += child_length * delta_squared / (((left + right) / 2) * left * right)
            self.assertEqual(energy, variance)
            self.assertLessEqual(carleson, variance / lower ** 3)

    def test_single_leaf_depth_sources_and_distinct_failure_padding(self):
        for weights in ([1, 1, 1, 1], [0.9, 1.1, 0.9, 1.1], [0.94, 1.02, 1.07, 0.97]):
            _, scaling, _, a = _gauge(weights)
            for sign in (1, -1):
                matrix = scaling * np.maximum(sign * a, 0)[None, :]
                alpha = 1 + 2 * np.linalg.norm(matrix, 2)
                row, column, u, v = _schur_sources(matrix, alpha)
                np.testing.assert_allclose(row.conj().T @ row, np.eye(4), atol=ATOL)
                np.testing.assert_allclose(column.conj().T @ column, np.eye(4), atol=ATOL)
                np.testing.assert_allclose(row.conj().T @ column, matrix / alpha, atol=ATOL)
                np.testing.assert_allclose(matrix @ v, alpha * (u - 1), atol=ATOL)
                np.testing.assert_allclose(matrix.T @ u, alpha * (v - 1), atol=ATOL)
                self.assertTrue(np.all(u >= 1 - ATOL) and np.all(v >= 1 - ATOL))
                # In particular the state-input column is a full normalized
                # failure column, not a missing or postselected zero vector.
                self.assertEqual(column[_source_index(4, 0, 0, 1), 0], 1)

    def test_biased_product_reference_gauge_weighted_energy_and_overlap(self):
        # q=(1/5,4/5) tensor (4/5,1/5), so its leaf weights differ by 16.
        # The native approximation of this product core is a separate
        # obligation; this fixture audits its exact reference geometry.
        reference = [F(4, 25), F(1, 25), F(16, 25), F(4, 25)]
        raw_density = [F(3, 4), F(5, 4), F(6, 5), F(4, 5)]
        normalization = sum(q * value for q, value in zip(reference, raw_density))
        density = [value / normalization for value in raw_density]
        self.assertEqual(sum(q * value for q, value in zip(reference, density)), 1)
        self.assertEqual(max(reference) / min(reference), 16)
        q, w = np.asarray(reference, dtype=float), np.asarray(density, dtype=float)
        core, scaling, c, a = _reference_gauge(q, w)
        target = _weighted_frame(q * w)
        np.testing.assert_allclose(np.sqrt(w)[:, None] * (core * c - scaling * a),
                                   target, atol=ATOL)
        np.testing.assert_allclose(target.T @ target, np.eye(4), atol=ATOL)
        self.assertTrue(np.any(a > 0) and np.any(a < 0))
        np.testing.assert_array_equal(_reference_gauge(q, np.ones(4))[3], 0)
        for sign in (1, -1):
            matrix = scaling * np.maximum(sign * a, 0)[None, :]
            alpha = 1 + 2 * np.linalg.norm(matrix, 2)
            row, column, _, _ = _schur_sources(matrix, alpha)
            np.testing.assert_allclose(row.conj().T @ row, np.eye(4), atol=ATOL)
            np.testing.assert_allclose(column.conj().T @ column, np.eye(4), atol=ATOL)
            np.testing.assert_allclose(row.conj().T @ column, matrix / alpha, atol=ATOL)
        nodes = list(_nodes(4))
        for _, _, start, length in nodes:
            mass = sum(reference[start:start + length])
            mean = sum(reference[x] * density[x] for x in range(start, start + length)) / mass
            variance = sum(reference[x] * (density[x] - mean) ** 2
                           for x in range(start, start + length))
            energy, carleson = F(0), F(0)
            for _, _, child_start, child_length in nodes:
                if not (start <= child_start and child_start + child_length <= start + length):
                    continue
                middle, stop = child_start + child_length // 2, child_start + child_length
                left_mass, right_mass = sum(reference[child_start:middle]), sum(reference[middle:stop])
                child_mass = left_mass + right_mass
                left = sum(reference[x] * density[x] for x in range(child_start, middle)) / left_mass
                right = sum(reference[x] * density[x] for x in range(middle, stop)) / right_mass
                child_mean = (left_mass * left + right_mass * right) / child_mass
                increment = left_mass * right_mass * (right - left) ** 2 / child_mass
                energy += increment
                carleson += increment / (child_mean * left * right)
            self.assertEqual(energy, variance)
            self.assertLessEqual(carleson, variance / min(density) ** 3)

    def test_complete_lcu_amplification_literal_frame_and_all_work_return(self):
        for weights in ([1, 1, 1, 1], [0.94, 1.02, 1.07, 0.97]):
            fixture = _lcu_fixture(weights)
            self.assertLessEqual(fixture[4], 2)  # Keep this fixture bounded.
            initial = _initial_columns(4, len(fixture[0][0]))
            frame = _weighted_frame(weights)
            raw = _lcu_word(initial, fixture)
            np.testing.assert_allclose(raw[0, 0, 0, 0], fixture[5] * frame, atol=ATOL)
            expected = np.zeros_like(initial)
            expected[0, 0, 0, 0] = frame
            actual = _amplified(initial, fixture)
            np.testing.assert_allclose(actual, expected, atol=ATOL)
            # Arbitrary occupied work and complex reference amplitudes on
            # the inactive branch see the actual prep/inverse cancellation.
            rng = np.random.default_rng(8091)
            occupied = rng.normal(size=initial.shape) + 1j * rng.normal(size=initial.shape)
            occupied /= np.linalg.norm(occupied)
            np.testing.assert_allclose(_amplified(occupied, fixture, active=False),
                                       occupied, atol=ATOL)

    def test_perturbed_source_uses_actual_inverse_and_retains_leakage(self):
        weights = [0.94, 1.02, 1.07, 0.97]
        exact = _lcu_fixture(weights)
        epsilon = 1e-3
        perturbed = _lcu_fixture(weights, perturbation=epsilon)
        initial = _initial_columns(4, len(exact[0][0]))
        ideal = _amplified(initial, exact)
        actual = _amplified(initial, perturbed)
        difference = np.linalg.norm((actual - ideal).reshape(-1, 4), 2)
        word_bound = (2 * exact[4] + 1) * 2 * np.sin(epsilon / 2)
        self.assertLessEqual(difference, word_bound + ATOL)
        leakage = actual.copy()
        leakage[0, 0, 0, 0] = 0
        self.assertGreater(np.linalg.norm(leakage), 1e-8)
        np.testing.assert_allclose(actual.reshape(-1, 4).conj().T @ actual.reshape(-1, 4),
                                   np.eye(4), atol=ATOL)


if __name__ == '__main__':
    unittest.main()
