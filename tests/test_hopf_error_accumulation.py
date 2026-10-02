"""Bounded Hopf angle and actual shared-flag source-error diagnostics.

The source fixtures use the literal Q and actual-inverse OAA identities.
One common Clifford-algebra representation contains every source mask in a
fixture, including inactive rows. This is an exact algebraic reduction of
the source operators, evaluated in complex128; it is not a new emitted
Clifford+T gate word. Dense matrices have dimension at most 64. Larger
small-n witnesses propagate a single state through 32-dimensional blocks.
The uniform asymptotic statements are proved in HOPF_ERROR_ACCUMULATION.md.
"""
from __future__ import annotations

import math
import unittest

import numpy as np

from compiler_robust_hopf.frames import direct_real_frame, hopf_ry


ATOL = 4e-12
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
H = (X + Z) / math.sqrt(2)


def _mask_vector(value, amplitudes):
    """Nearest grid point, with the chapter's literal endpoint convention."""
    width = len(amplitudes)
    scale = 1 << (width - 2)
    index = round((1 - value) * scale)
    if index == 2 * scale:
        signs = [-1] * width
    else:
        signs = [(-1) ** ((index >> (width - 2 - j)) & 1)
                 for j in range(width - 1)] + [1]
    return amplitudes * signs, 1 - index / scale


def _common_sources(width):
    """Compress all five vectors jointly; never reorient one layer alone."""
    delta = 2.0 ** (2 - width)
    amplitudes = np.array([2.0 ** (-(j + 1) / 2) for j in range(width - 1)]
                          + [2.0 ** (-(width - 1) / 2)])
    nzero = amplitudes.copy()
    nzero[0] *= -1
    rows = [amplitudes, nzero]
    angles, coefficients = [], []
    for factor in (0.5, 1.5):
        sine = math.sqrt(factor * delta)
        theta = math.asin(sine)
        nc, cosine = _mask_vector(math.cos(theta), amplitudes)
        ns, sine_rounded = _mask_vector(sine, amplitudes)
        if factor == 0.5:
            np.testing.assert_array_equal(nc, amplitudes)
            rows.append(ns)
        else:
            rows.extend([nc, ns])
        angles.append(theta)
        coefficients.append((cosine, sine_rounded))
    vectors = np.array(rows)
    basis, _ = np.linalg.qr(vectors.T, mode='reduced')
    coordinates = vectors @ basis
    np.testing.assert_allclose(coordinates @ coordinates.T, vectors @ vectors.T,
                               atol=ATOL, rtol=0)
    gammas = [np.kron(X, I), np.kron(Y, I), np.kron(Z, X),
              np.kron(Z, Y), np.kron(Z, Z)]
    sources = [sum(coefficient * gamma for coefficient, gamma in zip(row, gammas))
               for row in coordinates]
    m, nzero, ns_low, nc_high, ns_high = sources
    return delta, angles, coefficients, m, nzero, [(m, ns_low), (nc_high, ns_high)]


def _query(source, cosine_mask, sine_mask):
    """Physical order b,a,target,core; both flags are reused by later calls."""
    core = source.shape[0]
    lifted = []
    for mask in (cosine_mask, sine_mask):
        mn, nm = source @ mask, mask @ source
        scalar = np.block([[(mn + nm) / 2, (mn - nm) / 2],
                           [(mn - nm) / 2, (mn + nm) / 2]])
        lifted.append(np.block([
            [np.kron(I, scalar[a * core:(a + 1) * core,
                               aa * core:(aa + 1) * core]) for aa in range(2)]
            for a in range(2)]))
    zero = np.zeros_like(lifted[0])
    selected = np.block([[lifted[0], zero], [zero, lifted[1]]])
    target = np.kron(I, np.kron(-1j * Y, np.eye(core)))
    select = np.block([[np.eye(4 * core), zero], [zero, target]])
    hb = np.kron(H, np.eye(4 * core))
    return hb @ select @ selected @ hb


def _amplify(query):
    accepted = query.shape[0] // 4
    reflection = np.diag([-1] * accepted + [1] * (3 * accepted))
    return -query @ reflection @ query.conj().T @ reflection @ query, reflection


def _apply_layer(state, n, depth, active, inactive):
    """Apply genuine suffix-controlled blocks, retaining every flag sector."""
    core = active.shape[0] // 8
    result = np.empty_like(state)
    for prefix in range(1 << depth):
        for suffix in range(1 << (n - depth - 1)):
            indices = [(flag * (1 << n)
                        + (prefix << (n - depth))
                        + (target << (n - depth - 1)) + suffix) * core + dirty
                       for flag in range(4) for target in range(2)
                       for dirty in range(core)]
            result[indices] = (active if suffix == 0 else inactive) @ state[indices]
    return result


class HopfErrorAccumulationTests(unittest.TestCase):
    def test_equal_layer_magnitudes_have_base_and_sign_independent_spectrum(self):
        generator = np.random.default_rng(712)
        for n in range(1, 5):
            differences = 0.09 + 0.025 * np.arange(n)
            predicted = np.array([1 + 0j])
            for difference in differences:
                predicted = np.concatenate([
                    np.roots([1, -math.cos(difference) * (value + 1), value])
                    for value in predicted])
            expected_phases = np.sort(np.angle(predicted))
            lower_bound = math.sqrt(2 - 2 * np.prod(np.cos(differences)))
            for pattern in range(3):
                base = generator.uniform(0.4, 1.1, (1 << n) - 1)
                changes = np.concatenate([
                    difference * (np.ones(1 << depth) if pattern == 0 else
                                  generator.choice([-1, 1], 1 << depth))
                    for depth, difference in enumerate(differences)])
                original = direct_real_frame(n, base)
                perturbed = direct_real_frame(n, base + changes)
                relative = original.conj().T @ perturbed
                spectrum = np.linalg.eigvals(relative)
                np.testing.assert_allclose(np.sort(np.angle(spectrum)), expected_phases,
                                           atol=ATOL, rtol=0)
                error = np.linalg.norm(perturbed - original, 2)
                self.assertAlmostEqual(error, np.max(np.abs(predicted - 1)), delta=ATOL)
                self.assertGreaterEqual(error + ATOL, lower_bound)

    def test_ideal_angle_error_is_not_a_maximum_over_layers(self):
        root = math.sqrt(3)
        expected = np.array([[1 / 4, -root / 2, -root / 4, 0],
                             [root / 4, 1 / 2, -3 / 4, 0],
                             [root / 4, 0, 1 / 4, -root / 2],
                             [3 / 4, 0, root / 4, 1 / 2]])
        frame = direct_real_frame(2, [math.pi / 3] * 3)
        np.testing.assert_allclose(frame, expected, atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(hopf_ry(math.pi / 3) - I, 2), 1)
        self.assertAlmostEqual(np.linalg.norm(frame - np.eye(4), 2),
                               math.sqrt(5 + math.sqrt(13)) / 2)
        self.assertAlmostEqual(np.linalg.norm(frame[:, 0] - np.eye(4)[:, 0]),
                               math.sqrt(1.5))

    def test_zero_angle_amplification_flips_the_entire_rejected_space(self):
        source = np.kron(Z, I)
        for zero_mask in (np.kron(X, I), np.kron(Y, X), -np.kron(X, Z)):
            query = _query(source, source, zero_mask)
            amplified, reflection = _amplify(query)
            np.testing.assert_allclose(query.conj().T @ query, np.eye(32),
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(amplified, -reflection, atol=ATOL, rtol=0)
            # Merely returning accepted columns does not justify identity
            # on rejected flags: this sign controls subsequent accumulation.
            self.assertAlmostEqual(np.linalg.norm(amplified - np.eye(32), 2), 2)

    def test_alternating_actual_source_leakage_accumulates_coherently(self):
        dirty = np.array([1, 1j, -1, 1], dtype=complex) / 2
        for n, width in ((2, 20), (4, 24), (6, 28)):
            delta, angles, coefficients, source, nzero, masks = _common_sources(width)
            self.assertLessEqual(delta, (256 * n) ** -2)
            inactive, reflection = _amplify(_query(source, source, nzero))
            np.testing.assert_allclose(inactive, -reflection, atol=ATOL, rtol=0)
            layers = []
            for kind, ((cosine, sine), (nc, ns)) in enumerate(zip(coefficients, masks)):
                self.assertEqual(cosine, 1 - kind * delta)
                defect = 1 - cosine ** 2 - sine ** 2
                self.assertLessEqual(abs(defect - (2 * kind - 1) * delta / 2),
                                     2 * delta ** 1.5)
                query = _query(source, nc, ns)
                layer, _ = _amplify(query)
                embedding = np.eye(32)[:, :8]
                ideal = np.kron(hopf_ry(angles[kind]), np.eye(4))
                self.assertLessEqual(np.linalg.norm(layer @ embedding - embedding @ ideal, 2),
                                     4 * delta)
                # This also checks the exact branch coefficient and its sign.
                np.testing.assert_allclose(layer[16:24, :8],
                                           defect * np.kron(cosine * I + sine * 1j * Y,
                                                            np.eye(4)) / 2,
                                           atol=ATOL, rtol=0)
                self.assertLessEqual(np.linalg.norm(layer + reflection, 2),
                                     6 * math.sqrt(delta))
                layers.append(layer)
            initial = np.zeros(16 * (1 << n), dtype=complex)
            initial[:4] = dirty
            output = initial.copy()
            chosen = []
            for depth in range(n):
                kind = int((n - depth - 1) % 2 == 0)
                chosen.append(kind)
                output = _apply_layer(output, n, depth, layers[kind], inactive)
            witness = np.vdot(dirty, output[8 * (1 << n):8 * (1 << n) + 4])
            self.assertGreaterEqual(abs(witness), n * delta / 8)
            self.assertLessEqual(abs(witness - n * delta / 4),
                                 17 * n * n * delta ** 1.5)
            # Use each physical stage's actual inverse, including leaked work.
            restored = output.copy()
            for depth in reversed(range(n)):
                restored = _apply_layer(restored, n, depth, layers[chosen[depth]].conj().T,
                                        inactive.conj().T)
            np.testing.assert_allclose(restored, initial, atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
