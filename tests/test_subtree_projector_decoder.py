"""Finite algebraic audits of the subtree-projector spectral decoder.

Path products independently reconstruct the literal frame, including signed
and singular tuples. Small unitary completions are verification fixtures only:
this module emits no native loader, Boolean query, or scalable compiler and
does not establish their resource bounds. The matrix checks retain occupied
work and literal relative phases, rather than checking only probabilities.
"""
from __future__ import annotations

from fractions import Fraction as F
import unittest

import numpy as np

from compiler_robust_hopf.frames import real_frame_matrix


ATOL = 3e-11
REGULAR = ((F(3, 5), F(4, 5)), (F(-5, 13), F(12, 13)),
           (F(7, 25), F(-24, 25)), (F(-8, 17), F(-15, 17)))
SINGULAR = ((F(0), F(1)), (F(3, 5), F(-4, 5)),
            (F(1), F(0)), (F(-1), F(0)))


def _geometry(h, pattern=REGULAR):
    """Exact subtree path products, independently of the library recursion."""
    size = 1 << h
    pairs = tuple(pattern[i % len(pattern)] for i in range(size - 1))
    scaling = []
    for depth in range(h + 1):
        vectors = np.full((size, 1 << depth), F(0), dtype=object)
        for leaf in range(size):
            amplitude = F(1)
            for level in range(depth, h):
                node = (1 << level) - 1 + (leaf >> (h - level))
                branch = (leaf >> (h - level - 1)) & 1
                amplitude *= pairs[node][branch]
            vectors[leaf, leaf >> (h - depth)] = amplitude
        scaling.append(vectors)
    frame = np.full((size, size), F(0), dtype=object)
    frame[:, 0] = scaling[0][:, 0]
    tags = [0] * size
    for depth in range(h):
        for prefix in range(1 << depth):
            cosine, sine = pairs[(1 << depth) - 1 + prefix]
            # The marker is the computational word p 1 0...0, with the
            # prescribed -sin(left) + cos(right) orientation.
            marker = (2 * prefix + 1) << (h - depth - 1)
            frame[:, marker] = (-sine * scaling[depth + 1][:, 2 * prefix]
                                + cosine * scaling[depth + 1][:, 2 * prefix + 1])
            tags[marker] = depth + 1
    angles = [np.arctan2(float(sine), float(cosine)) for cosine, sine in pairs]
    return scaling, frame, tags, angles


def _loader_fixtures(scaling, *, leakage=0.0):
    """Same source columns, distant completions, and optional work leakage.

    The one-bit W register is a bounded stand-in for the complete initialized
    work bank. Its occupied subspace is retained in every operator comparison.
    Data-major ordering is (X, W); input anchors are p 0...0, W=0.
    """
    h = len(scaling) - 1
    size = len(scaling[-1])
    labels = 1 << h.bit_length()
    baseline, actual, anchors = [], [], []
    rotation = np.array([[np.cos(leakage), -np.sin(leakage)],
                         [np.sin(leakage), np.cos(leakage)]])
    leak = np.kron(np.eye(size), rotation)
    rng = np.random.default_rng(5071)
    for depth in range(labels):
        effective = min(depth, h)
        indices = [2 * (prefix << (h - effective))
                   for prefix in range(1 << effective)]
        target = np.zeros((2 * size, len(indices)), dtype=complex)
        target[::2] = np.asarray(scaling[effective], dtype=float)
        q, r = np.linalg.qr(np.column_stack((target, np.eye(2 * size))))
        for column in range(len(indices)):
            q[:, column] *= r[column, column] / abs(r[column, column])
        unused = [column for column in range(2 * size) if column not in indices]
        source = np.empty((2 * size, 2 * size), dtype=complex)
        source[:, indices] = q[:, :len(indices)]
        source[:, unused] = q[:, len(indices):]
        change, _ = np.linalg.qr(
            rng.normal(size=(len(unused), len(unused)))
            + 1j * rng.normal(size=(len(unused), len(unused)))
        )
        alternative = source.copy()
        alternative[:, unused] = source[:, unused] @ change
        baseline.append(source)
        actual.append(leak @ alternative)
        anchors.append(indices)
    return baseline, actual, anchors


def _physical_reflections(loaders, anchors, *, active):
    """CZ(a,c) U_D (I-2|ac=11><ac=11| P_anchor,D) U_D^dagger.

    Even disabled blocks multiply the actual loader and its actual inverse;
    they are not replaced with identity by the fixture.
    """
    dimension = len(loaders[0])
    blocks = np.empty((2, len(loaders), dimension, dimension), dtype=complex)
    for control in range(2):
        enabled = active and control
        for depth, loader in enumerate(loaders):
            central = np.eye(dimension)
            if enabled:
                central[anchors[depth], anchors[depth]] = -1
            sign = -1 if enabled else 1
            blocks[control, depth] = sign * loader @ central @ loader.conj().T
    return blocks


def _threshold(k, bit):
    return ((k >> (bit + 1)) << (bit + 1)) + (1 << bit) - 1


def _prepare_threshold(state, bit):
    result = np.empty_like(state)
    for k in range(state.shape[2]):
        for depth in range(state.shape[3]):
            result[:, :, k, depth ^ _threshold(k, bit)] = state[:, :, k, depth]
    return result


def _hadamard(state, bit):
    result = np.empty_like(state)
    mask = 1 << bit
    for k in range(state.shape[2]):
        if not k & mask:
            zero, one = state[:, :, k], state[:, :, k | mask]
            result[:, :, k] = (zero + one) / np.sqrt(2)
            result[:, :, k | mask] = (zero - one) / np.sqrt(2)
    return result


def _search(state, reflections, *, inverse=False):
    """Apply the physical search to axes (X,W,K,D,Z,batch)."""
    bits = (state.shape[2] - 1).bit_length()
    order = range(bits) if inverse else range(bits - 1, -1, -1)
    for bit in order:
        state = _hadamard(_prepare_threshold(state, bit), bit)
        reflected = np.empty_like(state)
        for k in range(state.shape[2]):
            for depth in range(state.shape[3]):
                block = state[:, :, k, depth]
                reflected[:, :, k, depth] = (
                    reflections[(k >> bit) & 1, depth]
                    @ block.reshape(state.shape[0] * state.shape[1], -1)
                ).reshape(block.shape)
        state = _prepare_threshold(_hadamard(reflected, bit), bit)
    return state


def _xor_tag(state, *, active):
    result = np.empty_like(state)
    for k in range(state.shape[2]):
        for tag in range(state.shape[4]):
            result[:, :, k, :, tag ^ (k if active else 0)] = state[:, :, k, :, tag]
    return result


def _initialized_columns(frame, tags, labels, *, extract=False, xor=False):
    size = len(frame)
    result = np.zeros((size, 2, labels, labels, labels, size * labels), dtype=complex)
    for column, sector in enumerate(tags):
        for tag in range(labels):
            result[:, 0, sector if extract else 0, 0, tag ^ (sector if xor else 0),
                   column * labels + tag] = frame[:, column]
    return result


class SubtreeProjectorDecoderTests(unittest.TestCase):
    def test_exact_nested_projectors_and_literal_marker_bands(self):
        for h in (1, 2, 3):
            for pattern in (REGULAR, SINGULAR):
                with self.subTest(h=h, singular=pattern is SINGULAR):
                    scaling, frame, tags, angles = _geometry(h, pattern)
                    np.testing.assert_allclose(
                        np.asarray(frame, dtype=float), real_frame_matrix(angles), atol=ATOL
                    )
                    identity = np.eye(1 << h, dtype=int)
                    np.testing.assert_array_equal(frame.T @ frame, identity)
                    projectors = [vectors @ vectors.T for vectors in scaling]
                    bands = [projectors[0]] + [projectors[d] - projectors[d - 1]
                                               for d in range(1, h + 1)]
                    np.testing.assert_array_equal(projectors[-1], identity)
                    for d, projector in enumerate(projectors):
                        for e, other in enumerate(projectors):
                            np.testing.assert_array_equal(projector @ other, projectors[min(d, e)])
                        for sector, band in enumerate(bands):
                            columns = frame[:, [j for j, tag in enumerate(tags) if tag == sector]]
                            np.testing.assert_array_equal(band, columns @ columns.T)
                            np.testing.assert_array_equal(
                                (2 * projector - identity) @ band,
                                (-1 if sector > d else 1) * band,
                            )
                            self.assertEqual(np.trace(band), 1 if sector == 0 else 1 << (sector - 1))
                    if pattern is SINGULAR and h > 1:
                        # A vanished incoming branch still has a unit marker.
                        np.testing.assert_array_equal(frame[:1 << (h - 1), 0], 0)
                        marker = frame[:, 1 << (h - 2)]
                        self.assertEqual(marker @ marker, F(1))

    def test_integer_search_has_reversible_history_and_padded_thresholds(self):
        def search(sector, word, h, depth_word=0, inverse=False):
            bits = h.bit_length()
            order = range(bits) if inverse else range(bits - 1, -1, -1)
            for bit in order:
                threshold = depth_word ^ _threshold(word, bit)
                word ^= int(sector > threshold) << bit
            return word

        for h in (1, 2, 4, 5, 7, 8, 15, 16, 31, 32):
            labels = 1 << h.bit_length()
            for sector in range(h + 1):
                self.assertEqual(search(sector, 0, h), sector)
                for initial in range(labels):
                    for occupied_depth in (0, labels - 1):
                        result = search(sector, initial, h, occupied_depth)
                        self.assertEqual(
                            search(sector, result, h, occupied_depth, inverse=True), initial
                        )
        self.assertEqual(_threshold(2, 0), 2)  # A valid h=2 branch requests d=h.

    def test_coherent_search_and_old_tag_xor_preserve_literal_phases(self):
        for pattern in (REGULAR, SINGULAR):
            scaling, rational_frame, tags, _ = _geometry(2, pattern)
            frame = np.asarray(rational_frame, dtype=float)
            _, loaders, anchors = _loader_fixtures(scaling)
            reflections = _physical_reflections(loaders, anchors, active=True)
            initial = _initialized_columns(frame, tags, 4)
            extracted = _search(initial, reflections)
            np.testing.assert_allclose(
                extracted, _initialized_columns(frame, tags, 4, extract=True), atol=ATOL
            )
            actual = _search(_xor_tag(extracted, active=True), reflections, inverse=True)
            expected = _initialized_columns(frame, tags, 4, xor=True)
            # All columns include every occupied old tag, so this also checks
            # arbitrary coherent sector/tag superpositions without phase quotients.
            np.testing.assert_allclose(actual, expected, atol=ATOL)
            no_sign = reflections.copy()
            no_sign[1] *= -1
            self.assertGreater(np.linalg.norm(_search(initial, no_sign) - extracted), 1)
            skipped_padding = reflections.copy()
            skipped_padding[1, 2:] = -np.eye(8)
            self.assertGreater(np.linalg.norm(_search(initial, skipped_padding) - extracted), 1)

    def test_inactive_identity_on_every_occupied_work_and_tag_column(self):
        scaling, _, _, _ = _geometry(1)
        baseline, loaders, anchors = _loader_fixtures(scaling, leakage=0.013)
        reflections = _physical_reflections(loaders, anchors, active=False)
        for block in reflections.reshape(-1, 4, 4):
            np.testing.assert_allclose(block, np.eye(4), atol=ATOL)
        # The complete 32-dimensional space includes arbitrary X,W,K,D,Z;
        # tensoring an untouched reference preserves this operator identity.
        initial = np.eye(32, dtype=complex).reshape(2, 2, 2, 2, 2, 32)
        searched = _search(initial, reflections)
        np.testing.assert_allclose(searched, initial, atol=ATOL)
        actual = _search(_xor_tag(searched, active=False), reflections, inverse=True)
        np.testing.assert_allclose(actual, initial, atol=ATOL)
        unguarded_xor = _search(_xor_tag(searched, active=True), reflections, inverse=True)
        self.assertAlmostEqual(np.linalg.norm((unguarded_xor - initial).reshape(32, 32), 2), 2)
        unguarded_sign = reflections.copy()
        unguarded_sign[1] *= -1
        self.assertGreater(np.linalg.norm(_search(initial, unguarded_sign) - initial), 1)
        # Replacing the actual inverse with an ideal completion breaks even a
        # disabled central gate, despite the small initialized-source error.
        self.assertGreater(np.linalg.norm(loaders[0] @ baseline[0].conj().T - np.eye(4), 2), 0.5)

    def test_full_projector_bound_with_leakage_and_arbitrary_completion(self):
        scaling, rational_frame, tags, _ = _geometry(2, SINGULAR)
        frame = np.asarray(rational_frame, dtype=float)
        angle = 0.009
        baseline, loaders, anchors = _loader_fixtures(scaling, leakage=angle)
        _, exact_loaders, _ = _loader_fixtures(scaling)
        deltas = []
        for depth, (loader, indices) in enumerate(zip(loaders, anchors)):
            effective = min(depth, 2)
            target = np.zeros((8, len(indices)), dtype=complex)
            target[::2] = np.asarray(scaling[effective], dtype=float)
            delta = np.linalg.norm(loader[:, indices] - target, 2)
            deltas.append(delta)
            self.assertAlmostEqual(delta, 2 * np.sin(angle / 2), places=12)
            self.assertGreater(np.linalg.norm(loader - baseline[depth], 2), 0.5)
            anchor = np.zeros((8, 8))
            anchor[indices, indices] = 1
            actual = loader @ (np.eye(8) - 2 * anchor) @ loader.conj().T
            vectors = np.asarray(scaling[effective], dtype=float)
            ideal = np.eye(8) - 2 * np.kron(vectors @ vectors.T, np.diag([1, 0]))
            error = np.linalg.norm(actual - ideal, 2)
            self.assertLessEqual(error, 4 * delta + ATOL)
            self.assertGreater(np.linalg.norm((actual - ideal)[:, 1::2], 2), 0.001)
            # With zero source error, the entire conjugated projector is
            # independent of the unconstrained occupied-input completion.
            completed = exact_loaders[depth]
            np.testing.assert_allclose(
                completed @ (np.eye(8) - 2 * anchor) @ completed.conj().T, ideal, atol=ATOL
            )
        # Padded d=3 preserves every data column only on initialized W. Its
        # actual loader is deliberately far from the full identity.
        self.assertGreater(np.linalg.norm(loaders[3] - np.eye(8), 2), 0.5)
        reflections = _physical_reflections(loaders, anchors, active=True)
        initial = _initialized_columns(frame, tags, 4)
        actual = _search(_xor_tag(_search(initial, reflections), active=True),
                         reflections, inverse=True)
        expected = _initialized_columns(frame, tags, 4, xor=True)
        error = np.linalg.norm((actual - expected).reshape(-1, 16), 2)
        self.assertGreater(error, 0.001)
        self.assertLessEqual(error, 8 * 2 * max(deltas) + ATOL)
        self.assertGreater(np.linalg.norm(actual[:, 1].reshape(-1, 16), 2), 0.001)


if __name__ == "__main__":
    unittest.main()
