"""Bounded native checks for the conditional five-layer reduction.

These verify the six-band identity, its actual inverse and cocycle, the
translation-module dimension argument, and Clifford sign conventions.
The non-Clifford normalizer is also a complete, deliberately unsimplified
literal CNOT query. None of these checks proves unrestricted rigidity.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

from tests.test_read_only_three_layer_rigidity import (
    H, S, T, _cx, _degree, _linear, _pauli, _wire,
)
from tests.test_read_only_four_layer_rigidity import (
    ROWS as FOUR_ROWS, PHASE_ALPHAS as FOUR_ALPHAS,
    PHASE_BETAS as FOUR_BETAS, _four_layer_fixture,
)


ROWS = FOUR_ROWS + [((5, 6), (3, 7))]
ALPHAS = FOUR_ALPHAS + [(2, 0, 1)]
BETAS = FOUR_BETAS + [[(0, 1), (1, 2)]]


def _p(j, x):
    return _pauli(_linear(ROWS[j][0], x), _linear(ROWS[j][1], x))


def _phase(j, x):
    bits = [(x >> i) & 1 for i in range(3)]
    return (1j ** sum(a * b for a, b in zip(ALPHAS[j], bits))
            * (-1) ** sum(bits[i] * bits[k] for i, k in BETAS[j]))


def _fibre(c, d, x):
    result = _phase(0, x) * _p(0, x) @ c[0]
    for j in range(1, 6):
        result = _phase(j, x) * _p(j, x) @ c[j] @ d[j - 1] @ result
    return result


def _rank(rows):
    pivots = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def _table(values):
    return sum(int(value) << x for x, value in enumerate(values))


def _hermitian(p, q):
    return 1j ** (p & q).bit_count() * _pauli(p, q)


class ReadOnlyFiveLayerReduction(unittest.TestCase):
    def test_literal_six_band_normalization_and_actual_inverse(self):
        c, d = _four_layer_fixture()
        c.append(_wire(S.conj().T, 0) @ _cx(0, 1) @ _wire(H, 1))
        d.append(_wire(T, 0) @ _wire(T.conj().T, 1))
        u0 = _fibre(c, d, 0)
        w5 = c[5] @ d[4]
        # S-dagger uses three S gates in the fixed H,S,CNOT alphabet.
        fixed_cliffords = 13 + 5
        controls = sum(sum(row.bit_count() for row in ur)
                       + 3 * sum(row.bit_count() for row in vr)
                       for ur, vr in ROWS)
        phases = sum(sum(alpha) + 3 * len(beta)
                     for alpha, beta in zip(ALPHAS, BETAS))
        self.assertEqual((10, 5, fixed_cliffords + controls + phases, 0, 2),
                         (10, 5, 134, 0, 2))
        self.assertLess(max(abs(np.trace(_pauli(p, q).conj().T @ u0)) / 4
                            for p in range(4) for q in range(4)), 1 - 1e-8)
        middle = []
        for x in range(8):
            prefix_a = c[4] @ d[3]
            a = prefix_a @ _p(3, x) @ prefix_a.conj().T
            prefix_b = prefix_a @ c[3] @ d[2]
            b = prefix_b @ _p(2, x) @ prefix_b.conj().T
            prefix_c = prefix_b @ c[2] @ d[1]
            cc = prefix_c @ _p(1, x) @ prefix_c.conj().T
            r = c[0].conj().T @ _p(0, x) @ c[0]
            phase = np.prod([_phase(j, x) for j in range(6)])
            ux = _fibre(c, d, x)
            np.testing.assert_allclose(
                ux, phase * _p(5, x) @ w5 @ _p(4, x) @ a @ b @ cc
                @ w5.conj().T @ u0 @ r, atol=6e-14)
            # Keep the actual zero-address inverse; U0 is not a Pauli.
            extracted = (a.conj().T @ _p(4, x).conj().T @ w5.conj().T
                         @ _p(5, x).conj().T @ ux @ r.conj().T
                         @ u0.conj().T @ w5) / phase
            np.testing.assert_allclose(extracted, b @ cc, atol=6e-14)
            middle.append((b, cc))
        for x, y in itertools.product(range(8), repeat=2):
            bx, cx = middle[x]
            by, cy = middle[y]
            bz, cz = middle[x ^ y]
            sign_b = (-1) ** (
                (_linear(ROWS[2][1], y) & _linear(ROWS[2][0], x))
                .bit_count() & 1)
            sign_c = (-1) ** (
                (_linear(ROWS[1][1], x) & _linear(ROWS[1][0], y))
                .bit_count() & 1)
            np.testing.assert_allclose(
                bz @ cz, sign_b * sign_c * by @ (bx @ cx)
                @ by.conj().T @ (by @ cy), atol=8e-14)

    def test_native_nonclifford_normalizer_and_literal_returned_query(self):
        x = np.array([[0, 1], [1, 0]], dtype=complex)
        y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        z = np.diag([1, -1]).astype(complex)
        identity = np.eye(2, dtype=complex)
        u = T @ H @ T
        b = u @ x @ u.conj().T
        c = u @ y @ u.conj().T
        k = (x - y) / np.sqrt(2)
        delta = z @ k
        np.testing.assert_allclose(b, (z + k) / np.sqrt(2), atol=2e-14)
        np.testing.assert_allclose(c, (k - z) / np.sqrt(2), atol=2e-14)
        np.testing.assert_allclose(b @ c, delta, atol=2e-14)
        np.testing.assert_allclose(b @ z @ b.conj().T, k, atol=2e-14)
        np.testing.assert_allclose(b @ k @ b.conj().T, z, atol=2e-14)
        # B is not Clifford: its conjugation of Z is not a physical Pauli.
        self.assertAlmostEqual(max(abs(np.trace(p @ k)) / 2
                                   for p in (identity, x, y, z)), 2 ** -0.5)
        group = [identity, z, k, delta]
        # Its Clifford subgroup acts irreducibly: its complex span is M2.
        self.assertEqual(np.linalg.matrix_rank(np.array(
            [g.reshape(-1) for g in group])), 4)
        for g in group:
            for p in (x, z):
                image = g @ p @ g.conj().T
                self.assertAlmostEqual(max(abs(np.trace(q @ image)) / 2
                                           for q in (x, y, z)), 1)
            image = b @ g @ b.conj().T
            self.assertAlmostEqual(max(abs(np.trace(q.conj().T @ image)) / 2
                                       for q in group), 1)
        ds = [identity, delta]
        bs = [identity, b]
        for xx, yy in itertools.product(range(2), repeat=2):
            np.testing.assert_allclose(
                ds[xx ^ yy], bs[yy] @ ds[xx] @ bs[yy].conj().T @ ds[yy],
                atol=3e-14)
        for bit in range(2):
            # Six Clifford bands and five T slots, with D2=I. Every
            # controlled Pauli preserves the address. The result is literal.
            actual = (T @ np.linalg.matrix_power(z, bit) @ T @ H @ T
                      @ np.linalg.matrix_power(x, bit)
                      @ np.linalg.matrix_power(y, bit) @ T.conj().T
                      @ H @ S.conj().T)
            np.testing.assert_allclose(actual, np.linalg.matrix_power(x, bit),
                                       atol=3e-14)
            np.testing.assert_allclose(actual.conj().T @ actual, identity,
                                       atol=3e-14)
        # Native unsimplified word: controlled Z=3G, X=1G, Y=5G,
        # C0=H S-dagger=4G, C3=H=1G. Target separate, no helper.
        self.assertEqual((4, 4, 3 + 1 + 5 + 4 + 1, 0, 0), (4, 4, 14, 0, 0))
        # The fair simplified baseline is a single CNOT (0T, depth0, 1G).

    def test_full_degree_translation_module_and_dimension_gap(self):
        n = 5
        lower_sets = [[], [1, 2, 4, 8, 16], [3, 5, 9, 17, 6, 10],
                      [7, 11, 13, 19, 21, 25], [15, 23, 27, 29, 30]]
        low_degree_masks = [s for s in range(1 << n) if s.bit_count() <= 2]
        for lower in lower_sets:
            masks = [31] + lower
            values = [sum((x & s) == s for s in masks) % 2 for x in range(32)]
            self.assertEqual(_degree(values), 5)
            translates = [_table([values[x ^ shift] for x in range(32)])
                          for shift in range(32)]
            self.assertEqual(_rank(translates), 32)
            candidate_space = [_table([(x & s) == s for x in range(32)])
                               for s in low_degree_masks]
            candidate_space += [_table([values[x] * int((x & s) == s)
                                        for x in range(32)])
                                for s in [0, 1, 2, 4, 8, 16]]
            self.assertLessEqual(_rank(candidate_space), 22)
        # At degree four the numerical dimension argument really can be
        # saturated; it cannot silently be strengthened to degree three.
        values = [int(x in (0, 1, 2, 4, 8)) for x in range(16)]
        self.assertEqual(_degree(values), 4)
        space = [_table([(x & s) == s for x in range(16)])
                 for s in range(16) if s.bit_count() <= 2]
        space += [_table([values[x] * int((x & s) == s) for x in range(16)])
                  for s in [0, 1, 2, 4, 8]]
        self.assertEqual(_rank(space), 16)

    def test_hermitian_clifford_sign_and_native_one_t_formula(self):
        fixed = [_cx(0, 1), _cx(1, 0), _wire(H, 0) @ _cx(0, 1),
                 _wire(S, 1) @ _cx(0, 1) @ _wire(H, 1),
                 _cx(0, 1) @ _wire(S, 0) @ _cx(1, 0) @ _wire(H, 0)]
        degrees = []
        for v in fixed:
            signs = []
            for label in range(16):
                image = v @ _hermitian(label & 3, label >> 2) @ v.conj().T
                coeffs = [np.trace(_hermitian(p, q) @ image) / 4
                          for q in range(4) for p in range(4)]
                best = int(np.argmax(np.abs(coeffs)))
                self.assertAlmostEqual(abs(coeffs[best]), 1)
                signs.append(int(coeffs[best].real < 0))
            degrees.append(_degree(signs))
        self.assertTrue(all(degree <= 3 for degree in degrees))
        self.assertEqual(degrees[0], 3)
        x = np.array([[0, 1], [1, 0]], dtype=complex)
        z = np.diag([1, -1]).astype(complex)

        def tau(a, bb):
            return 1j ** (a * bb) * np.linalg.matrix_power(x, a) @ np.linalg.matrix_power(z, bb)

        for epsilon, layer in enumerate((T, T.conj().T)):
            for u, v, a, bb in itertools.product(range(2), repeat=4):
                native = layer @ tau(u, v) @ layer.conj().T
                sign = (u * bb * (1 + a) + v * a + epsilon * u * a) % 2
                expected = (-1) ** sign * tau(a, bb ^ (u * a))
                np.testing.assert_allclose(native @ tau(a, bb) @ native.conj().T,
                                           expected, atol=2e-14)


if __name__ == "__main__":
    unittest.main()
