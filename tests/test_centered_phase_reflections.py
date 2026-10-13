"""Global ansatz identities and a literal finite member, not an endpoint proof."""
from fractions import Fraction
from functools import lru_cache
import math
import unittest

import numpy as np

from compiler_robust_hopf.centered_phase_reflections import (
    boolean_direction_words, centered_word, flat_walsh_direction,
    nearest_center_reflection, optimal_phase_table, ordered_hopf_floor,
    ordered_hopf_witness, plane_rotation, reflection, two_coordinate_floor,
    walsh_matrix,
)
from compiler_robust_hopf.native_coarse_fixture import gate_counts, inverse_word, simulate_columns
from verification.fault_tolerant.exact_arithmetic import CQ2
from tests.test_conjugate_recovery import exact_native


@lru_cache(maxsize=2)
def exact_fixture(inverse=False):
    _, word = boolean_direction_words()
    if inverse:
        word = inverse_word(word)
    return exact_native([(g.lower(), *qs) for g, qs in word], 5)


class CenterPhaseTests(unittest.TestCase):
    def test_reflection_identity_on_all_columns(self):
        rng = np.random.default_rng(10411)
        for size in (4, 8, 16, 32):
            d = np.exp(1j*rng.uniform(-np.pi, np.pi, size))
            w = walsh_matrix(size) @ d / np.sqrt(size); w[0] *= 1j
            r = np.eye(size); r[0, 0] = -1
            np.testing.assert_allclose(centered_word(d), reflection(w) @ r, atol=2e-14, rtol=0)
            np.testing.assert_allclose(centered_word(np.exp(.37j)*d), centered_word(d), atol=2e-14, rtol=0)

    def test_global_phase_optimum(self):
        rng = np.random.default_rng(10412)
        for size in (8, 16, 32):
            for real in (False, True):
                v = rng.normal(size=size).astype(complex)
                if not real:
                    v += 1j*rng.normal(size=size)
                v /= np.linalg.norm(v)
                r = np.eye(size); r[0, 0] = -1
                target = reflection(v) @ r
                d, minimum = optimal_phase_table(v)
                error = np.linalg.norm(centered_word(d) - target, 2)
                self.assertAlmostEqual(error, minimum, places=12)
                for _ in range(3):
                    other = np.exp(1j*rng.uniform(-np.pi, np.pi, size))
                    self.assertGreaterEqual(np.linalg.norm(centered_word(other) - target, 2) + 1e-12, minimum)

    def test_exact_real_family(self):
        rng = np.random.default_rng(10413)
        for size in (4, 8, 16):
            signs = np.array([1]*(size//2) + [-1]*(size//2))
            rng.shuffle(signs)
            b = flat_walsh_direction(signs)
            self.assertAlmostEqual(b[0], 0)
            for theta in (0, -.4, .7, np.pi):
                got = centered_word(np.exp(.5j*theta*signs))
                np.testing.assert_allclose(got, plane_rotation(b, theta), atol=2e-14, rtol=0)
                self.assertLess(np.max(abs(got.imag)), 2e-14)

    def test_nonflat_quadratic_floor(self):
        for size in (4, 8, 16):
            b = np.zeros(size); b[1] = b[2] = 1/np.sqrt(2)
            e = np.eye(size)[0]
            for theta in (.07, -.3, .8, np.pi/2):
                normal = np.cos(theta/2)*e + np.sin(theta/2)*b
                d, _ = optimal_phase_table(normal)
                actual = np.linalg.norm(centered_word(d) - plane_rotation(b, theta), 2)
                floor = two_coordinate_floor(theta)
                self.assertAlmostEqual(actual, floor, places=12)
                t = np.sin(theta/2)**2
                self.assertGreaterEqual(floor + 1e-15, t)
                self.assertLessEqual(floor, np.sqrt(2)*t + 1e-15)

    def test_optimum_stays_stable_near_zero(self):
        size = 8
        b = np.zeros(size); b[1] = b[2] = 1/np.sqrt(2)
        e = np.eye(size)[0]
        for theta in (1e-2, 1e-3, 1e-4):
            _, floor = optimal_phase_table(np.cos(theta/2)*e + np.sin(theta/2)*b)
            self.assertAlmostEqual(floor, two_coordinate_floor(theta), delta=3e-15)

    def test_ordered_hopf_exact_floor_and_attainer(self):
        rng = np.random.default_rng(10414)
        for size in (4, 8, 16):
            r = np.eye(size); r[0, 0] = -1
            for theta in (0, -.2, .6, np.pi/2):
                target = ordered_hopf_witness(theta, size)
                best = nearest_center_reflection(theta, size)
                self.assertAlmostEqual(np.linalg.norm(best - target, 2), ordered_hopf_floor(theta), places=12)
                for _ in range(2):
                    v = rng.normal(size=size) + 1j*rng.normal(size=size)
                    v /= np.linalg.norm(v)
                    self.assertGreaterEqual(np.linalg.norm(reflection(v)@r - target, 2) + 1e-13,
                                            ordered_hopf_floor(theta))

    def test_regular_optimal_table_denominators(self):
        rng = np.random.default_rng(10415)
        for size in (8, 16, 32):
            b = rng.normal(size=size); b[0] = 0; b /= np.linalg.norm(b)
            lam = np.sqrt(size)*(walsh_matrix(size)@b)
            e = np.eye(size)[0]
            for theta in (-np.pi/2, .01, np.pi/2):
                c, s = np.cos(theta/2), np.sin(theta/2)
                denominator = np.sqrt(c*c + s*s*lam*lam)
                self.assertGreaterEqual(denominator.min(), 1/np.sqrt(2) - 1e-14)
                d = (c+1j*s*lam)/denominator
                target = reflection(c*e+s*b)@reflection(e)
                _, optimum = optimal_phase_table(c*e+s*b)
                self.assertAlmostEqual(np.linalg.norm(centered_word(d)-target, 2), optimum, places=11)

    def test_star_tangent_does_not_imply_exact_real_coverage(self):
        rng = np.random.default_rng(10416)
        for size in (4, 8, 16):
            a = rng.normal(size=size); a[0] = 0
            circulant = np.array([[a[x ^ y] for y in range(size)] for x in range(size)])
            e = np.eye(size)[0]; r = reflection(e)
            s = np.eye(size, dtype=complex); s[0, 0] = 1j
            derivative = .5j*s@(circulant@r-r@circulant)@r@s.conj().T
            np.testing.assert_allclose(derivative, np.outer(a, e)-np.outer(e, a), atol=1e-13, rtol=0)

    def test_native_fixture_literal_counts(self):
        _, word = boolean_direction_words()
        self.assertEqual(len(word), 438)
        self.assertEqual(gate_counts(word)['T_count'], 144)
        self.assertLessEqual({g for g, _ in word}, {'H', 'S', 'SDG', 'T', 'TDG', 'CX'})
        self.assertEqual(max(q for _, qs in word for q in qs), 4)

    def test_native_fixture_all_initialized_columns_exact(self):
        full = exact_fixture()
        b = [Fraction(0), Fraction(1, 2), Fraction(0), Fraction(1, 2),
             Fraction(0), Fraction(1, 2), Fraction(0), Fraction(-1, 2)]
        for row in range(32):
            for col in range(8):
                expected = Fraction(0)
                if row < 8:
                    expected = (Fraction(row == col) - Fraction(row == 0 and col == 0)
                                - b[row]*b[col] + (b[row] if col == 0 else 0)
                                - (b[col] if row == 0 else 0))
                self.assertEqual(full[row][col], CQ2(expected))

    def test_native_fixture_full_inverse_exact(self):
        full, back = exact_fixture(), exact_fixture(True)
        for i in range(32):
            for j in range(32):
                self.assertEqual(back[i][j], full[j][i].conj())
        # A separate floating full-space multiplication checks the inverse order.
        _, word = boolean_direction_words()
        u = simulate_columns(word, nqubits=5)
        np.testing.assert_allclose(simulate_columns(inverse_word(word), u, nqubits=5),
                                   np.eye(32), atol=2e-13, rtol=0)

    def test_endpoint_floor_exponents(self):
        for n in range(8, 129):
            N = 1 << n
            root = math.isqrt(N)
            ell = root + int(root*root < N) + 6
            self.assertLess(2*ell + 12, N)

    def test_reject_invalid_diagnostics(self):
        for size in (True, 1, 3, 512):
            with self.assertRaises(ValueError):
                walsh_matrix(size)
        with self.assertRaises(ValueError):
            centered_word([1, 0])
        with self.assertRaises(ValueError):
            reflection([0, 0])
        with self.assertRaises(ValueError):
            flat_walsh_direction([1, 1, 1, -1])
        with self.assertRaises(ValueError):
            plane_rotation([1, 0], .3)
        for angle in (True, float('nan'), 1j):
            with self.assertRaises(ValueError):
                two_coordinate_floor(angle)
        with self.assertRaises(ValueError):
            ordered_hopf_floor(np.pi)


if __name__ == '__main__':
    unittest.main()
