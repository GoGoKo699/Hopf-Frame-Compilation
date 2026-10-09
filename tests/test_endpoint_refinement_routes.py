"""Scoped dense-mask, source-program and all-order tree identities.

These are bounded regressions for the retained route lemmas, not a search
for a compiler or a lower bound on unrestricted circuits. Exact Gaussian
rationals certify the eight-mode witness; native-source checks reuse the
literal gate-word fixtures and retain all physical columns.
"""
from dataclasses import dataclass
from fractions import Fraction as F
import itertools
import math
import unittest

import numpy as np

from tests.test_coarse_prefix_encoder import _prefix_frame
from tests.test_one_clean_compiler import (
    _amplification_word, _encoded_signs, _literal_mask, _paired_source_data,
    _sandwich_word,
)
from tests.test_operator_source_compiler import H, I2, X, Y, Z, _pauli, _word_matrix


ATOL = 5e-11


@dataclass(frozen=True)
class _Gaussian:
    real: F = F(0)
    imag: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'real', F(self.real))
        object.__setattr__(self, 'imag', F(self.imag))

    def __add__(self, other):
        other = _gaussian(other)
        return _Gaussian(self.real + other.real, self.imag + other.imag)

    __radd__ = __add__

    def __neg__(self):
        return _Gaussian(-self.real, -self.imag)

    def __sub__(self, other):
        return self + (-_gaussian(other))

    def __rsub__(self, other):
        return _gaussian(other) + (-self)

    def __mul__(self, other):
        other = _gaussian(other)
        return _Gaussian(self.real * other.real - self.imag * other.imag,
                         self.real * other.imag + self.imag * other.real)

    __rmul__ = __mul__

    def conjugate(self):
        return _Gaussian(self.real, -self.imag)

    def norm_squared(self):
        return self.real ** 2 + self.imag ** 2

    def __truediv__(self, other):
        other = _gaussian(other)
        product = self * other.conjugate()
        return _Gaussian(product.real / other.norm_squared(), product.imag / other.norm_squared())


def _gaussian(value):
    return value if isinstance(value, _Gaussian) else _Gaussian(value)


def _eye(n):
    return np.array([[_Gaussian(i == j) for j in range(n)] for i in range(n)], dtype=object)


def _dag(matrix):
    return np.array([[v.conjugate() for v in row] for row in matrix.T], dtype=object)


def _block(a, b):
    return np.block([[a, np.full((len(a), len(b)), _Gaussian())],
                     [np.full((len(b), len(a)), _Gaussian()), b]])


def _rotation(n, a, b, c, s):
    out = _eye(n)
    out[a, a] = out[b, b] = _Gaussian(c)
    out[a, b], out[b, a] = _Gaussian(-s), _Gaussian(s)
    return out


def _rank(matrix):
    matrix = matrix.copy()
    pivot = 0
    for column in range(matrix.shape[1]):
        row = next((r for r in range(pivot, len(matrix)) if matrix[r, column] != _Gaussian()), None)
        if row is None:
            continue
        matrix[[pivot, row]] = matrix[[row, pivot]]
        matrix[pivot] = [v / matrix[pivot, column] for v in matrix[pivot]]
        for row in range(len(matrix)):
            if row != pivot:
                scalar = matrix[row, column]
                matrix[row] = [v - scalar * w for v, w in zip(matrix[row], matrix[pivot])]
        pivot += 1
        if pivot == len(matrix):
            break
    return pivot


def _rational_frame(n, coins):
    size = 1 << n
    frame = np.array([[F(i == j) for j in range(size)] for i in range(size)], dtype=object)
    for depth in range(n):
        for p in range(1 << depth):
            c, s = coins[depth, p]
            a, b = p << (n - depth), (2 * p + 1) << (n - depth - 1)
            low, high = frame[a].copy(), frame[b].copy()
            frame[a] = c * low - s * high
            frame[b] = s * low + c * high
    return frame


def _coins(n):
    choices = [(F(3, 5), F(4, 5)), (F(5, 13), F(-12, 13)),
               (F(-8, 17), F(15, 17)), (F(7, 25), F(24, 25))]
    return {(d, p): choices[((1 << d) + p) % len(choices)]
            for d in range(n) for p in range(1 << d)}


class EndpointRefinementRouteTests(unittest.TestCase):
    def test_eight_mode_exact_gauge_spectrum_rank_four_and_gap(self):
        c = np.array([[_Gaussian(F((-1) ** ((x & y).bit_count()), 2))
                       for y in range(4)] for x in range(4)], dtype=object)
        eigenvalues = [_Gaussian(4, -3) / 5, _Gaussian(4, 3) / 5,
                       _Gaussian(5, 12) / 13, _Gaussian(-5, 12) / 13,
                       _Gaussian(15, -8) / 17, _Gaussian(15, 8) / 17,
                       _Gaussian(7, 24) / 25, _Gaussian(-7, 24) / 25]
        self.assertEqual(len(set(eigenvalues)), 8)
        self.assertTrue(all(v.norm_squared() == 1 for v in eigenvalues))
        diagonal = np.full((8, 8), _Gaussian(), dtype=object)
        for j, value in enumerate(eigenvalues):
            diagonal[j, j] = value
        child_basis = _block(c, c)
        gauge = child_basis @ diagonal @ _dag(child_basis)
        root = _rotation(8, 0, 4, F(3, 5), F(4, 5))
        transported = _dag(root) @ gauge @ root
        basis = _dag(root) @ child_basis
        np.testing.assert_array_equal(transported @ basis, basis @ diagonal)
        projector = basis[:, :1] @ _dag(basis[:, :1])
        self.assertEqual(projector[0, 0], _Gaussian(F(9, 100)))
        self.assertEqual(projector[0, 1], _Gaussian(F(3, 20)))
        self.assertEqual([v.norm_squared() for v in basis[:, 0]],
                         [F(9, 100), F(1, 4), F(1, 4), F(1, 4), F(4, 25), 0, 0, 0])
        left = (_rotation(4, 0, 1, F(4, 5), F(3, 5))
                @ _rotation(4, 2, 3, F(12, 13), F(5, 13))
                @ _rotation(4, 0, 2, F(8, 17), F(15, 17)))
        right = (_rotation(4, 0, 1, F(15, 17), F(8, 17))
                 @ _rotation(4, 2, 3, F(24, 25), F(7, 25))
                 @ _rotation(4, 0, 2, F(20, 29), F(21, 29)))
        frame = _block(left, right) @ root
        np.testing.assert_array_equal(_block(left, right) @ gauge @ root, frame @ transported)
        correction = transported @ _dag(gauge)
        residual = correction - _eye(8)
        support = np.column_stack([_eye(8)[:, 0], _eye(8)[:, 4], gauge[:, 0], gauge[:, 4]])
        self.assertEqual(_rank(residual), 4)
        self.assertEqual(_rank(np.column_stack([support, residual])), _rank(support))
        self.assertEqual(_rank(np.column_stack([support, _dag(residual)])), _rank(support))
        np.testing.assert_array_equal(correction @ _dag(correction), _eye(8))
        gap_squared = min((a - b).norm_squared() for j, a in enumerate(eigenvalues) for b in eigenvalues[j + 1:])
        self.assertEqual(gap_squared, F(52, 4225))
        magnitudes = [.5, .5, .5, .4, .3, 0, 0, 0]
        overlap = max(sum(magnitudes[:k]) / math.sqrt(k) for k in (1, 2, 4, 8))
        self.assertAlmostEqual(overlap, 19 / 20)
        distance = math.sqrt(39) / 20
        self.assertGreater(math.sqrt(float(gap_squared)) * distance / (1 + distance), .026)

    def test_projector_dyadic_divisibility_includes_degenerate_eigenspaces(self):
        def dyadic(rational):
            denominator = rational.denominator
            return denominator & (denominator - 1) == 0

        # Each Gaussian component of a two-qubit spectral projector is
        # a/4 with |a|<=4, even when multiple eigenvalues coincide.
        for numerator in range(-4, 5):
            self.assertEqual(dyadic(F(3 * numerator, 20)), numerator == 0)
        for left, right in itertools.product((0, 1), repeat=2):
            self.assertEqual(dyadic(F(9 * left + 16 * right, 25)), left == right)

    def test_three_point_hopf_row_has_uniform_flat_mixer_gap(self):
        for n in (3, 4):
            size = 1 << n
            coins = [np.eye(2) for _ in range(size - 1)]
            for node, angle in (((1 << (n - 2)) - 1, np.pi / 4),
                                ((1 << (n - 1)) - 1, math.asin(1 / math.sqrt(3)))):
                coins[node] = np.array([[math.cos(angle), -math.sin(angle)],
                                        [math.sin(angle), math.cos(angle)]])
            frame = _prefix_frame(n, coins)
            expected = np.zeros(size)
            expected[:3] = [1, -1, -1]
            np.testing.assert_allclose(frame[0], expected / math.sqrt(3), atol=ATOL, rtol=0)
            walsh = H
            for _ in range(n - 1):
                walsh = np.kron(H, walsh)
            # A finite phase grid checks the analytic cosine identity;
            # the dimension-independent proof supplies the uniform bound.
            for phases in itertools.product((0, np.pi / 3, np.pi / 2, np.pi), repeat=2):
                phase = np.array([0, *phases])
                cosines = np.cos([phase[0] - phase[1], phase[0] - phase[2], phase[1] - phase[2]])
                self.assertAlmostEqual(sum(cosines ** 2), (3 + abs(sum(np.exp(2j * phase))) ** 2) / 4)
                self.assertGreaterEqual(max(abs(cosines)), .5 - ATOL)
                row = np.zeros(size, complex)
                row[:3] = np.exp(1j * phase) / math.sqrt(3)
                mixed = row @ walsh
                variation = sum(abs(abs(mixed) ** 2 - 1 / size))
                flat_distance = np.linalg.norm(abs(mixed) - 1 / math.sqrt(size))
                self.assertGreaterEqual(variation, 1 / 3 - ATOL)
                self.assertGreaterEqual(flat_distance, 1 / 6 - ATOL)
                self.assertLessEqual(variation, 2 * flat_distance + ATOL)

    def test_signed_prefix_phase_gauge_on_every_column(self):
        for values in ((-.5, .25, 1.25), (-.5, .25, 1.25, 0, -1.2, .7, -.1),
                       (0,) * 7, (1.7, -1.7, 0, 1.6, -1.6, .01, -.01)):
            size = len(values) + 1
            n = size.bit_length() - 1
            row_phase = np.ones(2 * size, complex)
            marker_phase = np.ones(size, complex)
            complex_coins, real_coins = [], []
            for v, t in enumerate(values, start=1):
                z = (1 - 1j * t) / math.sqrt(1 + t * t)
                g, h = 1 / (1 + 1j * t), t / (1 + 1j * t)
                complex_coins.append(np.array([[1j * h, -g.conjugate()], [g, -1j * h.conjugate()]]))
                real_coins.append(np.array([[t, -1], [1, t]]) / math.sqrt(1 + t * t))
                row_phase[2 * v], row_phase[2 * v + 1] = 1j * z * row_phase[v], z * row_phase[v]
                marker_phase[v] = -1j * z.conjugate() ** 2 * row_phase[v].conjugate()
            for depth in range(n + 1):
                count = 1 << depth
                inputs = np.ones(count, complex)
                for j in range(depth):
                    for p in range(1 << j):
                        inputs[(2 * p + 1) << (depth - j - 1)] = marker_phase[(1 << j) + p]
                expected = np.diag(row_phase[count:2 * count]) @ _prefix_frame(depth, real_coins) @ np.diag(inputs)
                np.testing.assert_allclose(_prefix_frame(depth, complex_coins), expected, atol=ATOL, rtol=0)

    def test_signed_support_and_static_orthogonal_history_capacity(self):
        for n in range(1, 4):
            size = 1 << n
            frame = _rational_frame(n, _coins(n))
            self.assertEqual(np.count_nonzero(frame), size * (n + 1))
            self.assertTrue(all(np.count_nonzero(row) == n + 1 for row in frame))
            labels = (n + 1 - 1).bit_length()
            self.assertGreaterEqual((1 << labels) * size, size * (n + 1))
            self.assertLess((1 << (labels - 1)) * size, size * (n + 1))

    def test_exact_resolvent_sums_all_ancestor_insertion_orders(self):
        for n in (2, 3):
            size = (2 << n) - 1
            identity = np.array([[F(i == j) for j in range(size)] for i in range(size)], dtype=object)
            base, target = np.full((size, size), F(0)), np.full((size, size), F(0))
            for v in range(1, 1 << n):
                for b in (0, 1):
                    base[2 * v + b - 1, v - 1] = F((-1) ** (v + b), v + b + 2)
                    target[2 * v + b - 1, v - 1] = F(0) if v % 3 == 0 else F((-1) ** b, v + b + 3)

            def resolvent(shift):
                out, power = identity.copy(), identity.copy()
                for _ in range(n):
                    power = power @ shift
                    out += power
                np.testing.assert_array_equal(power @ shift, np.full((size, size), F(0)))
                return out

            ga, gb = resolvent(target), resolvent(base)
            difference = target - base
            np.testing.assert_array_equal(ga - gb, ga @ difference @ gb)
            np.testing.assert_array_equal(ga, resolvent(gb @ difference) @ gb)

    def test_damped_tree_chain_normalization_bound_in_exact_arithmetic(self):
        for n in range(1, 9):
            for rho in (F(1), F(3, 4), F(1, 2), F(n, n + 1)):
                sums = [sum(rho ** k for k in range(d + 1)) for d in range(n + 1)]
                test_norm_squared = sum(value * value for value in sums) / (n + 1)
                self.assertGreaterEqual(test_norm_squared / rho ** (2 * n), F((n + 1) ** 2, 3))


class SourceProgramTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.width, cls.core, cls.flag = 7, 3, 6
        cls.identity = np.eye(1 << cls.width, dtype=complex)
        source, gammas, weights, fixed = _paired_source_data(2)
        cls.source = np.kron(np.eye(16), source)
        cls.xflag, cls.zflag = _pauli(7, {6: X}), _pauli(7, {6: Z})
        cls.zero, cls.one = (cls.identity + cls.zflag) / 2, (cls.identity - cls.zflag) / 2
        cls.klein = [_pauli(7, {5: X}), _pauli(7, {5: Z})]
        fixed_mask = np.kron(np.eye(16), _word_matrix(3, _literal_mask(fixed)))
        df = cls.source @ fixed_mask @ cls.source @ fixed_mask.conj().T - .5 * cls.identity
        cls.fixed_scalar = .5 * cls.identity + cls.xflag @ df
        _, _, basis = np.linalg.svd(np.sqrt(weights)[None, :], full_matrices=True)
        perpendicular = [np.kron(np.eye(16), sum(a * g for a, g in zip(row, gammas))) for row in basis[1:]]
        cls.generators = perpendicular + [cls.xflag @ _pauli(7, {3 + j: p}) @ cls.klein[j] @ cls.source
                                         for j in range(2) for p in (X, Y, Z)]

    def _primitive(self, axis, rows, dressed=True):
        scalar = np.zeros_like(self.identity)
        for row, t in enumerate(rows):
            mask = np.kron(np.eye(16), _word_matrix(3, _literal_mask(_encoded_signs(2, np.arctan(t)))))
            ng = mask @ self.source @ mask.conj().T
            mean = np.trace(self.source @ ng).real / len(self.identity)
            projector = self.identity if len(rows) == 1 else (self.identity + (1 if row == 0 else -1) * _pauli(7, {3: Z})) / 2
            scalar += projector @ (mean * self.identity + self.xflag @ (self.source @ ng - mean * self.identity))
        first, second = {'X': (Y, Z), 'Y': (Z, X), 'Z': (X, Y)}[axis]
        klein = self.klein[1] if dressed else self.identity
        ra = self.zero + self.one @ _pauli(7, {4: first}) @ klein
        rb = self.zero + self.one @ _pauli(7, {4: second}) @ klein
        a, b = ra @ self.fixed_scalar @ ra, rb @ scalar @ rb
        query = a.conj().T @ b @ a
        return query @ self.zflag @ query.conj().T @ self.zflag @ query @ self.zflag @ query.conj().T @ self.zflag @ query

    def _defect(self, unitary):
        largest = 0
        dimension = len(unitary)
        for generator in self.generators:
            image = unitary @ generator @ unitary.conj().T
            projection = sum((np.trace(other @ image) / dimension) * other for other in self.generators)
            largest = max(largest, np.linalg.norm(image - projection) / math.sqrt(dimension))
        return largest

    def test_klein_dressing_preserves_the_literal_native_word_and_dirty_columns(self):
        native = _word_matrix(7, _amplification_word(_sandwich_word(
            2, _literal_mask(_encoded_signs(2, np.arctan(.25))), 4, 6), 6))
        np.testing.assert_allclose(native, self._primitive('Y', [.25], dressed=False), atol=ATOL, rtol=0)
        ordinary = self._primitive('Y', [.25, 1.25], dressed=False)
        dressed = self._primitive('Y', [.25, 1.25])
        change = self.zero + self.one @ self.klein[1]
        np.testing.assert_allclose(dressed, change @ ordinary @ change, atol=ATOL, rtol=0)
        np.testing.assert_allclose(self.zero @ dressed @ self.zero, self.zero @ ordinary @ self.zero, atol=ATOL, rtol=0)
        np.testing.assert_allclose(dressed.conj().T @ dressed, self.identity, atol=ATOL, rtol=0)

    def test_unequal_program_escapes_the_same_explicit_spin_algebra(self):
        for i, first in enumerate(self.generators):
            for j, second in enumerate(self.generators):
                np.testing.assert_allclose(first @ second + second @ first,
                                           2 * self.identity if i == j else 0, atol=ATOL, rtol=0)
        fixed = -1j * _pauli(7, {4: Y})
        uniform = fixed @ self._primitive('X', [.25, .25]) @ self._primitive('Z', [.25, .25])
        unequal = fixed @ self._primitive('X', [.25, 1.25]) @ self._primitive('Z', [.25, 1.25])
        self.assertLess(self._defect(uniform), ATOL)
        self.assertGreater(self._defect(unequal), .9)
        u, v = F(1, 4), F(5, 4)
        leakage_squared = (u - v) ** 2 / ((1 + u * u) * (1 + v * v))
        self.assertEqual(leakage_squared, F(256, 697))


if __name__ == '__main__':
    unittest.main()
