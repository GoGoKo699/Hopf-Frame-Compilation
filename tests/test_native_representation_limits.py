"""Exact gate identities and finite witnesses for native representation limits.

Integer multiplicities retain the full spin representation. Realification
checks use initialized isometries with clean leakage and arbitrary dirty data.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
import math
import unittest

import numpy as np

from tests.test_coarse_prefix_encoder import _prefix_frame


ATOL = 3e-12


@dataclass(frozen=True)
class _Cyclotomic:
    """Exact a + b sqrt(2) + i(c + d sqrt(2)) amplitudes."""
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)
    c: Fraction = Fraction(0)
    d: Fraction = Fraction(0)

    @staticmethod
    def of(value):
        return value if isinstance(value, _Cyclotomic) else _Cyclotomic(Fraction(value))

    @property
    def coefficients(self):
        return self.a, self.b, self.c, self.d

    def __add__(self, other):
        return _Cyclotomic(*(u + v for u, v in zip(self.coefficients, self.of(other).coefficients)))

    __radd__ = __add__

    def __neg__(self):
        return _Cyclotomic(*(-v for v in self.coefficients))

    def __sub__(self, other):
        return self + (-self.of(other))

    def __rsub__(self, other):
        return self.of(other) + (-self)

    def __mul__(self, other):
        a, b, c, d = self.coefficients
        e, f, g, h = self.of(other).coefficients
        return _Cyclotomic(a * e + 2 * b * f - c * g - 2 * d * h,
                           a * f + b * e - c * h - d * g,
                           a * g + 2 * b * h + c * e + 2 * d * f,
                           a * h + b * g + c * f + d * e)

    __rmul__ = __mul__

    def conjugate(self):
        return _Cyclotomic(self.a, self.b, -self.c, -self.d)


def _adjoint(matrix):
    return np.array([[_Cyclotomic.of(x).conjugate() for x in row] for row in matrix.T], dtype=object)


def _realification(matrix):
    return np.block([[matrix.real, -matrix.imag], [matrix.imag, matrix.real]])


def _rotation(theta):
    return np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])


class NativeRepresentationLimitTests(unittest.TestCase):
    def test_unit_circle_ring_points_by_integer_enumeration(self):
        for exponent in range(7):
            scale, solutions = 1 << exponent, set()
            for a, c in product(range(-scale, scale + 1), repeat=2):
                remainder = scale ** 2 - a * a - c * c
                if remainder < 0 or remainder % 2:
                    continue
                for b in range(-math.isqrt(remainder // 2), math.isqrt(remainder // 2) + 1):
                    d_squared = remainder // 2 - b * b
                    d_abs = math.isqrt(d_squared)
                    if d_abs * d_abs != d_squared:
                        continue
                    for d in {d_abs, -d_abs}:
                        if a * b + c * d == 0:
                            solutions.add((a, b, c, d))
            expected = {(scale, 0, 0, 0), (-scale, 0, 0, 0),
                        (0, 0, scale, 0), (0, 0, -scale, 0)}
            if exponent:
                expected |= {(0, x * (scale // 2), 0, y * (scale // 2))
                             for x, y in product((-1, 1), repeat=2)}
            self.assertEqual(solutions, expected)

    def test_exact_native_realification_gate_words_and_literal_orientation(self):
        q = _Cyclotomic
        halfroot, imaginary = q(0, Fraction(1, 2)), q(0, 0, 1)
        h = np.array([[halfroot, halfroot], [halfroot, -halfroot]], dtype=object)
        s, t = np.diag([q(1), imaginary]), np.diag([q(1), halfroot * (1 + imaginary)])
        z, x = np.diag([q(1), q(-1)]), np.array([[q(), q(1)], [q(1), q()]], dtype=object)
        identity, zero = np.eye(2, dtype=object), np.zeros((2, 2), dtype=object)
        cz = np.diag([q(1), q(1), q(1), q(-1)])
        cnot = np.block([[identity, zero], [zero, x]])

        def exact_equal(left, right):
            self.assertTrue(all(q.of(v) == q() for v in (left - right).flat))

        a = s @ h @ t @ h @ _adjoint(s)
        exact_equal(a @ z @ _adjoint(a), h)
        controlled_h = np.kron(identity, a) @ cz @ np.kron(identity, _adjoint(a))
        exact_equal(controlled_h, np.block([[identity, zero], [zero, h]]))
        real_t = controlled_h @ cz
        exact_equal(real_t, np.block([[identity, zero], [zero, h @ z]]))
        exact_equal(cnot @ cz, np.block([[identity, zero], [zero, x @ z]]))
        exact_equal(real_t.T @ real_t, np.eye(4, dtype=object))
        self.assertTrue(all(q.of(v).c == q.of(v).d == 0 for v in real_t.flat))
        # Linearity reduces the HS diagonalization to the scalars 1 and i.
        basis = h @ s
        for real, imag in ((1, 0), (0, 1)):
            block = np.array([[q(real), q(-imag)], [q(imag), q(real)]], dtype=object)
            exact_equal(basis @ block @ _adjoint(basis),
                        np.diag([q(real) + imaginary * imag, q(real) - imaginary * imag]))

    def test_realification_preserves_rectangular_and_full_isometry_error(self):
        rng = np.random.default_rng(606)
        rectangular = rng.normal(size=(20, 7)) + 1j * rng.normal(size=(20, 7))
        self.assertAlmostEqual(np.linalg.norm(_realification(rectangular), 2), np.linalg.norm(rectangular, 2), places=11)
        # Two logical, one dirty and one initialized clean qubit.
        raw = rng.normal(size=(16, 16)) + 1j * rng.normal(size=(16, 16))
        unitary, _ = np.linalg.qr(raw)
        insertion = np.zeros((16, 8))
        insertion[::2] = np.eye(8)
        frame = _prefix_frame(2, [_rotation(theta) for theta in rng.uniform(-3, 3, 3)]).real
        target = np.kron(frame, np.eye(2))
        error = unitary @ insertion - insertion @ target
        real_insertion, real_target = np.kron(np.eye(2), insertion), np.kron(np.eye(2), target)
        real_error = _realification(unitary) @ real_insertion - real_insertion @ real_target
        np.testing.assert_allclose(real_error, _realification(error), atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(real_error, 2), np.linalg.norm(error, 2), places=11)
        initialized_rebit = np.vstack((np.eye(8), np.zeros((8, 8))))
        self.assertLessEqual(np.linalg.norm(real_error @ initialized_rebit, 2), np.linalg.norm(error, 2) + ATOL)
        accepted = insertion.T @ unitary @ insertion
        gram = 2 * np.eye(8) - target.T @ accepted - accepted.conj().T @ target
        np.testing.assert_allclose(error.conj().T @ error, gram, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(unitary[1::2] @ insertion, 2), .1)

    def test_literal_hopf_marker_gap_from_every_exact_ring_rotation(self):
        for height in range(3, 7):
            coins = [np.eye(2) for _ in range((1 << height) - 1)]
            coins[(1 << (height - 1)) - 1] = _rotation(math.pi / 8)
            marker = _prefix_frame(height, coins)[:, 1]
            distances = []
            for index in range(8):
                candidate = np.zeros(1 << height)
                candidate[:2] = [-math.sin(index * math.pi / 4), math.cos(index * math.pi / 4)]
                distances.append(np.linalg.norm(marker - candidate))
            self.assertAlmostEqual(min(distances), 2 * math.sin(math.pi / 16), places=12)

    def test_even_sign_character_and_doubled_spin_multiplicities(self):
        for height in (3, 4):
            size, pairs = 1 << height, 1 << (height - 1)
            edges = [(p << (height - depth), (2 * p + 1) << (height - depth - 1))
                     for depth in range(height) for p in range(1 << depth)]
            basis = {}
            for first, second in edges:
                word = (1 << first) | (1 << second)
                while word:
                    pivot = word.bit_length() - 1
                    if pivot in basis:
                        word ^= basis[pivot]
                    else:
                        basis[pivot] = word
                        break
            self.assertEqual(len(basis), size - 1)
            full = (1 << size) - 1
            characters = Counter(min(word, full ^ word) for word in range(1 << size))
            self.assertEqual(len(characters), 1 << (size - 1))
            self.assertEqual(set(characters.values()), {2})
            modulus, counts = 3 ** pairs, Counter({0: 1})
            for index in range(pairs):
                next_counts = Counter()
                for value, multiplicity in counts.items():
                    for digit in (-1, 0, 0, 1):
                        next_counts[(value + digit * 3 ** index) % modulus] += multiplicity
                counts = next_counts
            self.assertEqual(sum(counts.values()), 1 << size)
            self.assertEqual(len(counts), modulus)
            self.assertEqual(max(counts.values()), 1 << pairs)
            residues = set()
            for digits in product((-1, 0, 1), repeat=pairs):
                value = sum(digit * 3 ** index for index, digit in enumerate(digits)) % modulus
                self.assertNotIn(value, residues)
                residues.add(value)
                self.assertEqual(counts[value], 1 << digits.count(0))
            for index in range(pairs):
                self.assertEqual(counts[3 ** index], 1 << (pairs - 1))
                self.assertEqual(counts[-3 ** index % modulus], 1 << (pairs - 1))

    def test_spin_target_dependent_gap_starts_at_height_four(self):
        for height in (3, 4):
            size, pairs = 1 << height, 1 << (height - 1)
            target_ratio = 2 ** (height + 1 - pairs)
            max_ratio = 2 ** (height + 2 - pairs)
            if height == 3:
                self.assertEqual(target_ratio, 1)
                self.assertEqual(max_ratio, 2)
            else:
                self.assertLess(target_ratio, 1)
                self.assertLess(max_ratio, 1)
                epsilon = 2 ** -size
                self.assertGreater(2 * math.sin(math.pi / (3 ** pairs)), epsilon)
                self.assertGreater(math.sin(1 / (8 * 3 ** pairs)), epsilon)


if __name__ == '__main__':
    unittest.main()
