"""Exact finite certificates for the global Haar research chapter.

These checks verify the literal witness identities and an anchored
nonalternating counterexample. The dimension-uniform norm theorem is
proved in STRUCTURAL_COMPILATION_LIMITS.md, not inferred from these tests.
"""
from dataclasses import dataclass
from fractions import Fraction
import unittest


@dataclass(frozen=True)
class Quadratic:
    """An exact element a + b sqrt(2)."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __add__(self, other):
        if not isinstance(other, Quadratic):
            other = Quadratic(Fraction(other))
        return Quadratic(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Quadratic(-self.a, -self.b)

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        if not isinstance(other, Quadratic):
            other = Quadratic(Fraction(other))
        return Quadratic(self.a * other.a + 2 * self.b * other.b,
                         self.a * other.b + self.b * other.a)

    __rmul__ = __mul__


ZERO = Quadratic()
ONE = Quadratic(Fraction(1))
HALF = Quadratic(Fraction(1, 2))
C = Quadratic(Fraction(0), Fraction(1, 2))


def identity(size):
    return [[ONE if i == j else ZERO for j in range(size)]
            for i in range(size)]


def block_diagonal(left, right):
    m, n = len(left), len(right)
    return ([row + [ZERO] * n for row in left]
            + [[ZERO] * m + row for row in right])


def transpose(matrix):
    return list(map(list, zip(*matrix)))


def multiply(left, right):
    columns = transpose(right)
    return [[sum((x * y for x, y in zip(row, col)), ZERO)
             for col in columns] for row in left]


def apply(matrix, vector):
    return [sum((x * y for x, y in zip(row, vector)), ZERO)
            for row in matrix]


def root_rotation(n):
    matrix = identity(1 << n)
    m = 1 << (n - 1)
    matrix[0][0] = matrix[m][m] = C
    matrix[0][m], matrix[m][0] = -C, C
    return matrix


def haar(n):
    if n == 0:
        return identity(1)
    child = haar(n - 1)
    return multiply(block_diagonal(child, child), root_rotation(n))


def inverse_sqrt_power_of_two(exponent):
    if exponent % 2:
        return Quadratic(Fraction(0), Fraction(1, 1 << ((exponent + 1) // 2)))
    return Quadratic(Fraction(1, 1 << (exponent // 2)))


class GlobalHaarCertificates(unittest.TestCase):
    def test_common_witness_conjugation_identity(self):
        for n in (3, 4):
            with self.subTest(n=n):
                size, m = 1 << n, 1 << (n - 1)
                q = haar(n)
                w = block_diagonal(identity(m), haar(n - 1))
                e0 = [ONE] + [ZERO] * (size - 1)
                em = [ZERO] * m + [ONE] + [ZERO] * (m - 1)
                a = inverse_sqrt_power_of_two(n - 1)
                ul, ur = [a] * m + [ZERO] * m, [ZERO] * m + [a] * m
                self.assertEqual(multiply(transpose(q), q), identity(size))
                self.assertEqual(apply(w, e0), e0)
                self.assertEqual(apply(transpose(w), e0), e0)
                self.assertEqual(apply(transpose(w), ur), em)
                self.assertEqual(apply(transpose(q), ur),
                                 [C * (x + y) for x, y in zip(e0, em)])
                conjugated = multiply(multiply(q, w), transpose(q))
                actual = apply(transpose(conjugated), ur)
                left_weight = (ONE - a) * HALF
                right_weight = HALF + a * (HALF - C)
                expected = [left_weight * x + C * y + right_weight * z
                            for x, y, z in zip(ul, em, ur)]
                self.assertEqual(actual, expected)

    def test_anchored_nonalt_candidate_has_asymmetric_moduli(self):
        q = haar(3)
        v = multiply(q, q)
        z = multiply(multiply(transpose(v), root_rotation(3)), v)
        self.assertEqual(z[0][1], Quadratic(Fraction(0), Fraction(-1, 16)))
        self.assertEqual(z[1][0], Quadratic(Fraction(-1, 8)))
        gap = z[0][1] * z[0][1] - z[1][0] * z[1][0]
        self.assertEqual(gap, Quadratic(Fraction(-1, 128)))


if __name__ == "__main__":
    unittest.main()
