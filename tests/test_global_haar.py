"""Exact finite certificates for the global Haar research chapter.

These checks verify the literal alternating witness, all signed path
witnesses at bounded sizes, and the absolute-matrix identities supporting
the weighted-norm proof. The dimension-uniform norm theorems are proved
in STRUCTURAL_COMPILATION_LIMITS.md, not inferred from these tests.
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


def signed_path_frame(n, leaf):
    """Literal full Hopf word with path angles zero or pi/2."""
    matrix = identity(1 << n)
    for depth in range(n):
        if leaf & (1 << (n - depth - 1)):
            prefix = leaf >> (n - depth)
            anchor = prefix << (n - depth)
            marker = anchor + (1 << (n - depth - 1))
            matrix[anchor], matrix[marker] = (
                [-entry for entry in matrix[marker]], matrix[anchor])
    return matrix


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

    def test_all_signed_path_witnesses_retain_literal_positive_root_amplitude(self):
        for n in (2, 3, 4):
            size = 1 << n
            e0 = [ONE] + [ZERO] * (size - 1)
            for leaf in range(size):
                with self.subTest(n=n, leaf=leaf):
                    frame = signed_path_frame(n, leaf)
                    self.assertEqual(
                        apply(frame, e0),
                        [ONE if j == leaf else ZERO for j in range(size)])
                    self.assertEqual(multiply(transpose(frame), frame), identity(size))
                    self.assertTrue(all(sum(x != ZERO for x in row) == 1
                                        for row in frame))

    def test_absolute_haar_gram_row_sums_and_connected_comparison_matrix(self):
        for n in (2, 3, 4):
            with self.subTest(n=n):
                q = haar(n)
                size = 1 << n
                root_entry = inverse_sqrt_power_of_two(n)
                self.assertEqual([row[0] for row in q], [root_entry] * size)
                # Haar entries are rational OR rational multiples of sqrt(2),
                # so absolute values are exact coefficientwise in this field.
                self.assertTrue(all(x.a * x.b == 0 for row in q for x in row))
                absolute = [[Quadratic(abs(x.a), abs(x.b)) for x in row]
                            for row in q]
                gram = multiply(absolute, transpose(absolute))
                self.assertEqual([sum(row, ZERO) for row in gram],
                                 [Quadratic(Fraction(n + 1))] * size)
                self.assertTrue(all(absolute[0][j] + absolute[j][0] != ZERO
                                    for j in range(size)))

    def test_free_five_mask_contradiction_has_exact_finite_margin(self):
        # At n=128 and epsilon<=1/2, eight Q occurrences would require
        # (2 sqrt(129))^9 >= 2^63. Both comparisons use integer arithmetic.
        self.assertLess(4 * 129, 32 ** 2)
        self.assertLess(32 ** 9, 2 ** 63)


if __name__ == "__main__":
    unittest.main()
