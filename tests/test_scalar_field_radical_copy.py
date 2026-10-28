"""Exact fixtures for the scalar-field radical copy screen.

The linear kernel calculation imposes all bilinear equivariance constraints;
it does not assume associativity, nilpotence, or a sparse product grammar.
"""
from itertools import product
import unittest


def _gf_mul(a, b, degree, polynomial):
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a >> degree:
            a ^= polynomial
    return value


def _row_basis(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return pivots


def _nullspace(rows, dimension):
    pivots = _row_basis(rows)
    answer = []
    for free in range(dimension):
        if free in pivots:
            continue
        vector = 1 << free
        for pivot in sorted(pivots):
            if (pivots[pivot] & vector).bit_count() & 1:
                vector |= 1 << pivot
        answer.append(vector)
    return answer, len(pivots)


def _equivariant_multiplications(degree, polynomial, banks, helpers, scalars):
    bank_bits = degree * banks
    dimension = 2 * bank_bits + helpers
    variables = dimension ** 3
    mask = (1 << degree) - 1

    def query(vector, scalar):
        result = vector
        for bank in range(banks):
            u = (vector >> (bank * degree)) & mask
            result ^= _gf_mul(scalar, u, degree, polynomial) << (
                bank_bits + bank * degree
            )
        return result

    def index(left, right, output):
        return (left * dimension + right) * dimension + output

    def product_expression(left, right, output):
        expression = 0
        for i in range(dimension):
            if (left >> i) & 1:
                for j in range(dimension):
                    if (right >> j) & 1:
                        expression ^= 1 << index(i, j, output)
        return expression

    constraints = []
    for scalar, i, j, output in product(
        scalars, range(dimension), range(dimension), range(dimension)
    ):
        row = product_expression(query(1 << i, scalar), query(1 << j, scalar), output)
        for bit in range(dimension):
            if (query(1 << bit, scalar) >> output) & 1:
                row ^= 1 << index(i, j, bit)
        constraints.append(row)
    basis, rank = _nullspace(constraints, variables)
    return dimension, bank_bits, constraints, basis, rank, index


class ScalarFieldRadicalCopyTests(unittest.TestCase):
    def test_complete_equivariant_product_spaces_have_square_zero_image(self):
        # Each tuple fixes the entire linear equation system, not a candidate
        # product list. The last two cases have, respectively, a noncanonical
        # address embedding and a proper ratio subfield GF(4) inside GF(16).
        cases = (
            (2, 0b111, 1, 0, (1, 2), 18),
            (2, 0b111, 1, 1, (1, 2), 43),
            (2, 0b111, 2, 0, (1, 2), 144),
            (3, 0b1011, 1, 1, (1, 2), 103),
            (4, 0b10011, 1, 0, (1, 8, 10), 108),
            (4, 0b10011, 1, 0, (1, 6), 144),
        )
        for degree, polynomial, banks, helpers, scalars, expected in cases:
            with self.subTest(degree=degree, banks=banks, helpers=helpers):
                n, bank, constraints, basis, rank, index = _equivariant_multiplications(
                    degree, polynomial, banks, helpers, scalars
                )
                self.assertEqual(len(basis), expected)
                self.assertEqual(rank + len(basis), n ** 3)
                self.assertEqual(len(_row_basis(basis)), len(basis))
                output_mask = (1 << n) - 1
                for multiplication in basis:
                    self.assertTrue(all(
                        ((row & multiplication).bit_count() & 1) == 0
                        for row in constraints
                    ))
                    for left, right in product(range(bank, 2 * bank), repeat=2):
                        self.assertEqual(
                            (multiplication >> index(left, right, 0)) & output_mask,
                            0,
                        )

    def test_one_address_hypothesis_is_essential(self):
        # Nilpotent algebra with basis a,d,b,c, products ab=d+c,
        # ba=d, b^2=c, and every product involving d or c zero.
        # Q(a)=a+b and Q(d)=d+c: U=(a,d), W=(b,c).
        def multiply(left, right):
            la, lb = left & 1, (left >> 2) & 1
            ra, rb = right & 1, (right >> 2) & 1
            return ((la * rb) ^ (lb * ra)) << 1 | (
                (la * rb) ^ (lb * rb)
            ) << 3

        def query(vector):
            return vector ^ ((vector & 1) << 2) ^ ((vector & 2) << 2)

        basis = (1, 2, 4, 8)
        for a, b, c in product(basis, repeat=3):
            self.assertEqual(multiply(multiply(a, b), c), 0)
            self.assertEqual(multiply(a, multiply(b, c)), 0)
        for a, b in product(range(16), repeat=2):
            self.assertEqual(query(multiply(a, b)), multiply(query(a), query(b)))
        for a in range(16):
            self.assertEqual(query(query(a)), a)
        self.assertEqual(multiply(1, 1), 0)
        self.assertEqual(multiply(4, 4), 8)  # W^2 is nonzero for p=1.

    def test_regular_affine_four_word_on_every_column(self):
        # Strictly upper-triangular 3x3 matrices over GF(4), encoded as
        # (u,w,h) = entries (23,13,12). Product has only entry h*u' at 13.
        def field(a, b):
            return _gf_mul(a, b, 2, 0b111)

        def multiply(left, right):
            h = (left >> 4) & 3
            u = right & 3
            return field(h, u) << 2

        def query(vector, scalar):
            return vector ^ (field(scalar, vector & 3) << 2)

        def affine(a, vector):
            return a ^ vector ^ multiply(a, vector)

        # Complete multiplication equivariance and generator conjugation.
        for scalar, left, right in product(range(4), range(64), range(64)):
            self.assertEqual(
                query(multiply(left, right), scalar),
                multiply(query(left, scalar), query(right, scalar)),
            )
        seed = 1
        generators = [query(seed, scalar) for scalar in range(4)]
        for scalar, vector in product(range(4), range(64)):
            self.assertEqual(affine(generators[scalar], affine(generators[scalar], vector)), vector)
            self.assertEqual(
                query(affine(seed, query(vector, scalar)), scalar),
                affine(generators[scalar], vector),
            )
        # All nonidentity units act freely, not just the four generators.
        for a in range(1, 64):
            self.assertTrue(all(affine(a, vector) != vector for vector in range(64)))
        planes = 0
        for origin, first, second in product(range(4), repeat=3):
            if len({0, first, second, first ^ second}) != 4:
                continue
            addresses = (origin, origin ^ first, origin ^ second, origin ^ first ^ second)
            for vector in range(64):
                result = vector
                for address in addresses:
                    result = affine(generators[address], result)
                self.assertEqual(result, vector)
            planes += 1
        self.assertEqual(planes, 24)


if __name__ == '__main__':
    unittest.main()
