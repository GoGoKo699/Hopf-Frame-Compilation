"""Small exact identities used by the fixed Clifford-sum source screen.

The source fixtures are deliberately nonflat. They check direct complete
compressions, physical normalization and actual inverses, not the uniform
gap theorem. Sparse indicator and exponent checks are bounded arithmetic.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_majorana_walsh_query import _inverse


def _conjugate(unitary, matrix):
    return multiply(multiply(unitary, matrix), adjoint(unitary))


def _normalized_hs_squared(matrix):
    return sum(entry.abs2() for row in matrix for entry in row) / len(matrix)


class FixedCliffordSourceGapTests(unittest.TestCase):
    def test_nonflat_complete_source_defect_and_literal_clifford_expansion(self):
        for clean in (1, 2):
            with self.subTest(clean=clean):
                logical_wire = clean
                width = clean + 1
                size = 1 << width
                sectors = 1 << clean
                # Each actual encoder has exactly one T band. All physical
                # columns are evaluated; only J's declared clean input is used.
                ua = ([('CX', logical_wire, a) for a in range(clean)]
                      + [('H', logical_wire), ('T', logical_wire)])
                ub = ([('CX', logical_wire, 0)]
                      + [('H', a) for a in range(clean)]
                      + [('H', logical_wire), ('T', logical_wire),
                         ('H', logical_wire)])
                source = ua + _inverse(ub)
                uam, ubm = native_columns(width, ua), native_columns(width, ub)
                unitary = native_columns(width, source)
                self.assertEqual(native_columns(width, _inverse(source)),
                                 adjoint(unitary))
                for word, matrix in ((ua, uam), (ub, ubm)):
                    self.assertEqual(native_columns(width, _inverse(word)),
                                     adjoint(matrix))

                code = [0, sectors]
                j = [[CO if row == column * sectors else CZ
                      for column in range(2)] for row in range(size)]
                pi = multiply(j, adjoint(j))
                pi_a, pi_b = _conjugate(uam, pi), _conjugate(ubm, pi)
                accepted = [[unitary[row][column] for column in code]
                            for row in code]
                alpha_squared = CO / 2
                gram = multiply(adjoint(accepted), accepted)
                self.assertNotEqual(gram, scale(eye(2), alpha_squared))
                self.assertNotEqual(gram, scale(eye(2), (gram[0][0] + gram[1][1]) / 2))

                physical_p = native_columns(width, [('Z', logical_wire)])
                logical_p = native_columns(1, [('Z', 0)])
                c, d = _conjugate(ubm, physical_p), _conjugate(uam, physical_p)
                direct = add(multiply(multiply(multiply(pi_a, pi_b), c), pi_a),
                             scale(multiply(d, pi_a), -alpha_squared))
                logical_defect = add(
                    multiply(multiply(adjoint(accepted), logical_p), accepted),
                    scale(logical_p, -alpha_squared))
                embedded = multiply(multiply(j, logical_defect), adjoint(j))
                self.assertEqual(direct, _conjugate(uam, embedded))
                self.assertEqual(_normalized_hs_squared(direct),
                                 _normalized_hs_squared(logical_defect) / sectors)
                self.assertNotEqual(direct, zeros(size))

                # Every stabilizer image below is an actual Hermitian
                # Clifford (one physical T band conjugating a Pauli). Keep
                # all ordered factors and literal signs in the expansion.
                images_a, images_b = [], []
                for mask in range(sectors):
                    parity = native_columns(width, [('Z', a) for a in range(clean)
                                                    if mask & (1 << a)])
                    images_a.append(_conjugate(uam, parity))
                    images_b.append(_conjugate(ubm, parity))
                expansion = zeros(size)
                count = 0
                for left, middle, right in product(range(sectors), repeat=3):
                    term = multiply(multiply(multiply(images_a[left],
                                                       images_b[middle]), c),
                                    images_a[right])
                    expansion = add(expansion, scale(term, CO / sectors ** 3))
                    count += 1
                for image_a in images_a:
                    expansion = add(expansion, scale(multiply(d, image_a),
                                                     -alpha_squared / sectors))
                    count += 1
                self.assertEqual(count, 2 ** (3 * clean) + 2 ** clean)
                self.assertEqual(expansion, direct)

    def test_sparse_indicator_probe_has_small_hs_but_full_operator_defect(self):
        # The bank-Z sign and every arbitrary dirty/reference multiplicity
        # cancel in the squared normalized trace. Only the address fraction
        # remains, so no exponentially wider indicator matrix is built.
        address_bits = 3
        rows = 1 << address_bits
        for selected in range(rows):
            ratios = [1 - 2 * (address == selected) for address in range(rows)]
            defects = [ratio - 1 for ratio in ratios]
            self.assertEqual(Fraction(sum(value * value for value in defects), rows),
                             Fraction(4, rows))
            self.assertEqual(max(abs(value) for value in defects), 2)

    def test_explicit_gap_exponents_satisfy_the_required_arithmetic_budgets(self):
        # Arithmetic verification only: the uniform separation theorem is
        # analytic. Keep huge gap powers symbolic as integer exponents.
        coefficient_bound = 2
        exponent = 0
        expected = [8, 50, 325, 2641, 26460, 317579, 4446174,
                    71138878, 1280499909]
        for terms, expected_exponent in zip(range(2, 11), expected):
            ceiling_log = (terms * coefficient_bound - 1).bit_length()
            threshold = 4 * exponent + 4 * ceiling_log + 2
            denominator_exponent = (terms - 1) * threshold // 2 + 1
            next_exponent = denominator_exponent + ceiling_log
            self.assertEqual(next_exponent, expected_exponent)
            self.assertEqual(threshold % 2, 0)
            # Cross-component error <= half the preceding squared gap.
            self.assertLessEqual((terms * coefficient_bound) ** 2,
                                 1 << (threshold // 2 - 2 * exponent - 1))
            # Connected-component algebraic norm denominator is paid.
            self.assertLessEqual(terms * coefficient_bound,
                                 1 << (next_exponent - denominator_exponent))
            self.assertGreaterEqual(next_exponent, exponent + 1)
            exponent = next_exponent
        # The direct k=1 defect is an integer ten-term sum after scaling by8.
        self.assertEqual(2 * exponent + 9, 2560999827)


if __name__ == '__main__':
    unittest.main()
