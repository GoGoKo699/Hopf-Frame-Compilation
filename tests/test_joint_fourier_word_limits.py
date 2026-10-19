"""Bounded exact receipts for two scoped joint-word screens.

The Walsh check is exhaustive only at three bits. Native fixtures retain
literal phases, all occupied columns and actual inverses. They do not
assert an unrestricted native-depth obstruction or a faster source.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _ccz


def _walsh_numerators(values):
    return [sum(value * (-1 if (t & w).bit_count() & 1 else 1)
                for w, value in enumerate(values))
            for t in range(len(values))]


def _degree(table):
    coefficients = list(table)
    width = (len(table) - 1).bit_length()
    for bit in range(width):
        for mask in range(len(table)):
            if mask & (1 << bit):
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return max((mask.bit_count() for mask, value in enumerate(coefficients)
                if value), default=0)


class JointFourierWordLimitsTests(unittest.TestCase):
    def test_all_three_bit_tables_and_a_cubic_exact_closure(self):
        size = 8
        for truth in range(1 << size):
            table = [(truth >> w) & 1 for w in range(size)]
            fourier = _walsh_numerators([1 - 2 * bit for bit in table])
            support = [t for t in range(size) if fourier[t]]
            degree = max(1, _degree(table))
            self.assertEqual(sum(value * value for value in fourier), size * size)
            for h in (0, 1):
                mismatch = [sum(table[w ^ t] ^ table[w] ^ h for w in range(size))
                            for t in range(size)]
                residual_squared = Fraction(
                    4 * sum(fourier[t] ** 2 * mismatch[t] for t in range(size)),
                    size ** 3,
                )
                closes = all((s & t).bit_count() % 2 == h
                             for s in support for t in support)
                self.assertEqual(residual_squared == 0, closes)
                if closes:
                    self.assertEqual(h, table[0] ^ table[size - 1])
                else:
                    self.assertGreaterEqual(residual_squared,
                                            Fraction(2) ** (4 - 2 * degree))

        # A nonquadratic phase really can close; the theorem restricts
        # its scalar sign, not the internal phase to degree two.
        cubic = [((w ^ (w >> 1)) & 1)
                 & (((w >> 2) ^ (w >> 3)) & 1)
                 & (((w >> 4) ^ (w >> 5)) & 1) for w in range(64)]
        self.assertEqual(_degree(cubic), 3)
        support = [t for t, value in enumerate(
            _walsh_numerators([1 - 2 * bit for bit in cubic])) if value]
        self.assertTrue(all((s & t).bit_count() % 2 == 0
                            for s in support for t in support))

    def test_native_fourth_words_keep_the_literal_sign_and_inverse(self):
        # Chronological D,H realizes HD. A linear one-bit phase gives -I.
        linear = ([('S', 0), ('S', 0), ('H', 0)]) * 4
        self.assertEqual(native_columns(1, linear), scale(eye(2), -CO))
        self.assertEqual(native_columns(1, _inverse(linear)), scale(eye(2), -CO))

        cubic = _ccz((0, 1, 2))
        word = (cubic + [('H', wire) for wire in range(3)]) * 4
        matrix = native_columns(3, word)
        self.assertEqual(native_columns(3, _inverse(word)), adjoint(matrix))
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 28)
        for sign in (-CO, CO):
            residual = add(matrix, scale(eye(8), -sign))
            squared = sum(entry.abs2() for row in residual for entry in row) / 8
            self.assertGreaterEqual(squared, Fraction(1, 4))

    def test_compressed_spectrum_and_allowed_sqrt_two_normalization(self):
        cin = [('H', 0), ('H', 1), ('S', 2), ('CX', 1, 0), ('CX', 0, 2)]
        cout = [('H', 2), ('S', 0), ('CX', 2, 0), ('H', 1)]
        word = cin + _ccz((0, 1, 2)) + cout
        unitary = native_columns(3, word)
        self.assertEqual(native_columns(3, _inverse(word)), adjoint(unitary))
        accepted = [[unitary[row][column] for column in range(0, 8, 2)]
                    for row in range(0, 8, 2)]
        gram = multiply(adjoint(accepted), accepted)
        polynomial = multiply(gram, multiply(
            add(gram, scale(eye(4), -CO / 2)), add(gram, scale(eye(4), -CO))))
        self.assertEqual(polynomial, zeros(4))
        self.assertNotEqual(gram, scale(eye(4), CO / 4))

        native_band = cin + [('T', 0), ('TDG', 1), ('T', 2)] + cout
        unitary = native_columns(3, native_band)
        self.assertEqual(native_columns(3, _inverse(native_band)), adjoint(unitary))
        accepted = [[unitary[row][column] for column in range(0, 8, 2)]
                    for row in range(0, 8, 2)]
        gram = multiply(adjoint(accepted), accepted)
        # The extra two allowed squared singular values are the roots
        # of z^2-z+1/8; this check remains exact over Q(sqrt(2),i).
        quadratic = add(add(multiply(gram, gram), scale(gram, -CO)),
                        scale(eye(4), CO / 8))
        polynomial = multiply(gram, multiply(
            add(gram, scale(eye(4), -CO / 2)), multiply(
                add(gram, scale(eye(4), -CO)), quadratic)))
        self.assertEqual(polynomial, zeros(4))
        self.assertNotEqual(gram, scale(eye(4), CO / 4))

        # A nontrivial flat one-band source has both active classes.
        # The Bell Clifford transports Z_a,Z_t to X_a X_t,Z_a Z_t.
        bell = [('H', 0), ('CX', 0, 1)]
        source = _inverse(bell) + [('T', 0), ('T', 1)] + bell + [('H', 0)]
        unitary = native_columns(2, source)
        accepted = [[unitary[row][column] for column in (0, 2)] for row in (0, 2)]
        logical = native_columns(1, [('T', 0), ('H', 0), ('T', 0), ('H', 0)])
        inverse_sqrt_two = native_columns(1, [('H', 0)])[0][0]
        self.assertEqual(accepted, scale(logical, inverse_sqrt_two))
        self.assertEqual(multiply(adjoint(accepted), accepted), scale(eye(2), CO / 2))
        self.assertEqual(native_columns(2, _inverse(source)), adjoint(unitary))
        omega = native_columns(1, [('T', 0)])[1][1]
        self.assertEqual(native_columns(1, [('H', 0), ('S', 0)] * 3),
                         scale(eye(2), omega))

        # The screen cannot exclude every constant normalization.
        # One signal H supplies CCZ/sqrt(2), and the two-S recovery
        # returns that signal with the literal common phase exp(i*pi/4).
        source = [('H', 0)] + _ccz((1, 2, 3))
        u = native_columns(4, source)
        accepted = [[u[row][column] for column in range(0, 16, 2)]
                    for row in range(0, 16, 2)]
        self.assertEqual(multiply(adjoint(accepted), accepted),
                         scale(eye(8), CO / 2))
        recovered = source + [('S', 0)] + _inverse(source) + [('S', 0)] + source
        actual = native_columns(4, recovered)
        for column in range(0, 16, 2):
            phase = -omega if column == 14 else omega
            for row in range(16):
                self.assertEqual(actual[row][column], phase if row == column else CZ)
        self.assertEqual(native_columns(4, _inverse(recovered)), adjoint(actual))


if __name__ == '__main__':
    unittest.main()
