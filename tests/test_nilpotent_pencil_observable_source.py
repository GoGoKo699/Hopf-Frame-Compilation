"""Exact matrix fixtures for the nilpotent-pencil source and its interfaces.

Matrices act on arbitrary binary bank labels, so their equalities cover
all carrier and second-bank inputs without Hilbert-space enumeration.
This does not emit the unpaid in-place source or a native schedule.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import CI, CO, CQ2, CZ, Q2, adjoint, native_columns
from tests.test_majorana_walsh_query import _inverse as _word_inverse


def _mul(a, b):
    """GF(4)=F_2[t]/(t^2+t+1), binary basis (1,t)."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def _power(a, exponent):
    result = 1
    for _ in range(exponent):
        result = _mul(result, a)
    return result


def _identity(n):
    return tuple(1 << i for i in range(n))


def _product(a, b):
    """Binary matrices stored as integer row masks; return a b."""
    result = []
    for row in a:
        value = 0
        for i, source in enumerate(b):
            if row & (1 << i):
                value ^= source
        result.append(value)
    return tuple(result)


def _transpose(a):
    return tuple(sum(((a[j] >> i) & 1) << j for j in range(len(a)))
                 for i in range(len(a)))


def _inverse(a):
    n = len(a)
    left, right = list(a), list(_identity(n))
    for column in range(n):
        pivot = next(i for i in range(column, n)
                     if left[i] & (1 << column))
        left[column], left[pivot] = left[pivot], left[column]
        right[column], right[pivot] = right[pivot], right[column]
        for i in range(n):
            if i != column and left[i] & (1 << column):
                left[i] ^= left[column]
                right[i] ^= right[column]
    if tuple(left) != _identity(n):
        raise AssertionError('Binary elimination did not reach identity')
    return tuple(right)


def _field_matrix(a):
    """Expand a GF(4) matrix into binary basis (1,t) in every block."""
    n = len(a)
    rows = [0] * (2 * n)
    for i, row in enumerate(a):
        for j, coefficient in enumerate(row):
            for bit in range(2):
                image = _mul(coefficient, 1 << bit)
                for output_bit in range(2):
                    rows[2 * i + output_bit] |= (
                        ((image >> output_bit) & 1) << (2 * j + bit))
    return tuple(rows)


def _double(a):
    n = len(a)
    return tuple(a) + tuple(row << n for row in a)


def _shear(a):
    n = len(a)
    return _identity(n) + tuple(row | (1 << (n + i))
                                for i, row in enumerate(a))


class NilpotentPencilObservableSourceTests(unittest.TestCase):
    def test_paired_pauli_sum_is_half_rank_swap_not_a_thin_code(self):
        # The numerator of P(a,b) tensor P(a,b) on |c,d> is the
        # integer sign (-1)^(a.b+b.(c+d)), with output |c+a,d+a>.
        # Summing exact integers keeps every phase and every bank column.
        for n in (1, 2, 3):
            size = 1 << n
            for nu in range(1, size):
                nonzero_columns = 0
                for c in range(size):
                    for d in range(size):
                        numerator = {}
                        for a in range(size):
                            if (a & nu).bit_count() % 2:
                                continue
                            output = (c ^ a, d ^ a)
                            for b in range(size):
                                sign = (-1) ** ((b & (a ^ c ^ d)).bit_count() % 2)
                                numerator[output] = numerator.get(output, 0) + sign
                        numerator = {row: value for row, value in numerator.items()
                                     if value}
                        same_sector = not ((nu & (c ^ d)).bit_count() % 2)
                        expected = {(d, c): size} if same_sector else {}
                        self.assertEqual(numerator, expected)
                        nonzero_columns += bool(numerator)
                self.assertEqual(nonzero_columns, size * size // 2)

    def test_partial_gram_phase_s_readout_returns_arbitrary_private_fiber(self):
        # Logical banks are c=0,d=1, private h=2, clean signal a=3.
        # The source acts as H_h tensor H_a on c=d, and H_h tensor X_a
        # otherwise. Its accepted Gram is Pi_equal/2; its nontrivial
        # private Hadamard fiber need not be returned by the source.
        def cz(control, target):
            return [('H', target), ('CX', control, target), ('H', target)]

        rotation = [('S', 3)] * 3 + [('H', 3), ('T', 3),
                                     ('H', 3), ('S', 3)]
        controlled_h = _word_inverse(rotation) + cz(0, 3) + rotation
        source = ([('H', 2), ('CX', 1, 0), ('H', 3)]
                  + controlled_h + [('CX', 0, 3), ('CX', 1, 0)])
        source_matrix = native_columns(4, source)
        self.assertEqual(native_columns(4, _word_inverse(source)), adjoint(source_matrix))
        half_root = Q2(0, Fraction(1, 2))
        omega = CQ2(half_root, half_root)
        half = CQ2(Q2(Fraction(1, 2)))
        for column in range(8):
            equal = (column & 1) == ((column >> 1) & 1)
            for row in range(8):
                sign = -1 if (row & 4 and column & 4) else 1
                expected = half * sign if equal and row & 3 == column & 3 else CZ
                self.assertEqual(source_matrix[row][column], expected)

        recovered = source + [('S', 3)] + _word_inverse(source) + [('S', 3)] + source
        recovered_matrix = native_columns(4, recovered)
        self.assertEqual(native_columns(4, _word_inverse(recovered)),
                         adjoint(recovered_matrix))
        # Exact clean-code pass phase and whole rejected amplitude; the
        # arbitrary h fiber remains H until the actual inverse readout.
        for column in range(8):
            equal = (column & 1) == ((column >> 1) & 1)
            for row in range(16):
                sign = -1 if (row & 4 and column & 4) else 1
                scalar = omega if equal else CI
                signal = 0 if equal else 8
                expected = (scalar * CQ2(half_root) * sign
                            if row & 11 == (column & 3) | signal else CZ)
                self.assertEqual(recovered_matrix[row][column], expected)

        # V=W^dagger Z_a W is exactly Hermitian. The emitted actual
        # inverse is the same word, not an independently compiled target.
        readout = recovered + [('S', 3), ('S', 3)] + _word_inverse(recovered)
        readout_matrix = native_columns(4, readout)
        self.assertEqual(readout_matrix, adjoint(readout_matrix))
        self.assertEqual(native_columns(4, _word_inverse(readout)), readout_matrix)
        for column in range(16):
            parity = ((column >> 0) ^ (column >> 1) ^ (column >> 3)) & 1
            for row in range(16):
                expected = (CO if not parity else -CO) if row == column else CZ
                self.assertEqual(readout_matrix[row][column], expected)

    def test_gf4_complete_bank_and_quadratic_xor_matrices(self):
        p, m, n = 2, 4, 8
        # Columns are v_u, then binary expansion orders beta_0 v_u,
        # beta_1 v_u. Physical output Y_u is therefore binary slot 2u.
        vandermonde = [[_power(u, j) for u in range(m)]
                       for j in range(m)]
        basis = _field_matrix(vandermonde)
        basis_inverse = _inverse(basis)
        self.assertEqual(_product(basis, basis_inverse), _identity(n))
        self.assertEqual(_product(basis_inverse, basis), _identity(n))
        to_dual = _transpose(basis_inverse)
        from_dual = _transpose(basis)

        shift = [[int(i == j + 1) for j in range(m)] for i in range(m)]
        star = [[0] * m for _ in range(m)]
        for u in range(1, m):
            reciprocal = _power(u, m - 2)
            star[0][u] = star[u][0] = star[u][u] = reciprocal
        self.assertEqual(_product(basis_inverse,
                                 _product(_field_matrix(shift), basis)),
                         _field_matrix(star))

        queries, labels = [], []
        for x in range(m):
            standard_field = [[int(i == j) ^ _mul(x, shift[i][j])
                               for j in range(m)] for i in range(m)]
            standard = _field_matrix(standard_field)
            pencil = _product(basis_inverse, _product(standard, basis))
            inverse = _inverse(pencil)
            self.assertEqual(_product(pencil, inverse), _identity(n))
            star_pencil = [[int(i == j) ^ _mul(x, star[i][j])
                            for j in range(m)] for i in range(m)]
            self.assertEqual(pencil, _field_matrix(star_pencil))

            # J=U_{N^{-T}}: its probe coordinate is Y_x on every
            # arbitrary carrier. Every other output/private coefficient
            # of this row is zero; no carrier is set to a constant.
            source_data = _transpose(inverse)
            self.assertEqual(source_data[0], 1 << (p * x))
            inverse_column = sum(((row >> 0) & 1) << i
                                 for i, row in enumerate(inverse))
            self.assertEqual(inverse_column, 1 << (p * x))
            labels.append(inverse_column)

            # Build the paid query separately from upper-neighbor field
            # multiplication and the identity COPY. The block matrix
            # includes BOTH arbitrary principal banks.
            multiplication = _field_matrix([[x]])
            sparse_transpose = list(_identity(n))
            for block in range(m - 1):
                for out_bit in range(p):
                    for in_bit in range(p):
                        if multiplication[in_bit] & (1 << out_bit):
                            sparse_transpose[p * block + out_bit] ^= (
                                1 << (p * (block + 1) + in_bit))
            self.assertEqual(tuple(sparse_transpose), _transpose(standard))
            implemented = _product(
                _double(from_dual),
                _product(_shear(tuple(sparse_transpose)), _double(to_dual)))
            expected = _shear(_transpose(pencil))
            self.assertEqual(implemented, expected)
            self.assertEqual(_product(implemented, implemented),
                             _identity(2 * n))
            queries.append(implemented)

            zeros = [u for u in range(1, m)
                     if 1 ^ _mul(x, _power(u, m - 2)) == 0]
            self.assertEqual(zeros, [] if x == 0 else [x])

        # The whole XOR query is affine in the address matrix entries,
        # while its target inverse-column labels have quadratic degree.
        self.assertEqual(queries[3], tuple(a ^ b ^ c for a, b, c in
                                          zip(queries[0], queries[1], queries[2])))
        self.assertNotEqual(labels[0] ^ labels[1] ^ labels[2] ^ labels[3], 0)


if __name__ == '__main__':
    unittest.main()
