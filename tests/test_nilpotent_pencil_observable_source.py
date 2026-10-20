"""One exact GF(4) matrix fixture for the nilpotent-pencil source.

Matrices act on arbitrary binary bank labels, so their equalities cover
all carrier and second-bank inputs without Hilbert-space enumeration.
This does not emit the unpaid in-place source or a native schedule.
"""
import unittest


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
