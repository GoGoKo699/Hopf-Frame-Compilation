"""Exact guards for paid signal interleaving and its failed lookup sectors."""
from fractions import Fraction
import unittest

import numpy as np

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, native_columns,
    scale, zeros,
)
from tests.test_bent_gram_source import ry_t
from tests.test_half_bilinear_gram_source import (
    _binary_image, _ccs_word, _private_toffoli, _walsh_integer,
)
from tests.test_majorana_walsh_query import _inverse
from tests.test_nilpotent_pencil_observable_source import _field_matrix, _power
from tests.test_phase_signal_query_limits import _counts
from tests.test_t_wall_xor_gram_source import _p1_mixer


def _ct_word():
    query = _private_toffoli(0, 1, 2, 5)
    return query + [('T', 2)] + query + [('TDG', 2)] + _ccs_word()


def _root_powers():
    omega = CQ2(Q2(0, Fraction(1, 2)), Q2(0, Fraction(1, 2)))
    powers = [CO]
    for _ in range(7):
        powers.append(powers[-1] * omega)
    return powers


def _p1_parts(address):
    powers, diagonal = _root_powers(), zeros(16)
    for value in range(16):
        z, t = value & 3, value >> 2
        shifted = t ^ (z & 2)
        exponent = (shifted ^ (address * (z >> 1))).bit_count() - shifted.bit_count()
        diagonal[value][value] = powers[exponent % 8]
    mixer = native_columns(4, _p1_mixer())
    e = multiply(multiply(adjoint(mixer), diagonal), mixer)
    f = multiply(native_columns(4, [('Z', 0), ('Z', 2)]), adjoint(diagonal))
    return e, f


def _control(matrix):
    size = len(matrix)
    result = zeros(2 * size)
    for i in range(size):
        result[i][i] = CO
        for j in range(size):
            result[size + i][size + j] = matrix[i][j]
    return result


def _one_clean_block(e, f):
    a, b = scale(e, CI), scale(f, CI)
    alpha = CO + CI * CQ2(Q2(-1, 1))
    beta = CQ2(Q2(1, 1)) + CI
    numerator = add(add(add(scale(eye(len(e)), alpha), scale(a, alpha.conj())),
                        scale(b, beta)), scale(multiply(b, a), -beta.conj()))
    return scale(numerator, (CO + CI) * CQ2(Fraction(1, 8)))


def _integer_receipt(p, family):
    """Exact accepted norm at x=1,z=e_p,t=0 in Z[omega].

    Two full Walsh walls give denominator 2^(2n). Coefficients are
    kept in int64; at p<=2 their bounds are far below 2^63. Convert
    to Python integers before the quadratic norm expressions.
    """
    m, n = 1 << p, p * (1 << p)
    size, dimension = 1 << n, 1 << (2 * n)
    basis = ((3, 2) if p == 1 else
             _field_matrix([[_power(u, j) for u in range(m)] for j in range(m)]))
    images = np.array([_binary_image(basis, z) for z in range(size)], dtype=np.int64)
    labels = np.arange(dimension, dtype=np.int64)
    z, t = labels & (size - 1), labels >> n
    permutation = (images[z] ^ t) | (z << n)
    weights = np.array([value.bit_count() for value in range(size)], dtype=np.int64)
    shifted = t ^ (z & ~1)
    h = weights[shifted ^ (z >> p)] - weights[shifted]
    basis_v = 1 << p
    flip = 0 if family == 'one' else (1 if family == 'two' else 1 | (1 << p))
    values = np.zeros((4, dimension), dtype=np.int64)
    values[0, basis_v ^ flip] = 1

    def phase(old, powers):
        new = np.zeros_like(old)
        for power in range(8):
            mask = powers == power
            for degree in range(4):
                target = (degree + power) % 8
                new[target % 4, mask] = ((-1 if target >= 4 else 1)
                                        * old[degree, mask])
        return new

    def polynomial(old, coefficients):
        new = np.zeros_like(old)
        for degree in range(4):
            for shift, coefficient in enumerate(coefficients):
                target = (degree + shift) % 8
                new[target % 4] += ((-1 if target >= 4 else 1)
                                    * coefficient * old[degree])
        return new

    def norm(old, denominator):
        rational = radical = 0
        for column in old.T:
            a, b, c, d = map(int, column)
            rational += a*a + b*b + c*c + d*d
            radical += a*b + b*c + c*d - a*d
        return Fraction(rational, denominator**2), Fraction(radical, denominator**2)

    routed = np.zeros_like(values)
    routed[:, permutation] = values
    values = routed
    for row in values:
        _walsh_integer(row)
    values = phase(values, h % 8)
    for row in values:
        _walsh_integer(row)
    values = values[:, permutation]
    denominator = 1 << (2 * n)
    assert norm(values, denominator) == (1, 0)
    fpower = (-h + 4 * ((z ^ t) & 1)) % 8
    fvalues = phase(values, fpower)
    basis_column = np.zeros_like(values)
    basis_column[0, basis_v] = denominator
    if family != 'one':
        accepted = (basis_column + phase(basis_column, fpower)
                    + values - fvalues)
        return norm(accepted, 4 * denominator)
    # alpha I + i conj(alpha) E + i beta F + conj(beta) F E,
    # followed by the literal prefactor (1+i)/8.
    accepted = (polynomial(basis_column, [1, 1, -1, 1])
                + polynomial(polynomial(values, [1, -1, 1, -1]), [0, 0, 1, 0])
                + polynomial(polynomial(phase(basis_column, fpower),
                                        [1, 1, 1, -1]), [0, 0, 1, 0])
                + polynomial(fvalues, [1, 1, -1, -1]))
    accepted = polynomial(accepted, [1, 0, 1, 0])
    return norm(accepted, 8 * denominator)


class SignalInterleavedKernelSourceTests(unittest.TestCase):
    def test_literal_controlled_t_returns_all_four_dirty_helpers(self):
        word = _ct_word()
        self.assertEqual(_counts(word), (70, 88))
        self.assertEqual(_counts(_inverse(word)), (70, 88))
        self.assertTrue(all(gate[1] in (2, 3, 4, 5) for gate in word
                            if gate[0] in ('T', 'TDG')))
        actual = native_columns(6, word)
        for column in range(64):
            scalar = _root_powers()[1] if column & 3 == 3 else CO
            for row in range(64):
                self.assertEqual(actual[row][column], scalar if row == column else CZ)
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))

    def test_controlled_phase_keeps_coherent_address_and_every_dirty_column(self):
        # signal,z1,t0,address, then four private helpers. Unaffected
        # bank coordinates are literal identity factors in this fixture.
        mapping = {0: 0, 1: 2, 2: 4, 3: 5, 4: 6, 5: 7}
        ct = [tuple([gate[0]] + [mapping[wire] for wire in gate[1:]])
              for gate in _ct_word()]
        query = _private_toffoli(3, 1, 2, 7)
        word = _inverse(ct) + query + ct + _inverse(query)
        self.assertEqual(_counts(word), (156, 196))
        actual = native_columns(8, word)
        powers = _root_powers()
        for column in range(256):
            signal, z1, t0, address = [(column >> bit) & 1 for bit in range(4)]
            exponent = signal * ((t0 ^ (address * z1)) - t0)
            for row in range(256):
                self.assertEqual(actual[row][column],
                                 powers[exponent % 8] if row == column else CZ)

    def test_complete_signal_word_matches_polynomial_and_base_normalization(self):
        signal = 4
        rotation = native_columns(5, ry_t(signal))
        s_gate = native_columns(5, [('S', signal)])
        h_gate = native_columns(5, [('H', signal)])
        for address in (0, 1):
            e, f = _p1_parts(address)
            source = eye(32)
            for operation in (rotation, s_gate, h_gate, _control(scale(e, CI)),
                              h_gate, _control(scale(f, CI)), h_gate, rotation):
                source = multiply(operation, source)
            accepted = [row[:16] for row in source[:16]]
            expected = _one_clean_block(e, f)
            self.assertEqual(accepted, expected)
            gram = multiply(adjoint(accepted), accepted)
            if address == 0:
                for row in range(16):
                    for column in range(16):
                        equal = not ((column ^ (column >> 2)) & 1)
                        scalar = CQ2(Fraction(1, 2)) if row == column and equal else CZ
                        self.assertEqual(gram[row][column], scalar)
            else:
                self.assertEqual(gram[2][2], CQ2(Fraction(1, 2)))

    def test_exact_small_sector_receipts(self):
        expected = {
            (1, 'two'): (Fraction(1, 4), Fraction(1, 8)),
            (1, 'paired'): (Fraction(3, 8), Fraction(1, 16)),
            (1, 'one'): (Fraction(1, 2), Fraction(0)),
            (2, 'two'): (Fraction(7991, 32768), Fraction(1391, 16384)),
            (2, 'paired'): (Fraction(9305, 32768), Fraction(1, 16)),
            (2, 'one'): (Fraction(1197, 2048), Fraction(-155, 8192)),
        }
        for inputs, receipt in expected.items():
            self.assertEqual(_integer_receipt(*inputs), receipt)


if __name__ == '__main__':
    unittest.main()
