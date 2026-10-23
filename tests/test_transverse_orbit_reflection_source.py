"""Exact guards for broad orbit dressings and a paid transverse trial.

The native p=1 word checks arbitrary address and helper inputs.  The
small exact fibers challenge the analytic all-size obstruction; neither
these fixtures nor the normalization screen is a general depth bound.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_boolean_reflection_gram_source import _cs
from tests.test_half_bilinear_gram_source import _private_toffoli
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts
from tests.test_signal_interleaved_kernel_source import _ct_word, _root_powers


HALF = CQ2(Fraction(1, 2))
QUARTER = CQ2(Fraction(1, 4))


def _block_diagonal(first, second):
    size = len(first)
    result = zeros(2 * size)
    for row in range(size):
        for column in range(size):
            result[row][column] = first[row][column]
            result[size + row][size + column] = second[row][column]
    return result


def _real(matrix):
    return scale(add(matrix, adjoint(matrix)), HALF)


def _gram(square):
    return scale(add(eye(len(square)), _real(square)), QUARTER)


def _p1_parts(address):
    # w coordinates z1,t0,t1; the orbit bit l is absent.
    powers, diagonal = _root_powers(), zeros(8)
    for value in range(8):
        z1, t0 = value & 1, (value >> 1) & 1
        diagonal[value][value] = powers[((t0 ^ (address * z1)) - t0) % 8]
    mixer = native_columns(3, [('CX', 2, 0), ('H', 1)])
    e = multiply(multiply(adjoint(mixer), diagonal), mixer)
    a = multiply(multiply(multiply(adjoint(e),
                                   native_columns(3, [('S', 1)])), e), diagonal)
    return diagonal, e, a


def _p1_native_a():
    # x,z1,t0,t1,g.  The private phase helper g is never initialized.
    query = _private_toffoli(0, 1, 2, 4)
    diagonal = [('TDG', 2)] + query + [('T', 2)] + _inverse(query)
    mixer = [('CX', 3, 1), ('H', 2)]
    e = mixer + diagonal + _inverse(mixer)
    return diagonal + e + [('S', 2)] + _inverse(e)


def _p1_native_control_a():
    # l,z1,t0,x,t1,d,e,f,g.  Four returned CT helpers suffice serially.
    mapping = {0: 0, 1: 2, 2: 5, 3: 6, 4: 7, 5: 8}
    ct = [tuple([gate[0]] + [mapping[wire] for wire in gate[1:]])
          for gate in _ct_word()]
    query = _private_toffoli(3, 1, 2, 8)
    diagonal = _inverse(ct) + query + ct + _inverse(query)
    mixer = [('CX', 4, 1), ('H', 2)]
    e = mixer + diagonal + _inverse(mixer)
    return diagonal + e + _cs(0, 2) + _inverse(e)


class TransverseOrbitReflectionSourceTests(unittest.TestCase):
    def test_broad_noncommuting_dressing_normalization_and_transverse_square(self):
        # No commutation with a target, or between the two branches, is used.
        v0 = native_columns(1, [('H', 0), ('S', 0)])
        v1 = native_columns(1, [('T', 0), ('H', 0)])
        self.assertNotEqual(multiply(v0, v1), multiply(v1, v0))
        a = multiply(adjoint(v0), v1)
        v = _block_diagonal(v0, v1)
        x = native_columns(2, [('X', 1)])
        for anchor, weight in (([('H', 1)], HALF), ([('X', 1)], CO)):
            r = multiply(multiply(adjoint(v), native_columns(2, anchor)), v)
            self.assertEqual(r, adjoint(r))
            self.assertEqual(multiply(r, r), eye(4))
            rs = multiply(r, x)
            square = multiply(rs, rs)
            bank_gram = scale(add(eye(2), _real(multiply(a, a))),
                              weight * QUARTER)
            self.assertEqual(_gram(square), _block_diagonal(bank_gram, bank_gram))
            self.assertEqual(bank_gram,
                             scale(multiply(_real(a), _real(a)), weight * HALF))
            if anchor == [('X', 1)]:
                self.assertEqual(square, _block_diagonal(multiply(a, a),
                                                        multiply(adjoint(a), adjoint(a))))

    def test_literal_p1_a_keeps_coherent_address_helper_and_actual_reverse(self):
        word = _p1_native_a()
        self.assertEqual(_counts(word), (54, 69))
        actual = native_columns(5, word)
        self.assertEqual(native_columns(5, _inverse(word)), adjoint(actual))
        parts = [_p1_parts(address)[2] for address in (0, 1)]
        for column in range(32):
            address, bank, helper = column & 1, (column >> 1) & 7, column >> 4
            for row in range(32):
                expected = (parts[address][(row >> 1) & 7][bank]
                            if (row & 1) == address and (row >> 4) == helper else CZ)
                self.assertEqual(actual[row][column], expected)
        # Literal composition is counted separately. Its helper-return
        # primitives are checked by the retained controlled-D native fixture.
        controlled = _p1_native_control_a()
        self.assertEqual(_counts(controlled), (471, 598))
        self.assertEqual(_counts(_inverse(controlled)), (471, 598))

    def test_full_p1_square_base_sector_and_broken_blind_tail_symmetry(self):
        x_orbit = native_columns(4, [('X', 3)])
        blind = native_columns(3, [('X', 2)])
        for address in (0, 1):
            _, e, a = _p1_parts(address)
            v = _block_diagonal(eye(8), a)
            r = multiply(multiply(adjoint(v), x_orbit), v)
            self.assertEqual(r, adjoint(r))
            self.assertEqual(multiply(r, r), eye(16))
            rs = multiply(r, x_orbit)
            square = multiply(rs, rs)
            expected = _block_diagonal(multiply(a, a), multiply(adjoint(a), adjoint(a)))
            self.assertEqual(square, expected)
            charge_word = [('Z', 1)] + ([('Z', 0), ('Z', 2)] if address else [])
            charge = native_columns(3, charge_word)
            target = _block_diagonal(charge, charge)
            if not address:
                self.assertEqual(square, target)
            else:
                self.assertNotEqual(multiply(e, blind), multiply(blind, e))
                self.assertNotEqual(multiply(a, blind), multiply(blind, a))
                self.assertTrue(any(a[row][column] != CZ
                                    for row in range(8) for column in range(8) if row != column))
                self.assertEqual(_gram(square)[4][4], QUARTER)
                self.assertEqual(_gram(target)[4][4], CZ)
                delta = [square[row][4] - target[row][4] for row in range(16)]
                self.assertEqual(sum(value.conj() * value for value in delta), CQ2(2))

    def test_exact_two_dimensional_tail_fiber_at_three_sizes(self):
        powers = _root_powers()
        h = native_columns(1, [('H', 0)])
        phase = [[powers[1], CZ], [CZ, powers[7]]]
        expected_e = multiply(multiply(h, phase), h)
        minus_y = [[CZ, CI], [-CI, CZ]]
        minus_z = [[-CO, CZ], [CZ, CO]]
        x_bit = native_columns(1, [('X', 0)])
        s_bit = native_columns(1, [('S', 0)])
        self.assertEqual(multiply(multiply(adjoint(s_bit), x_bit), s_bit), minus_y)
        for p in (1, 2, 3):
            n = p * (1 << p)
            tail = 1 << (n - p)
            self.assertNotEqual(tail, 1)
            exponents = []
            for b in (0, 1):
                z, t = 0, tail | b
                # B's fan-in omits b and a=z_p, and includes the tail.
                bank_parity = ((z >> 1).bit_count() + (t & ~1).bit_count()) & 1
                mixed_z = z ^ (bank_parity << p)
                self.assertEqual(mixed_z, 1 << p)
                shifted = t ^ (mixed_z & ~1)
                exponent = (shifted ^ (mixed_z >> p)).bit_count() - shifted.bit_count()
                exponents.append(exponent)
                self.assertEqual(exponent, 1 - 2 * b)
            e = multiply(multiply(h, [[powers[exponents[0] % 8], CZ],
                                      [CZ, powers[exponents[1] % 8]]]), h)
            self.assertEqual(e, expected_e)
            a = multiply(multiply(adjoint(e), s_bit), e)
            self.assertEqual(multiply(a, a), minus_y)
            self.assertEqual(multiply(multiply(adjoint(a), x_bit), a), minus_z)
            defect = add(multiply(a, a), scale(minus_z, -CO))
            self.assertEqual(multiply(adjoint(defect), defect), scale(eye(2), CQ2(2)))
            gram_defect = add(_gram(multiply(a, a)), scale(_gram(minus_z), -CO))
            self.assertEqual(multiply(gram_defect, gram_defect),
                             scale(eye(2), CQ2(Fraction(1, 8))))


if __name__ == '__main__':
    unittest.main()
