"""Exact finite identities for the bounded-clean two-band Gram screen.

These fixtures do not prove the dimension-independent rank--norm bound
or exclude a source at every specified moderate Gram tolerance.
"""
from fractions import Fraction
from functools import reduce
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_majorana_walsh_query import _inverse


def _product(*matrices):
    return reduce(multiply, matrices)


def _conjugate(unitary, matrix):
    return _product(unitary, matrix, adjoint(unitary))


def _subtract(left, right):
    return add(left, scale(right, -CO))


def _diagonal(entries):
    return [[value if i == j else CZ for j in range(len(entries))]
            for i, value in enumerate(entries)]


class TwoBandGramScreenTests(unittest.TestCase):
    def test_complete_bank_encoded_defect_and_clifford_sum(self):
        # Two clean inputs, one address bit, and both arbitrary bits of one
        # bank. This is the single-bank Clifford-conjugate paired target;
        # the paired lookup support is checked separately below.
        # The encoder networks move the address and mix it with the clean
        # inputs. Each displayed encoder has just one native T band.
        clean, width = 2, 5
        sectors, size, logical_size = 1 << clean, 1 << width, 8
        ua = [('CX', 2, 0), ('H', 0), ('H', 1),
              ('T', 0), ('T', 1), ('T', 2),
              ('H', 2), ('CX', 1, 4), ('CX', 3, 0)]
        ub = [('CX', 3, 1), ('H', 0), ('H', 2),
              ('T', 0), ('T', 2), ('T', 4),
              ('H', 4), ('CX', 0, 2), ('CX', 4, 1)]
        uam, ubm = native_columns(width, ua), native_columns(width, ub)
        for word, matrix in ((ua, uam), (ub, ubm)):
            self.assertEqual(native_columns(width, _inverse(word)), adjoint(matrix))
        source = ua + _inverse(ub)
        self.assertEqual(native_columns(width, _inverse(source)),
                         adjoint(native_columns(width, source)))

        j = [[CO if row == column * sectors else CZ
              for column in range(logical_size)] for row in range(size)]
        va = multiply(uam, j)
        pa = multiply(va, adjoint(va))
        pb = _product(ubm, j, adjoint(j), adjoint(ubm))
        physical_l = native_columns(width, [('X', 3)])
        logical_l = _product(adjoint(j), physical_l, j)
        d = _conjugate(uam, physical_l)
        accepted = _product(adjoint(j), adjoint(ubm), uam, j)
        gram = multiply(adjoint(accepted), accepted)
        defect = _product(pa, _subtract(pb, _conjugate(d, pb)), pa)
        self.assertEqual(defect, _product(
            va, _subtract(gram, _conjugate(logical_l, gram)), adjoint(va)))
        self.assertEqual(multiply(d, pa), multiply(pa, d))

        images_a, images_b = [], []
        for mask in range(sectors):
            parity = native_columns(width, [('Z', a) for a in range(clean)
                                            if mask & (1 << a)])
            images_a.append(_conjugate(uam, parity))
            images_b.append(_conjugate(ubm, parity))
        sum_a = reduce(add, images_a, zeros(size))
        sum_b = reduce(add, images_b, zeros(size))
        self.assertEqual(scale(sum_a, CO / sectors), pa)
        self.assertEqual(scale(sum_b, CO / sectors), pb)
        self.assertEqual(_subtract(images_b[0], _conjugate(d, images_b[0])),
                         zeros(size))
        middle = reduce(add, [_subtract(b, _conjugate(d, b))
                             for b in images_b[1:]], zeros(size))
        # Distributivity expands this exact identity into 4 * 6 * 4
        # literal Clifford summands; ordered products are not collapsed.
        self.assertEqual(2 * sectors ** 2 * (sectors - 1), 96)
        self.assertEqual(scale(_product(sum_a, middle, sum_a),
                               CO / sectors ** 3), defect)

        omega = _diagonal([CO * (-1) ** ((i >> (1 + (i & 1))) & 1)
                           for i in range(logical_size)])
        ideal_gram = scale(add(eye(logical_size), omega), CO / 4)
        ideal_defect = _product(
            va, _subtract(ideal_gram, _conjugate(logical_l, ideal_gram)), adjoint(va))
        address_zero = _diagonal([CO * (1 - (i & 1))
                                  for i in range(logical_size)])
        self.assertEqual(multiply(ideal_defect, ideal_defect),
                         scale(_product(va, address_zero, adjoint(va)), CO / 4))
        self.assertEqual(sum(entry.abs2() for row in ideal_defect for entry in row)
                         / size, Fraction(1, 32))

    def test_paired_lookup_defect_and_full_physical_trace(self):
        for address_bits in (1, 2):
            rows = 1 << address_bits
            for selected in range(rows):
                squared_sum = Fraction(0)
                columns = 0
                for address in range(rows):
                    for bank_c in range(1 << rows):
                        for bank_d in range(1 << rows):
                            omega = (-1) ** (((bank_c ^ bank_d) >> address) & 1)
                            flipped = (-1) ** (((bank_c ^ bank_d ^ (1 << selected))
                                               >> address) & 1)
                            defect = Fraction(omega - flipped, 4)
                            expected = Fraction(omega, 2) if address == selected else 0
                            self.assertEqual(defect, expected)
                            squared_sum += defect * defect
                            columns += 1
                logical_hs_squared = squared_sum / columns
                self.assertEqual(logical_hs_squared, Fraction(1, 4 * rows))
                for clean in (1, 2, 3):
                    self.assertEqual(logical_hs_squared / (1 << clean),
                                     Fraction(1, 1 << (address_bits + clean + 2)))


if __name__ == '__main__':
    unittest.main()
