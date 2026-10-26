"""Exact complete-space receipts for the lookup Clifford symmetry normal form.

The full normal form is analytic; these fixtures keep its phase, arbitrary
helper and genuine address-mixing conclusions separate from depth claims.
"""
import unittest

from tests.exact_matrices import CZ, adjoint, native_columns
from tests.test_joint_lookup_pauli_source import _inverse


def _lookup(value, p):
    x = value & ((1 << p) - 1)
    return sum((value >> (p + subset)) & 1
               for subset in range(1 << p) if subset & x == subset) & 1


def _coefficients(values):
    result = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for subset in range(len(values)):
            if subset >> bit & 1:
                result[subset] ^= result[subset ^ (1 << bit)]
    return result


def _degree(values):
    return max((mask.bit_count() for mask, value in enumerate(_coefficients(values))
                if value), default=-1)


def _classical(value, word):
    for name, *wires in word:
        if name == 'CX':
            value ^= ((value >> wires[0]) & 1) << wires[1]
        elif name == 'X':
            value ^= 1 << wires[0]
        else:
            raise AssertionError(name)
    return value


def _native_classical(word):
    result = []
    for name, *wires in word:
        if name == 'X':
            result += [('H', wires[0]), ('S', wires[0]), ('S', wires[0]), ('H', wires[0])]
        else:
            result.append((name, *wires))
    return result


def _quadratic_phase(values):
    """Literal phase (-1)^g from a certified quadratic truth table."""
    word = []
    for mask, coefficient in enumerate(_coefficients(values)):
        if not coefficient:
            continue
        wires = [wire for wire in range((len(values) - 1).bit_length()) if mask >> wire & 1]
        if len(wires) == 0:
            # Chronological X,Z,X,Z is literally -I on wire zero.
            x = [('H', 0), ('S', 0), ('S', 0), ('H', 0)]
            word += x + [('S', 0), ('S', 0)] + x + [('S', 0), ('S', 0)]
        elif len(wires) == 1:
            word += [('S', wires[0]), ('S', wires[0])]
        elif len(wires) == 2:
            word += [('H', wires[1]), ('CX', wires[0], wires[1]), ('H', wires[1])]
        else:
            raise AssertionError('nonquadratic exponent')
    return word


def _low_word(p):
    word = [('S', p), ('S', p)]
    for i in range(p):
        wire = p + (1 << i)
        word += [('H', wire), ('CX', i, wire), ('H', wire)]
    return word


def _sandwich_lookup(matrix, p):
    phases = [_lookup(value, p) for value in range(len(matrix))]
    return [[-entry if phases[row] ^ phases[column] else entry
             for column, entry in enumerate(entries)] for row, entries in enumerate(matrix)]


class LookupCliffordSymmetriesTests(unittest.TestCase):
    def test_complete_affine_derivative_space_and_small_address_exceptions(self):
        # Includes one fully arbitrary helper. Exhaustive p=2 checks every
        # displacement and all physical input values, not a clean sector.
        for p in (0, 1, 2):
            width = p + (1 << p) + 1
            values = [_lookup(value, p) for value in range(1 << width)]
            hard = ((1 << p) - 1) | sum(1 << (p + subset)
                                                   for subset in range(1 << p)
                                                   if subset.bit_count() >= 2)
            for displacement in range(1 << width):
                affine = _degree([values[value] ^ values[value ^ displacement]
                                  for value in range(1 << width)]) <= 1
                self.assertEqual(affine, p <= 1 or not displacement & hard)
        # At p=3 every address direction has a genuinely higher-order bank
        # term; an arbitrary helper displacement cannot remove it.
        p, width = 3, 12
        values = [_lookup(value, p) for value in range(1 << width)]
        for displacement in (1, 2, 4, 1 | (1 << 11), 1 << 6, 1 << 10,
                             1 << 3, 1 << 4, 1 << 5, 1 << 7, 1 << 11):
            hard = 7 | sum(1 << (p + subset) for subset in range(8)
                           if subset.bit_count() >= 2)
            self.assertEqual(_degree([values[value] ^ values[value ^ displacement]
                                      for value in range(1 << width)]) <= 1,
                             not displacement & hard)

    def test_complete_native_symmetry_with_address_bank_and_helper_mixing(self):
        p, width = 2, 7
        # Wires x1,x2,c_empty,c1,c2,c12,h. The hard swap plus affine
        # hard mixing deliberately moves the address algebra outside itself.
        monomial = [('CX', 0, 5), ('CX', 5, 0), ('CX', 0, 5),
                    ('CX', 0, 1), ('X', 5)]
        native_m = _native_classical(monomial)
        residual = [('H', 3), ('H', 6), ('CX', 1, 6), ('S', 6)]
        low = _low_word(p)
        differences = [_lookup(value, p) ^ _lookup(_classical(value, monomial), p)
                       for value in range(1 << width)]
        self.assertLessEqual(_degree(differences), 2)
        r_word = native_m + residual
        # QSMQ=(QSQ)(QMQ), with QSQ=L S L and QMQ=M D_g.
        conjugate_word = _quadratic_phase(differences) + native_m + low + residual + low
        r = native_columns(width, r_word)
        target = _sandwich_lookup(r, p)
        self.assertEqual(native_columns(width, conjugate_word), target)
        self.assertEqual(native_columns(width, _inverse(r_word)), adjoint(r))
        self.assertEqual(native_columns(width, _inverse(conjugate_word)), adjoint(target))
        # Every column is nonmonomial but stays within its fixed hard values.
        hard_wires = (0, 1, 5)
        for column in range(1 << width):
            rows = [row for row in range(1 << width) if r[row][column] != CZ]
            self.assertEqual(len(rows), 4)
            expected = _classical(column, monomial)
            for row in rows:
                self.assertTrue(all(((row ^ expected) >> wire) & 1 == 0
                                    for wire in hard_wires))
        # Address output x1 really reads the old high-degree coefficient.
        self.assertEqual((_classical(1 << 5, monomial) >> 0) & 1, 1)
        self.assertEqual((_classical(1 << 0, monomial) >> 0) & 1, 0)

    def test_hard_hadamards_have_a_nonpauli_conjugated_generator(self):
        p, width = 2, 7
        values = [_lookup(value, p) for value in range(1 << width)]
        for wire in range(width):
            derivative = [values[value] ^ values[value ^ (1 << wire)]
                          for value in range(1 << width)]
            # (Q H_j Q) Z_j (Q H_j Q) is X_j times this sign.
            # Its phase is Pauli exactly when the derivative is affine.
            self.assertEqual(_degree(derivative) <= 1, wire not in (0, 1, 5))
            if wire not in (0, 1, 5):
                r = native_columns(width, [('H', wire)])
                low = _low_word(p)
                self.assertEqual(_sandwich_lookup(r, p),
                                 native_columns(width, low + [('H', wire)] + low))

    def test_ordinary_bank_affine_masks_match_empty_and_singleton_coefficients(self):
        p, m = 3, 8
        for constant in (0, 1):
            for linear in (0, 1, 3, 7):
                truth = [constant ^ ((x & linear).bit_count() & 1) for x in range(m)]
                coefficients = _coefficients(truth)
                self.assertEqual(coefficients[0], constant)
                for i in range(p):
                    self.assertEqual(coefficients[1 << i], (linear >> i) & 1)
                self.assertTrue(all(not value for subset, value in enumerate(coefficients)
                                    if subset.bit_count() >= 2))


if __name__ == '__main__':
    unittest.main()
