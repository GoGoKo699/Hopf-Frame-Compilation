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


def _polynomial_product(left, right):
    result = set()
    for a in left:
        for b in right:
            monomial = a | b  # Boolean squares reduce before collecting.
            if monomial in result:
                result.remove(monomial)
            else:
                result.add(monomial)
    return result


def _polynomial_substitute(polynomial, forms):
    result = set()
    for monomial in polynomial:
        term = {0}
        for wire, form in enumerate(forms):
            if monomial >> wire & 1:
                term = _polynomial_product(term, form)
        result.symmetric_difference_update(term)
    return result


def _lookup_polynomial(p):
    return {(1 << (p + subset)) | subset for subset in range(1 << p)}


def _polynomial_degree(polynomial):
    return max((monomial.bit_count() for monomial in polynomial), default=-1)


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

    def test_top_symbol_alone_misses_the_two_address_exception(self):
        # The address/top-coefficient swap preserves the highest exterior
        # form for every p. The next derivative layer detects it only when
        # that layer survives a quadratic defect (p at least three).
        for p in (2, 3, 4):
            width = p + (1 << p) + 1
            top_wire = p + (1 << p) - 1
            forms = [{1 << wire} for wire in range(width)]
            forms[0], forms[top_wire] = forms[top_wire], forms[0]
            polynomial = _lookup_polynomial(p)
            top_symbol = {((1 << p) - 1) | (1 << top_wire)}
            self.assertEqual(_polynomial_substitute(top_symbol, forms), top_symbol)
            defect = polynomial ^ _polynomial_substitute(polynomial, forms)
            self.assertEqual(_polynomial_degree(defect), p)
            self.assertEqual(_polynomial_degree(defect) <= 2, p == 2)
            address_mask = (1 << p) - 1
            derivative_symbols = [{address_mask ^ (1 << i)} for i in range(p)]
            # These are exact derivatives by c_[p]\{i}, before substitution.
            for i, symbol in enumerate(derivative_symbols):
                direction_wire = p + (address_mask ^ (1 << i))
                shifted = [{1 << wire} for wire in range(width)]
                shifted[direction_wire].add(0)
                self.assertEqual(polynomial ^ _polynomial_substitute(polynomial, shifted),
                                 symbol)
            essential_span = 0
            moved_span = 0
            for symbol in derivative_symbols:
                essential_span |= next(iter(symbol))
                moved = _polynomial_substitute(symbol, forms)
                self.assertEqual(len(moved), 1)
                moved_span |= next(iter(moved))
            self.assertEqual(essential_span, address_mask)
            self.assertEqual(moved_span, (address_mask ^ 1) | (1 << top_wire))
            self.assertNotEqual(moved_span, address_mask)

    def test_affine_address_relabeling_and_soft_translations_remain_allowed(self):
        for p in (3, 4):
            width = p + (1 << p) + 1
            # Address involution x1 -> x1+x0+1, with its induced dual
            # coefficient transformation. This is an actual affine map,
            # including the Boolean reductions x0*x0=x0.
            address_forms = [{1 << i} for i in range(p)]
            address_forms[1] = {1 << 1, 1 << 0, 0}
            forms = [{1 << wire} for wire in range(width)]
            forms[:p] = address_forms
            for subset in range(1 << p):
                forms[p + subset] = set()
            for subset in range(1 << p):
                transformed = _polynomial_substitute({subset}, address_forms)
                for target in transformed:
                    forms[p + target].add(1 << (p + subset))
            # Certify actual invertibility of the complete affine map.
            for wire in range(width):
                self.assertEqual(_polynomial_substitute(forms[wire], forms), {1 << wire})
            polynomial = _lookup_polynomial(p)
            self.assertEqual(_polynomial_substitute(polynomial, forms), polynomial)
            # Its hard outputs have no soft translation dependence.
            hard = ((1 << p) - 1) | sum(1 << (p + subset)
                                       for subset in range(1 << p)
                                       if subset.bit_count() >= 2)
            for wire, form in enumerate(forms):
                if hard >> wire & 1:
                    self.assertTrue(all(not monomial & ~hard for monomial in form))
            for wire in range(p):
                self.assertTrue(all(not monomial & ~((1 << p) - 1)
                                    for monomial in forms[wire]))
            # Moving a soft translation into the address violates the
            # characteristic hard algebra, even with an arbitrary helper.
            leaked = [{1 << wire} for wire in range(width)]
            leaked[0].add(1 << p)  # c_empty -> address x0
            defect = polynomial ^ _polynomial_substitute(polynomial, leaked)
            self.assertGreater(_polynomial_degree(defect), 2)
            # A constant pair-coefficient shift, in contrast, has exactly
            # a quadratic defect and is an allowed Clifford symmetry.
            pair_shift = [{1 << wire} for wire in range(width)]
            pair_shift[p + 3].add(0)
            self.assertEqual(polynomial ^ _polynomial_substitute(polynomial, pair_shift),
                             {3})

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
