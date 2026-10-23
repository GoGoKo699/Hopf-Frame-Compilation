"""Exact fixtures for shared address motion and its bank-moment limit.

The native fixture keeps every row control and helper arbitrary.  The
moment collisions challenge complete lookup, not just a selected state.
No synchronized scalar filter or complete indicator is supplied here.
"""
from fractions import Fraction
from collections import Counter
import unittest

from tests.exact_matrices import CI, CO, CZ, adjoint, native_columns
from tests.test_affine_orbit_reflection_source import _private_c3x, _private_ccz
from tests.test_half_bilinear_gram_source import _private_toffoli
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts


def _moments(bank, p):
    parity = label_sum = 0
    for label in range(1 << p):
        if (bank >> label) & 1:
            parity ^= 1
            label_sum ^= label
    return parity, label_sum


def _matching(value, pairs):
    for control, target in pairs:
        value ^= ((value >> control) & 1) << target
    return value


def _anchored(value, label, pairs):
    return _matching(value ^ label, pairs) ^ label


def _stage_word():
    # Four row controls, two address bits, arbitrary parity/phase helpers.
    rows, address_control, address_target, parity, phase = range(4), 4, 5, 6, 7
    fan = [('CX', row, parity) for row in rows]
    controlled = _private_toffoli(parity, address_control, address_target, phase)
    # Exactly labels 1 and 3 have address-control coordinate one.
    translate = [('CX', 1, address_target), ('CX', 3, address_target)]
    return fan + controlled + list(reversed(fan)) + controlled + translate


def _parallel_private_toffolis(terms, helpers):
    """Eight native phase layers; outputs are not controls in this wall."""
    targets = sorted({term[2] for term in terms})
    word = [('H', target) for target in targets]
    old = 0
    for vertex in (0, 1, 3, 2, 6, 7, 5, 4):
        difference = old ^ vertex
        if difference:
            bit = difference.bit_length() - 1
            word.extend(('CX', term[bit], helper)
                        for term, helper in zip(terms, helpers))
        word.extend(('TDG' if vertex.bit_count() % 2 else 'T', helper)
                    for helper in helpers)
        old = vertex
    word.extend(('CX', term[2], helper)
                for term, helper in zip(terms, helpers))
    word.extend(('H', target) for target in targets)
    return word


def _selector_stage_word():
    # Rows 0 and 1 of a two-bit matching; all other row products vanish.
    # The two selectors, two bank bits and all three helpers are arbitrary.
    s0, s1, y0, y1, xc, xt, parity, d0, d1 = range(9)
    bilinear = _parallel_private_toffolis(
        ((s0, y0, parity), (s1, y1, parity)), (d0, d1))
    controlled = _private_toffoli(parity, xc, xt, d0)
    translate = _private_toffoli(s1, y1, xt, d0)
    return (bilinear + controlled + _inverse(bilinear)
            + _inverse(controlled) + translate)


def _selector_stage(address, selectors, bank, p, pairs):
    parity, label_sum = _moments(selectors & bank, p)
    if parity:
        address = _matching(address, pairs)
    return address ^ _matching(label_sum, pairs) ^ label_sum


def _rational_product(left, right):
    return [[sum(a * b for a, b in zip(row, column))
             for column in zip(*right)] for row in left]


def _mixed_selector_matrix():
    # Selector is the low bit. Opposite anchored CNOTs act on two address bits.
    g = lambda x: x ^ ((x & 1) << 1)
    h = lambda x: x ^ ((x >> 1) & 1)
    a = [[Fraction((i & 1) == (j & 1)
                   and (i >> 1) == (g(j >> 1) if j & 1 else j >> 1))
          for j in range(8)] for i in range(8)]
    b = [[Fraction(int((i >> 1) == (j >> 1))
                   + (1 if (i & 1) == (j & 1) else -1)
                   * int((i >> 1) == h(j >> 1)), 2)
          for j in range(8)] for i in range(8)]
    return _rational_product(_rational_product(_rational_product(a, b), a), b)


def _private_cs(a, b, dirty, phase):
    second_difference = [
        ('T', dirty), ('CX', a, dirty), ('TDG', dirty), ('CX', b, dirty),
        ('T', dirty), ('CX', a, dirty), ('TDG', dirty), ('CX', b, dirty)]
    return _inverse(second_difference) + _private_ccz(dirty, a, b, phase)


def _private_ccs(a, b, c, dirty, auxiliary, phase):
    toggle = _private_toffoli(a, b, dirty, phase)
    sign = _private_cs(dirty, c, auxiliary, phase)
    return (_inverse(sign) + toggle + sign + toggle
            + _private_c3x(a, b, c, dirty, auxiliary, phase)[1:-1])


class SynchronizedAddressStageTests(unittest.TestCase):
    def test_native_stage_returns_all_helpers_with_literal_phase(self):
        word = _stage_word()
        self.assertEqual(_counts(word), (16, 30))
        self.assertEqual(_counts(_inverse(word)), (16, 30))
        actual = native_columns(8, word)
        for column in range(256):
            bank, address = column & 15, (column >> 4) & 3
            parity, label_sum = _moments(bank, 2)
            output = (_matching(address, ((0, 1),)) if parity else address)
            output ^= _matching(label_sum, ((0, 1),)) ^ label_sum
            target = (column & ~(3 << 4)) | (output << 4)
            for row in range(256):
                self.assertEqual(actual[row][column], CO if row == target else CZ)
        self.assertEqual(native_columns(8, _inverse(word)), adjoint(actual))

    def test_anchored_product_matches_shared_moments(self):
        fixtures = (
            (2, ((0, 1),), range(16)),
            (3, ((2, 0),), range(256)),
            (4, ((0, 1), (2, 3)), (0, 1, 3, 0x55, 0x1234, 0x8001, 0xffff)),
        )
        for p, pairs, banks in fixtures:
            for bank in banks:
                parity, label_sum = _moments(bank, p)
                for address in range(1 << p):
                    expected = _matching(address, pairs) if parity else address
                    expected ^= _matching(label_sum, pairs) ^ label_sum
                    for labels in (range(1 << p), reversed(range(1 << p))):
                        actual = address
                        for label in labels:
                            if (bank >> label) & 1:
                                actual = _anchored(actual, label, pairs)
                        self.assertEqual(actual, expected)

    def test_four_corner_collision_flips_required_lookup_sector(self):
        for p in (2, 3):
            corners = (0, 1, 2, 3)
            delta = sum(1 << label for label in corners)
            self.assertEqual(_moments(delta, p), (0, 0))
            for bank in range(1 << (1 << p)):
                other = bank ^ delta
                self.assertEqual(_moments(bank, p), _moments(other, p))
                for address in range(1 << p):
                    left = 1 - 2 * ((bank >> address) & 1)
                    right = 1 - 2 * ((other >> address) & 1)
                    self.assertEqual(left == right, address not in corners)
                    if address in corners:
                        self.assertEqual(abs(left - right), 2)
                        gram_left, gram_right = Fraction(1 + left, 4), Fraction(1 + right, 4)
                        self.assertEqual(abs(gram_left - gram_right), Fraction(1, 2))

    def test_different_generator_stages_cannot_be_reordered_by_row(self):
        first, second = ((0, 1),), ((1, 0),)
        labels, address = (0, 1), 0
        by_row = address
        for label in labels:
            by_row = _anchored(_anchored(by_row, label, first), label, second)
        by_stage = address
        for pairs in (first, second):
            for label in labels:
                by_stage = _anchored(by_stage, label, pairs)
        self.assertNotEqual(by_row, by_stage)

    def test_native_independent_selector_stage_returns_every_helper(self):
        word = _selector_stage_word()
        self.assertEqual(_counts(word), (56, 66))
        actual = native_columns(9, word)
        for column in range(512):
            selectors, bank = column & 3, (column >> 2) & 3
            address = (column >> 4) & 3
            output = _selector_stage(address, selectors, bank, 2, ((0, 1),))
            target = (column & ~(3 << 4)) | (output << 4)
            for row in range(512):
                self.assertEqual(actual[row][column], CO if row == target else CZ)
        self.assertEqual(native_columns(9, _inverse(word)), adjoint(actual))

    def test_dirty_selector_flip_echo_collapses_on_every_small_fiber(self):
        for p, banks, selectors in (
                (2, range(16), range(16)),
                (3, (0, 1, 0x55, 0x96, 0xff), (0, 1, 0x55, 0xaa, 0xff))):
            pairs, mask = ((0, 1),), (1 << (1 << p)) - 1
            for bank in banks:
                for selector in selectors:
                    for address in range(1 << p):
                        # Rightmost selector flip, inverse stage, flip, stage.
                        moved = _selector_stage(address, selector ^ mask, bank,
                                                p, pairs)
                        actual = _selector_stage(moved, selector, bank, p, pairs)
                        expected = _selector_stage(address, mask, bank, p, pairs)
                        self.assertEqual(actual, expected)

    def test_unmixed_selector_acceptance_retains_only_support_span(self):
        for p in (2, 3):
            full = (1 << (1 << p)) - 1
            banks = (full, full ^ 1)
            distributions = []
            for bank in banks:
                counts = Counter(_moments(selector & bank, p)
                                 for selector in range(1 << (1 << p)))
                distributions.append({key: Fraction(count, 1 << (1 << p))
                                      for key, count in counts.items()})
            self.assertEqual(distributions[0], distributions[1])
            self.assertEqual(len(distributions[0]), 1 << (p + 1))
            # A concrete noncommuting two-stage accepted address matrix.
            matrices = []
            for bank in banks:
                columns = []
                for address in range(1 << p):
                    outputs = Counter()
                    for selector in range(1 << (1 << p)):
                        moved = _selector_stage(address, selector, bank,
                                                p, ((0, 1),))
                        moved = _selector_stage(moved, selector, bank,
                                                p, ((1, 0),))
                        outputs[moved] += 1
                    columns.append(tuple(Fraction(outputs[row], 1 << (1 << p))
                                         for row in range(1 << p)))
                matrices.append(tuple(columns))
            self.assertEqual(matrices[0], matrices[1])
            self.assertNotEqual(banks[0] & 1, banks[1] & 1)

    def test_native_mixed_selector_commutator_has_complete_literal_action(self):
        for width, primitive, mask, counts in (
                (4, _private_cs(0, 1, 2, 3), 3, (12, 12)),
                (6, _private_ccs(0, 1, 2, 3, 4, 5), 7, (72, 80))):
            self.assertEqual(_counts(primitive), counts)
            self.assertEqual(_counts(_inverse(primitive)), counts)
            matrix = native_columns(width, primitive)
            self.assertEqual(native_columns(width, _inverse(primitive)), adjoint(matrix))
            # Logical inputs are used only as physical Clifford controls.
            self.assertTrue(all(gate[1] >= mask.bit_length()
                                for gate in primitive if gate[0] in ('T', 'TDG')))
            for column in range(1 << width):
                for row in range(1 << width):
                    expected = (CI if column & mask == mask else CO) if row == column else CZ
                    self.assertEqual(matrix[row][column], expected)
        # s,y,x0,x1,d,g: one row bank bit and two helpers, all arbitrary.
        a = _private_c3x(0, 1, 2, 3, 4, 5)
        b = [('H', 0)] + _private_c3x(0, 1, 3, 2, 4, 5) + [('H', 0)]
        word = _inverse(b) + _inverse(a) + b + a
        self.assertEqual(_counts(word), (128, 156))
        self.assertEqual(_counts(_inverse(word)), (128, 156))
        actual = native_columns(6, word)
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))
        logical = _mixed_selector_matrix()
        for column in range(64):
            for row in range(64):
                if (row ^ column) & (2 | 48):
                    expected = CZ
                elif not column & 2:
                    expected = CO if row == column else CZ
                else:
                    i = (row & 1) | ((row >> 2 & 3) << 1)
                    j = (column & 1) | ((column >> 2 & 3) << 1)
                    expected = CO * logical[i][j]
                self.assertEqual(actual[row][column], expected)

    def test_mixed_selector_leakage_power_polynomial_and_short_echoes(self):
        w = _mixed_selector_matrix()
        g = lambda x: x ^ ((x & 1) << 1)
        h = lambda x: x ^ ((x >> 1) & 1)
        leakage = [[Fraction(int(i == g(j)) + int(i == g(h(j)))
                             - int(i == h(g(j))) - int(i == h(g(h(j)))), 4)
                    for j in range(4)] for i in range(4)]
        self.assertEqual([[w[2*i][2*j+1] for j in range(4)] for i in range(4)],
                         leakage)
        gram = _rational_product(list(zip(*leakage)), leakage)
        self.assertEqual(sum(gram[i][i] for i in range(4)), Fraction(3, 4))
        self.assertEqual(_rational_product(gram, gram),
                         [[Fraction(3, 4) * value for value in row] for row in gram])
        w2, identity = _rational_product(w, w), [
            [Fraction(i == j) for j in range(8)] for i in range(8)]
        w3 = _rational_product(w2, w)
        self.assertEqual(sum(w[i][i] for i in range(8)), Fraction(25, 4))
        for i in range(8):
            for j in range(8):
                self.assertEqual(w3[i][j] - Fraction(5, 4) * w2[i][j]
                                 + Fraction(5, 4) * w[i][j] - identity[i][j], 0)
        # Both selector states over the marked address, and over the uniform
        # unmarked vector, are fixed. Normalization is unnecessary here.
        for selector in (0, 1):
            for support in ((0,), (1, 2, 3)):
                vector = [Fraction((i & 1) == selector and (i >> 1) in support)
                          for i in range(8)]
                self.assertEqual([sum(a*b for a, b in zip(row, vector)) for row in w],
                                 vector)
        x = [[Fraction(i == (j ^ 1)) for j in range(8)] for i in range(8)]
        z = [[Fraction((-1 if i & 1 else 1) * (i == j))
              for j in range(8)] for i in range(8)]
        for pauli in (x, z, _rational_product(x, z)):
            echo = _rational_product(_rational_product(_rational_product(w, pauli),
                                                       list(zip(*w))),
                                     list(zip(*pauli)))
            returns = all(echo[2*i][2*j+1] == echo[2*i+1][2*j] == 0
                          and echo[2*i][2*j] == echo[2*i+1][2*j+1]
                          for i in range(4) for j in range(4))
            self.assertFalse(returns)
        # Quarter-dressed opposite directions have only the marked common
        # +1 eigenspace. This one-row projector identity is not a full source.
        for label in range(4):
            projectors = []
            for control, target in ((0, 1), (1, 0)):
                r = [[CZ for _ in range(4)] for _ in range(4)]
                for column in range(4):
                    active = ((column ^ label) >> control) & 1
                    r[column ^ (active << target)][column] = CI if active else CO
                inverse = adjoint(r)
                projector = [[(r[i][j] + inverse[i][j]) * Fraction(1, 2)
                              for j in range(4)] for i in range(4)]
                self.assertEqual(projector, [
                    [CO if i == j and not ((i ^ label) >> control & 1) else CZ
                     for j in range(4)] for i in range(4)])
                projectors.append(projector)
            self.assertEqual(sum(projectors[0][i][i] * projectors[1][i][i]
                                 for i in range(4)), CO)


if __name__ == '__main__':
    unittest.main()
