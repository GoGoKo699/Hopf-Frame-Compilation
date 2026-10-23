"""Exact fixtures for shared address motion and its bank-moment limit.

The native fixture keeps every row control and helper arbitrary.  The
moment collisions challenge complete lookup, not just a selected state.
No synchronized scalar filter or complete indicator is supplied here.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
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


if __name__ == '__main__':
    unittest.main()
