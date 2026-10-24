"""Exact fixtures for complete coefficient assembly with a serial shear.

The source is complete at every bank size, but its emitted depth is linear
in the address length.  These checks do not imply a sublogarithmic source
or an unrestricted lower bound.
"""
from collections import Counter
import unittest

from tests.exact_matrices import CI, CO, CZ, adjoint, native_columns
from tests.test_joint_lookup_pauli_source import (
    _bands, _inverse, _mobius_bank_word, _normal_form_source,
)


def _coefficient_shear(control, left, right):
    """Product exp(i pi X_left[j] control right[j]/2), two native bands."""
    assert len(left) == len(right) and left
    assert len({control, *left, *right}) == 1 + 2 * len(left)
    hadamards = [('H', target) for target in left]
    pairs = [('CX', a, b) for a, b in zip(left, right)]
    control_word = [('CX', control, target) for target in left]
    first = [gate for a, b in zip(left, right)
             for gate in [('TDG', a), ('T', b)]]
    second = [gate for a, b in zip(left, right)
              for gate in [('T', a), ('TDG', b)]]
    return (hadamards + pairs + first + list(reversed(pairs))
            + control_word + pairs + second + list(reversed(pairs))
            + list(reversed(control_word)) + hadamards)


def _coefficient_source(p):
    """Complete ANF lookup source; arbitrary inputs, no helper wires."""
    assert p >= 2
    mapping = (0, 1, p, p + 3, p + 1, p + 2)
    word = [(gate[0], *(mapping[wire] for wire in gate[1:]))
            for gate in _normal_form_source()]
    for address in range(2, p):
        half = 1 << address
        word += _coefficient_shear(
            address, list(range(p, p + half)),
            list(range(p + half, p + 2 * half)))
    return word


def _raw_phase_source(p):
    wall = _mobius_bank_word(p)
    return wall + _coefficient_source(p) + list(reversed(wall))


def _inverse_coefficient_rows(p, address):
    """Full source inverse on bank-linear rows, ignoring only its phases."""
    rows = [1 << index for index in range(1 << p)]
    for bit in reversed(range(2, p)):
        if (address >> bit) & 1:
            half = 1 << bit
            for index in range(half):
                rows[index] ^= rows[index + half]
    a, b = address & 1, (address >> 1) & 1
    rows[0] ^= (a * rows[1]) ^ (b * rows[2]) ^ (a * b * rows[3])
    return rows


class CompleteCoefficientShearTests(unittest.TestCase):
    def test_two_band_shear_all_small_columns_and_actual_inverse(self):
        for pairs in (1, 2):
            width = 1 + 2 * pairs
            left = list(range(1, 1 + pairs))
            right = list(range(1 + pairs, width))
            word = _coefficient_shear(0, left, right)
            actual = native_columns(width, word)
            phases = (CO, CI, -CO, -CI)
            for column in range(1 << width):
                active = column & 1
                bank = column >> (pairs + 1)
                row = column ^ ((active * bank) << 1)
                phase = phases[(active * bank.bit_count()) % 4]
                for output in range(1 << width):
                    self.assertEqual(actual[output][column],
                                     phase if output == row else CZ)
            self.assertEqual(native_columns(width, _inverse(word)),
                             adjoint(actual))

    def test_native_ledgers_and_distinct_physical_targets(self):
        for pairs in (1, 2, 4, 8):
            word = _coefficient_shear(
                0, list(range(1, pairs + 1)),
                list(range(pairs + 1, 2 * pairs + 1)))
            for native in (word, _inverse(word)):
                bands = _bands(native)
                self.assertEqual(list(map(len, bands)), [2 * pairs] * 2)
                for band in bands:
                    self.assertEqual(len({gate[1] for gate in band}), len(band))
                counts = Counter(gate[0] for gate in native)
                self.assertEqual(counts['T'] + counts['TDG'], 4 * pairs)
                self.assertEqual(counts['H'] + counts['S'] + counts['CX'],
                                 8 * pairs)
        for p in (2, 3, 4, 6):
            size = 1 << p
            for word, cliffords in (
                    (_coefficient_source(p), 8 * size + 36),
                    (_raw_phase_source(p), p * size + 8 * size + 36)):
                for native in (word, _inverse(word)):
                    bands = _bands(native)
                    self.assertEqual(len(bands), 2 * p - 1)
                    for band in bands:
                        self.assertEqual(len({gate[1] for gate in band}), len(band))
                    counts = Counter(gate[0] for gate in native)
                    self.assertEqual(counts['T'] + counts['TDG'], 4 * size - 4)
                    self.assertEqual(counts['H'] + counts['S'] + counts['CX'],
                                     cliffords)
                    self.assertEqual({wire for gate in native for wire in gate[1:]},
                                     set(range(p + size)))

    def test_complete_inverse_row_for_every_address_through_six_bits(self):
        for p in range(2, 7):
            for address in range(1 << p):
                expected = sum(1 << subset for subset in range(1 << p)
                               if subset & address == subset)
                self.assertEqual(_inverse_coefficient_rows(p, address)[0], expected)

    def test_three_address_native_echo_on_independent_and_dense_banks(self):
        p, size = 3, 8
        source = _raw_phase_source(p)
        inverse = _inverse(source)
        # Every address and every independent raw-bank direction, plus dense
        # banks, challenge the former five missing coefficients and return.
        banks = [0, *(1 << bit for bit in range(size)), 0x55, 0xAA, 0xFF]
        inputs = [address | (bank << p)
                  for address in range(size) for bank in banks]
        echo = native_columns(p + size,
                              inverse + [('S', p), ('S', p)] + source,
                              inputs)
        returned = native_columns(p + size, source + inverse, inputs[::12])
        for index, column in enumerate(inputs):
            address, bank = column & (size - 1), column >> p
            sign = -CO if (bank >> address) & 1 else CO
            for row in range(1 << (p + size)):
                self.assertEqual(echo[row][index], sign if row == column else CZ)
        for index, column in enumerate(inputs[::12]):
            for row in range(1 << (p + size)):
                self.assertEqual(returned[row][index], CO if row == column else CZ)


if __name__ == '__main__':
    unittest.main()
