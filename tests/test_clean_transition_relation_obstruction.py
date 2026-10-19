"""Exact easy targets with an unbounded single-list relation requirement.

The native words keep the clean transition and its actual inverse.  The
target is a product of disjoint CCZ gates, already realizable with no
clean work in three native bands; no hard-query conclusion is intended.
"""
from collections import Counter
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse


def _transition_words(blocks, with_clean):
    if blocks < 2 or blocks % 2:
        raise ValueError('the cancellation family uses a positive even block count')
    offset = int(with_clean)
    changes = [[], [], []]
    phases = [[], [], []]
    lists = [[], [], []]
    for block in range(blocks):
        t = 1 if block < blocks // 2 else -1
        first, second, third = [offset + 3 * block + j for j in range(3)]
        clean_mask = 1 if with_clean else 0
        local = [1 << bit for bit in (first, second, third)]
        lists[0].extend(clean_mask | axis for axis in local)
        lists[1].extend((clean_mask | local[0] | local[1],
                         clean_mask | local[0] | local[2]))
        lists[2].extend((clean_mask | local[1] | local[2],
                         clean_mask | local[0] | local[1] | local[2]))
        if with_clean:
            changes[0].extend(('CX', 0, bit) for bit in (first, second, third))
            changes[1].extend([('CX', 0, second), ('CX', 0, third)])
        changes[1].extend([('CX', first, second), ('CX', first, third)])
        changes[2].append(('CX', second, first))
        if with_clean:
            changes[2].append(('CX', 0, second))
        changes[2].append(('CX', third, second))
        phases[0].extend(('T' if t == 1 else 'TDG', bit)
                         for bit in (first, second, third))
        phases[1].extend(('TDG' if t == 1 else 'T', bit)
                         for bit in (second, third))
        phases[2].extend([('T' if t == 1 else 'TDG', second),
                          ('TDG' if t == 1 else 'T', first)])
    bands = [_inverse(change) + phase + change
             for change, phase in zip(changes, phases)]
    endpoint = [('H', 0)] if with_clean else []
    word = endpoint + bands[0] + bands[1] + _inverse(bands[2]) + endpoint
    return word, lists


def _sign_bit(value, blocks):
    return sum(((value >> (3 * block)) & 7) == 7
               for block in range(blocks)) % 2


class CleanTransitionRelationObstructionTests(unittest.TestCase):
    def test_native_three_band_words_all_columns_and_actual_inverses(self):
        for with_clean in (False, True):
            width = 6 + int(with_clean)
            word, lists = _transition_words(2, with_clean)
            actual = native_columns(width, word)
            self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
            for column in range(1 << width):
                data = column >> int(with_clean)
                sign = -CO if _sign_bit(data, 2) else CO
                for row in range(1 << width):
                    self.assertEqual(actual[row][column], sign if row == column else CZ)
            self.assertEqual([_rank(axes) for axes in lists], [6, 4, 4])
            if with_clean:
                self.assertTrue(all(axis & 1 for axes in lists for axis in axes))
        for blocks in (2, 4, 8):
            for with_clean in (False, True):
                word, lists = _transition_words(blocks, with_clean)
                expected_clifford = 20 * blocks + 2 if with_clean else 8 * blocks
                for actual_word in (word, _inverse(word)):
                    counts = Counter(gate[0] for gate in actual_word)
                    self.assertEqual(counts['T'] + counts['TDG'], 7 * blocks)
                    self.assertEqual(counts['H'] + counts['CX'], expected_clifford)
                self.assertEqual([len(band) for band in _native_t_bands(word)],
                                 [3 * blocks, 2 * blocks, 2 * blocks])
                self.assertEqual([_rank(axes) for axes in lists],
                                 [3 * blocks, 2 * blocks, 2 * blocks])

    def test_merged_moments_and_every_mixed_translation_direction(self):
        for with_clean, expected_rank in ((False, 6), (True, 7)):
            _, lists = _transition_words(2, with_clean)
            axes = sum(lists, [])
            width = 6 + int(with_clean)
            rows = [sum(((axis >> bit) & 1) << index
                        for index, axis in enumerate(axes)) for bit in range(width)]
            self.assertTrue(all((left & right).bit_count() % 2 == 0
                                for left in rows for right in rows))
            self.assertEqual(_rank(axes), expected_rank)
            self.assertGreaterEqual(len(axes) - _rank(axes), _rank(axes))
        # Mixed directions cannot cancel the quadratic derivatives from
        # separate cubic blocks.  This checks all 63 directions at once.
        truth = [_sign_bit(value, 2) for value in range(64)]
        for direction in range(1, 64):
            coefficients = [truth[value] ^ truth[value ^ direction]
                            for value in range(64)]
            for bit in range(6):
                for mask in range(64):
                    if mask >> bit & 1:
                        coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
            self.assertEqual(max(mask.bit_count() for mask, value in enumerate(coefficients)
                                 if value), 2)


if __name__ == '__main__':
    unittest.main()
