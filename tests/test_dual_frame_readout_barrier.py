"""Exact endpoint and native-word fixtures for the dual-frame readout barrier.

The denominator checks concern separately implemented W and Q. They do
not charge these endpoints after resynthesis of the complete indicator.
Sector simulation retains every dirty Y value and every carrier column.
"""
from fractions import Fraction
from itertools import combinations
from math import comb
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, native_columns, scale,
)
from tests.test_majorana_walsh_query import (
    _controlled_pauli, _inverse, _majoranas, _mask_label, _pauli_word,
    _shift, _walsh_schedule, _word,
)
from tests.test_two_layer_galois_return import _conjugate, _product, _tau


def _determinant(rows):
    """Integer Bareiss elimination, independent of the subset factorization."""
    if not rows:
        return 1
    matrix = [list(row) for row in rows]
    previous, sign = 1, 1
    for column in range(len(matrix) - 1):
        pivot_row = next((row for row in range(column, len(matrix))
                          if matrix[row][column]), None)
        if pivot_row is None:
            return 0
        if pivot_row != column:
            matrix[column], matrix[pivot_row] = matrix[pivot_row], matrix[column]
            sign = -sign
        pivot = matrix[column][column]
        for row in range(column + 1, len(matrix)):
            for target in range(column + 1, len(matrix)):
                numerator = (matrix[row][target] * pivot
                             - matrix[row][column] * matrix[column][target])
                quotient, remainder = divmod(numerator, previous)
                assert remainder == 0
                matrix[row][target] = quotient
            matrix[row][column] = 0
        previous = pivot
    return sign * matrix[-1][-1]


def _minor(rows, columns):
    return _determinant([[(-1) ** ((row & column).bit_count())
                          for column in columns] for row in rows])


def _gamma_product(indices, width):
    indices = tuple(indices)
    factors = [_tau(_majoranas(width)[index], width) for index in indices]
    product = _product(eye(1 << width), *factors)
    phase = (CO, CI, -CO, -CI)[(len(indices) * (len(indices) - 1) // 2) % 4]
    return scale(product, phase)


def _overlap(left, right, normalization=None):
    product = _product(left, right)
    return sum((product[j][j] for j in range(len(product))), CZ) / (
        len(product) if normalization is None else normalization)


def _counts(word):
    return (sum(gate[0] in ('T', 'TDG') for gate in word),
            sum(gate[0] not in ('T', 'TDG') for gate in word))


class DualFrameReadoutBarrierTests(unittest.TestCase):
    def test_downset_determinants_telescope_and_attain_each_grade_ceiling(self):
        for p in range(2, 7):
            m, radius = 1 << p, (p - 1) // 2
            ordered = sorted(range(m), key=lambda value: (value.bit_count(), value))
            cumulative, ceilings = 0, []
            for k in range(m + 1):
                if k:
                    cumulative += ordered[k - 1].bit_count()
                downset = ordered[:k]
                # Every initial segment in weight order is a downset.
                self.assertTrue(all(subset in downset for value in downset
                                    for subset in range(m) if subset & value == subset))
                determinant = _minor(downset, downset)
                self.assertEqual(abs(determinant), 1 << cumulative)
                ceiling = k * p - 2 * cumulative
                ceilings.append(ceiling)
                self.assertGreaterEqual(ceiling, 0)
                # Clearing exactly this many sqrt(2) factors leaves an odd
                # integer numerator for the normalized determinant.
                self.assertEqual(k * p - 2 * (abs(determinant).bit_length() - 1),
                                 ceiling)
            selected = [value for value in range(m) if value.bit_count() <= radius]
            telescoping = sum((p - 2 * j) * comb(p, j) for j in range(radius + 1))
            self.assertEqual(telescoping, p * comb(p - 1, radius))
            self.assertEqual(max(ceilings), telescoping)
            self.assertEqual(ceilings[len(selected)], telescoping)
            self.assertEqual(Fraction(_minor(selected, selected) ** 2,
                                      1 << (p * len(selected))),
                             Fraction(1, 1 << telescoping))

        # Off-chain minors also obey the ceiling in the smallest full case.
        weights = [0, 1, 1, 2]
        for k in range(5):
            divisor = 1 << sum(weights[:k])
            for rows in combinations(range(4), k):
                for columns in combinations(range(4), k):
                    self.assertEqual(_minor(rows, columns) % divisor, 0)

    def test_native_walsh_eight_minor_and_all_dirty_character_coefficients(self):
        p, width, selected = 3, 4, (0, 1, 2, 4)
        word = _word(_walsh_schedule(p))
        native = native_columns(width, word)
        inverse = native_columns(width, _inverse(word))
        self.assertEqual(inverse, adjoint(native))
        self.assertEqual(_counts(word)[0], 12)
        source = _gamma_product(selected, width)
        image = _conjugate(native, source)
        inverse_image = _conjugate(inverse, source)
        self.assertEqual(_overlap(source, image), -CO * Fraction(1, 8))
        self.assertEqual(_overlap(source, inverse_image), -CO * Fraction(1, 8))

        weights = {}
        for subset in combinations(range(8), len(selected)):
            determinant = Fraction(_minor(subset, selected), 64)
            coefficient = _overlap(_gamma_product(subset, width), image)
            self.assertEqual(coefficient, CO * determinant)
            reverse_coefficient = _overlap(_gamma_product(subset, width), inverse_image)
            self.assertEqual(reverse_coefficient, CO * determinant)
            mask = sum(1 << index for index in subset)
            weights[mask] = determinant * determinant
        self.assertEqual(sum(weights.values()), 1)
        # In each of all 256 complete Y sectors C_Y supplies the character
        # (-1)^(Y dot subset). Fourier extraction isolates a single subset.
        sector_coefficients = [sum(weight * (-1) ** ((y & mask).bit_count())
                                   for mask, weight in weights.items())
                               for y in range(256)]
        selected_mask = sum(1 << index for index in selected)
        coefficient = sum((-1) ** ((y & selected_mask).bit_count()) * value
                          for y, value in enumerate(sector_coefficients)) / 256
        self.assertEqual(coefficient, Fraction(1, 64))

    def test_six_native_readout_calls_give_complete_indicator_with_scalar_gauge(self):
        p, m, width = 2, 4, 2
        gamma = _majoranas(width)
        all_labels = 0
        for label in gamma:
            all_labels ^= label
        dual = [all_labels ^ label for label in gamma]
        walsh = _word(_walsh_schedule(p))
        walsh_inverse = _inverse(walsh)
        # The leftmost dual is delta_0: chronological controlled factors
        # therefore run in descending order. This order matters in extraction.
        controlled = [gate for j in reversed(range(m))
                      for gate in _controlled_pauli(j, dual[j], width, m)]
        q_native = _shift(walsh_inverse, m) + controlled + _shift(walsh, m)
        self.assertEqual(len(controlled), 18)
        self.assertEqual(_counts(q_native), (8, 78))
        self.assertEqual(_counts(_inverse(q_native)), (8, 78))
        sectors, sector_inverses = [], []
        for y in range(1 << m):
            center = [gate for j in reversed(range(m)) if y >> j & 1
                      for gate in _pauli_word(dual[j], width)]
            sector_word = walsh_inverse + center + walsh
            matrix = native_columns(width, sector_word)
            inverse = native_columns(width, _inverse(sector_word))
            self.assertEqual(inverse, adjoint(matrix))
            sectors.append(matrix)
            sector_inverses.append(inverse)

        masks, mask_inverses = [], []
        for x in range(1 << p):
            word = [gate for bit in range(p) if x >> bit & 1
                    for gate in _pauli_word(_mask_label(p, bit), width)]
            masks.append(native_columns(width, word))
            mask_inverses.append(native_columns(width, _inverse(word)))
        mask_word = [gate for bit in range(p)
                     for gate in _controlled_pauli(bit, _mask_label(p, bit), width, p + m)]
        self.assertEqual(_counts(mask_word), (0, 9))
        x_word = _shift(_pauli_word(1, 1), p)
        hadamards = [('H', p + j) for j in range(m)]
        parity = [('S', p + j) for j in range(m) for _ in range(2)]
        d_events = ('PDG', 'X', 'QDG', 'X', 'Q', 'P')
        d_inverse_events = ('PDG', 'QDG', 'X', 'Q', 'X', 'P')
        events = ('QDG',) + d_inverse_events + ('Q',) + d_events

        for gauge in (0, 1):
            source_word = [('S', 0)] * gauge + q_native
            self.assertEqual(_counts(source_word), (8, 78 + gauge))
            self.assertEqual(_counts(_inverse(source_word)), (8, 78 + 3 * gauge))
            source = _shift(source_word, p)
            source_inverse = _inverse(source)
            d_word = (_inverse(mask_word) + x_word + source_inverse
                      + x_word + source + mask_word)
            indicator = (hadamards + source_inverse + _inverse(d_word)
                         + source + d_word + parity + hadamards)
            for orientation in (indicator, _inverse(indicator)):
                self.assertEqual(_counts(orientation), (48, 536 + 12 * gauge))
                self.assertTrue(all(gate[0] in ('H', 'S', 'CX', 'T', 'TDG')
                                    for gate in orientation))
                self.assertLess(max(wire for gate in orientation for wire in gate[1:]),
                                p + m + width)
            common_g = max(_counts(source_word)[1], _counts(_inverse(source_word))[1])
            self.assertLessEqual(_counts(indicator)[1],
                                 6 * common_g + 20 * p * width + 4 * m + 16)

            for x in range(1 << p):
                phases = []
                for initial_y in range(1 << m):
                    y, actual = initial_y, eye(1 << width)
                    for event in events:
                        if event == 'X':
                            y ^= 1
                            continue
                        if event in ('Q', 'QDG'):
                            phase = CI if gauge and y & 1 else CO
                            operation = (scale(sectors[y], phase) if event == 'Q'
                                         else scale(sector_inverses[y], phase.conj()))
                        else:
                            operation = masks[x] if event == 'P' else mask_inverses[x]
                        actual = _product(operation, actual)
                    self.assertEqual(y, initial_y)
                    self.assertEqual(actual, scale(eye(1 << width),
                                                  (-1) ** (y.bit_count() + ((y >> x) & 1))))
                    corrected = scale(actual, (-1) ** y.bit_count())
                    self.assertEqual(corrected, scale(eye(1 << width), (-1) ** ((y >> x) & 1)))
                    phases.append((-1) ** ((y >> x) & 1))
                # The endpoint Hadamards turn the verified phase into the
                # complete arbitrary-Y XOR, with all four carrier columns fixed.
                for initial_y in range(1 << m):
                    for final_y in range(1 << m):
                        amplitude = Fraction(sum(phases[y] * (-1) ** (
                            ((initial_y ^ final_y) & y).bit_count())
                            for y in range(1 << m)), 1 << m)
                        self.assertEqual(amplitude, int(final_y == initial_y ^ (1 << x)))

    def test_clean_projector_contraction_sums_full_ptm_entries_without_extra_average(self):
        word = [('CX', 0, 1), ('T', 1), ('CX', 0, 1)]
        native = native_columns(2, word)
        self.assertEqual(native_columns(2, _inverse(word)), adjoint(native))
        self.assertEqual(_counts(word), (1, 2))
        single_t = native_columns(1, [('T', 0)])
        for column in range(2):
            for row in range(4):
                self.assertEqual(native[row][column], single_t[row][column] if row < 2 else CZ)
        pauli_in, pauli_out, clean_z = _tau(1, 2), _tau(5, 2), _tau(8, 2)
        projector = scale(add(eye(4), clean_z), Fraction(1, 2))
        contracted = _overlap(pauli_out,
                              _conjugate(native, _product(pauli_in, projector)),
                              normalization=2)
        full_entries = [_overlap(pauli_out, _conjugate(native, term))
                        for term in (pauli_in, _product(pauli_in, clean_z))]
        expected = CQ2(Q2(0, Fraction(1, 2)))
        self.assertEqual(full_entries, [CZ, expected])
        self.assertEqual(contracted, expected)
        self.assertEqual(contracted, sum(full_entries, CZ))
        self.assertNotEqual(contracted, sum(full_entries, CZ) / 2)


if __name__ == '__main__':
    unittest.main()
