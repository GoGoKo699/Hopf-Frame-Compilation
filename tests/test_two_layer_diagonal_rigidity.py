"""Exact scope checks for native two-layer diagonal rigidity.

The one-clean example is diagonal and non-C3 only on its occupied code.
Its input-dependent fourth-root phase prevents interpreting it as a
literal Boolean XOR query.  All arithmetic is in Q(sqrt(2), i).
"""
import unittest

from tests.exact_matrices import CZ, CO, CI, adjoint, native_columns
from tests.test_dirty_indicator_depth import _ccz_two_layers
from tests.test_t_depth import _word


def _inverse(word):
    names = {'H': 'H', 'CX': 'CX', 'S': 'SDG', 'SDG': 'S',
             'T': 'TDG', 'TDG': 'T'}
    return [(names[gate[0]], *gate[1:]) for gate in reversed(word)]


def _primitive_inverse(word):
    return [expanded for gate in _inverse(word)
            for expanded in ([('S', gate[1])] * 3
                             if gate[0] == 'SDG' else [gate])]


def _reduce_adjacent_s(word):
    """Only reduce literal consecutive S gates on the same wire, modulo four."""
    result = []
    position = 0
    while position < len(word):
        gate = word[position]
        stop = position + 1
        if gate[0] == 'S':
            while stop < len(word) and word[stop] == gate:
                stop += 1
            result.extend([gate] * ((stop - position) % 4))
        else:
            result.append(gate)
        position = stop
    return result


def _counterexample():
    # Little endian: t=0, x_1,...,x_4=1,...,4, and initialized c=5.
    coordinates = ([('CX', 0, j) for j in range(1, 5)]
                   + [('CX', j, 5) for j in range(1, 5)]
                   + [('CX', 0, 5)])
    band = [('T', j) for j in range(6)]
    ez = coordinates + band + list(reversed(coordinates))
    ex = [('H', 0)] + ez + [('H', 0)]
    # K^dagger=A^dagger S^dagger, with A=H S^dagger H.
    final_clifford = [('S', 0)] * 3 + [('H', 0), ('S', 0), ('H', 0)]
    return coordinates, band, ex + _inverse(ez) + final_clifford


def _quadratic(column):
    bits = [(column >> j) & 1 for j in range(1, 5)]
    return (sum(bits) + sum(bits[i] * bits[j]
                           for i in range(4) for j in range(i + 1, 4))) % 2


def _mobius(values, modulus):
    coefficients = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for mask in range(len(values)):
            if mask >> bit & 1:
                coefficients[mask] = (coefficients[mask]
                                      - coefficients[mask ^ (1 << bit)]) % modulus
    return coefficients


class TwoLayerDiagonalRigidityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.coordinates, cls.band, cls.word = _counterexample()
        cls.actual = native_columns(6, cls.word)

    def test_native_coordinates_and_two_physical_band_ledger(self):
        self.assertEqual(len(self.coordinates), 9)
        self.assertEqual([gate[1] for gate in self.band], list(range(6)))
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in self.word), 12)
        self.assertEqual(sum(g[0] == 'CX' for g in self.word), 36)
        self.assertEqual(sum(g[0] == 'H' for g in self.word), 4)
        self.assertEqual(sum(g[0] == 'S' for g in self.word), 4)
        self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in self.word), 44)
        # Locate the two contiguous native bands, not just their total count.
        bands = []
        for gate in self.word:
            if gate[0] in ('T', 'TDG'):
                if not bands or bands[-1] is None:
                    bands.append([])
                bands[-1].append(gate)
            elif bands and bands[-1] is not None:
                bands.append(None)
        bands = [band for band in bands if band is not None]
        self.assertEqual(len(bands), 2)
        for band in bands:
            self.assertEqual(sorted(gate[1] for gate in band), list(range(6)))
        for initial in range(64):
            current = initial
            for _, control, target in self.coordinates:
                current ^= ((current >> control) & 1) << target
            t, c = initial & 1, initial >> 5
            xs = [(initial >> j) & 1 for j in range(1, 5)]
            expected = (t + sum((t ^ x) << j for j, x in enumerate(xs, 1))
                        + ((c ^ t ^ (sum(xs) % 2)) << 5))
            self.assertEqual(current, expected)

    def test_every_occupied_column_literal_phase_and_full_sector_return(self):
        for column in range(32):
            q = _quadratic(column)
            phase = CO if not q else CI * (-1 if column & 1 else 1)
            for row in range(64):
                self.assertEqual(self.actual[row][column],
                                 phase if row == column else CZ)
        for row in range(64):
            for column in range(64):
                if (row >> 5) != (column >> 5):
                    self.assertEqual(self.actual[row][column], CZ)
        # Its unused clean-one sector is not diagonal, so the full-space
        # diagonal theorem has not been contradicted.
        self.assertTrue(any(self.actual[row][column] != CZ
                            for row in range(32, 64)
                            for column in range(32, 64) if row != column))

    def test_literal_inverse_and_explicit_clifford_reduction(self):
        raw = _primitive_inverse(self.word)
        reduced = _reduce_adjacent_s(raw)
        for word, cliffords in ((raw, 52), (reduced, 44)):
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 12)
            self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in word), cliffords)
            self.assertEqual(native_columns(6, word), adjoint(self.actual))

    def test_quartic_phase_is_detected_from_the_exact_occupied_matrix(self):
        phases = (CO, CI, -CO, -CI)
        exponents = [2 * phases.index(self.actual[column][column])
                     for column in range(32)]
        coefficients = _mobius(exponents, 8)
        self.assertEqual(coefficients[0b11110], 4)
        self.assertEqual(max(mask.bit_count() for mask, value
                             in enumerate(coefficients) if value), 4)
        # A fourth coordinate difference of a weighted cubic phase vanishes.
        # Here both target-bit slices have the Boolean sign -1 instead.
        for t in (0, 1):
            fourth_difference = sum((-1) ** x.bit_count()
                                    * exponents[t | (x << 1)] for x in range(16))
            self.assertEqual(fourth_difference % 8, 4)

    def test_output_hadamards_retain_an_unpaid_nonclifford_input_phase(self):
        transformed = native_columns(6, [('H', 0)] + self.word + [('H', 0)])
        for column in range(32):
            q = _quadratic(column)
            output = column ^ q
            phase = CI if q else CO
            for row in range(64):
                self.assertEqual(transformed[row][column],
                                 phase if row == output else CZ)
        self.assertEqual(transformed[3][2], CI)
        self.assertNotEqual(transformed[3][2], CO)
        # q itself is only quadratic. Its i^q lift, required for correction,
        # has a quartic mod-eight coefficient and cannot be a free Clifford.
        q_values = [_quadratic(x << 1) for x in range(16)]
        q_coefficients = _mobius(q_values, 2)
        self.assertEqual(max(mask.bit_count() for mask, value
                             in enumerate(q_coefficients) if value), 2)
        self.assertEqual(_mobius([2 * q for q in q_values], 8)[15], 4)

    def test_zero_clean_two_layer_ccz_with_arbitrary_returned_helper(self):
        schedule = _ccz_two_layers(0, 1, 2, 3)
        self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 2)
        word = _word(schedule)
        actual = native_columns(4, word)
        for column in range(16):
            phase = -CO if column & 7 == 7 else CO
            for row in range(16):
                self.assertEqual(actual[row][column], phase if row == column else CZ)
        self.assertEqual(native_columns(4, _primitive_inverse(word)), adjoint(actual))

    def test_formal_c3_correction_is_real_quartic_but_not_a_supplied_band(self):
        # This multiplies exact phase values only. It supplies NO native
        # implementation of the correction within the original two bands.
        phases = (CO, CI, -CO, -CI)
        real_exponents = []
        for column in range(32):
            xs = [(column >> j) & 1 for j in range(1, 5)]
            pairs = sum(xs[i] * xs[j] for i in range(4)
                        for j in range(i + 1, 4))
            correction = phases[(-sum(xs) + pairs) % 4]
            actual = correction * self.actual[column][column]
            sign = (int(all(xs)) + (column & 1) * _quadratic(column)) % 2
            self.assertEqual(actual, -CO if sign else CO)
            real_exponents.append(4 * sign)
        self.assertEqual(_mobius(real_exponents, 8)[0b11110], 4)

    def test_exceptional_relation_has_nondegenerate_even_quadratic_polar(self):
        # The r-axis relation reduces to d=r-2 independent central bits.
        # Check the exact half-weight formula and polar rank for bounded d.
        for dimension in (2, 4, 6, 8):
            values = [(x.bit_count() + (x.bit_count() % 2)) // 2 % 2
                      for x in range(1 << dimension)]
            coefficients = _mobius(values, 2)
            self.assertTrue(all(not value or mask.bit_count() <= 2
                                for mask, value in enumerate(coefficients)))
            columns = []
            for bit in range(dimension):
                columns.append(sum((values[(1 << bit) ^ (1 << other)]
                                    ^ values[1 << bit] ^ values[1 << other])
                                   << other for other in range(dimension)))
            images = {0}
            for column in columns:
                images |= {value ^ column for value in images}
            self.assertEqual(len(images), 1 << dimension)


if __name__ == '__main__':
    unittest.main()
