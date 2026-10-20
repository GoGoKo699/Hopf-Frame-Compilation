"""Literal four-band representation guards, not depth lower bounds.

The identity and nontrivial Boolean fixtures protect genuine original
physical-pair provenance and phases. Both have shorter complete native
realizations, so their displayed pair structure implies no depth floor.
"""
import unittest

from tests.exact_matrices import CO, CZ, adjoint, eye, multiply, native_columns, scale
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_dirty_indicator_depth import _parity_basis
from tests.test_original_middle_pair_insertions import _ghz, _pair_coordinates
from tests.test_quadratic_native_query import _odd_parity_band
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _conjugate, _symplectic, _tau


class FourBandOriginalPairScoutTests(unittest.TestCase):
    def test_proper_original_pair_in_a_literal_four_band_identity(self):
        width = 3
        k_word = _ghz(width)
        ea_word = [('T', bit) for bit in range(width)]
        coordinates = _pair_coordinates(2)
        eb_word = _inverse(coordinates) + [('T', 0), ('T', 1)] + coordinates
        ed_word = _inverse(eb_word)
        word = (k_word + ea_word + eb_word + ed_word
                + _inverse(ea_word) + _inverse(k_word))
        k, ea, eb, ed, encoder = [native_columns(width, source) for source in
                                  (k_word, ea_word, eb_word, ed_word, coordinates)]

        clean, pa = _tau(8, width), _tau(7, width)
        b, d = _tau(3, width), _tau(27, width)
        gamma = scale(_tau(24, width), -1)
        self.assertEqual(_conjugate(k, clean), pa)
        self.assertEqual(_conjugate(encoder, _tau(8, width)), b)
        self.assertEqual(_conjugate(encoder, _tau(16, width)), d)
        self.assertEqual(multiply(b, d), gamma)
        j = _conjugate(ea, pa)
        self.assertEqual(_conjugate(j, b), d)
        self.assertEqual(_conjugate(j, d), b)
        self.assertEqual(_conjugate(eb, j), j)
        self.assertEqual(_conjugate(ed, j), j)
        for band in (ea, eb, ed):
            self.assertEqual(_conjugate(band, gamma), gamma)

        active, pair = [8, 16, 32], [3, 27]
        self.assertEqual(_rank(active), 3)
        self.assertEqual(_rank(pair), 2)
        self.assertEqual(pair[0] ^ pair[1], active[0] ^ active[1])
        self.assertNotEqual(pair[0] ^ pair[1], active[0] ^ active[1] ^ active[2])
        self.assertEqual([sum(_symplectic(axis, a, width) << bit
                              for bit, a in enumerate(active)) for axis in pair],
                         [3, 3])

        actual = native_columns(width, word)
        # Equality checks all eight physical columns, including both clean
        # sectors. It also gives every arbitrary dirty/reference extension.
        self.assertEqual(actual, eye(8))
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(ed, adjoint(eb))
        for source in (word, _inverse(word)):
            self.assertEqual([len(band) for band in _native_t_bands(source)],
                             [3, 2, 2, 3])
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in source), 10)

    def test_disjoint_actual_pairs_complete_a_nontrivial_boolean_source(self):
        # Physical order x1,x2,z,d,w; x1 is the initialized clean qubit.
        width, size = 5, 32
        x_second = [('H', 1), ('S', 1), ('S', 1), ('H', 1)]
        endpoint_coordinates = _pair_coordinates(2) + [('S', 0)] + x_second
        ea_word = (_inverse(endpoint_coordinates)
                   + [('T', 0), ('TDG', 1), ('T', 4)]
                   + endpoint_coordinates)
        first = [('CX', 3, wire) for wire in (0, 1, 2)]
        second = ([('CX', 2, 3)]
                  + [('CX', 3, wire) for wire in (0, 1)]
                  + [('CX', wire, 2) for wire in (0, 1)]
                  + [('CX', 3, 2)]
                  + [('CX', wire, 3) for wire in (0, 1)])
        eb_word = _odd_parity_band(first, [(0, 1), (1, 1), (2, 1), (3, -1)])
        ed_word = _odd_parity_band(second, [(0, -1), (1, -1), (2, -1), (3, 1)])
        k_word = [('CX', 4, 0), ('H', 4)]
        word = (k_word + ea_word + eb_word + ed_word
                + _inverse(ea_word) + _inverse(k_word))
        reduced = k_word + eb_word + ed_word + _inverse(k_word)
        encoder, ea, eb, ed, k = [native_columns(width, source) for source in
                                (endpoint_coordinates, ea_word, eb_word,
                                 ed_word, k_word)]
        active = [35, 67, 512]  # Y1 X2, X1 Y2, Z_w, with literal positive signs.
        for bit, axis in zip((0, 1, 4), active):
            self.assertEqual(_conjugate(encoder, _tau(1 << (width + bit), width)),
                             _tau(axis, width))
        pa = _tau(48, width)  # Z1 X_w.
        self.assertEqual(_conjugate(k, _tau(32, width)), pa)
        self.assertEqual(_rank(active), 3)
        self.assertEqual([_symplectic(axis, 48, width) for axis in active],
                         [1, 1, 1])

        b_axes = [9 << width, 10 << width, 12 << width, 8 << width]
        d_axes = [13 << width, 14 << width, 11 << width, 15 << width]
        self.assertFalse(set(b_axes) & set(d_axes))
        self.assertEqual((_rank(b_axes), _rank(d_axes)), (4, 4))
        j = _conjugate(ea, pa)
        for axes in (b_axes, d_axes):
            for before, after in zip(axes, [axes[1], axes[0], *axes[2:]]):
                self.assertEqual(_conjugate(j, _tau(before, width)),
                                 _tau(after, width))
            self.assertEqual(multiply(_tau(axes[0], width), _tau(axes[1], width)),
                             _tau(96, width))
        self.assertEqual(b_axes[0] ^ b_axes[1], active[0] ^ active[1])
        self.assertNotEqual(b_axes[0] ^ b_axes[1],
                            active[0] ^ active[1] ^ active[2])
        self.assertEqual([sum(_symplectic(axis, a, width) << bit
                              for bit, a in enumerate(active)) for axis in d_axes],
                         [3, 3, 0, 0])
        for band in (ea, eb, ed):
            self.assertEqual(_conjugate(band, _tau(96, width)), _tau(96, width))
        # There are no negative J-singletons: W=J in the source normal form.
        self.assertEqual(_conjugate(eb, j), j)
        self.assertEqual(_conjugate(ed, j), j)

        middle = multiply(ed, eb)
        expected_middle = [[CZ for _ in range(size)] for _ in range(size)]
        expected = [[CZ for _ in range(size)] for _ in range(size)]
        for column in range(size):
            x1, x2, z, w = [(column >> bit) & 1 for bit in (0, 1, 2, 4)]
            expected_middle[column][column] = -CO if x1 * x2 * z else CO
            expected[column][column] = -CO if (x1 ^ w) * x2 * z else CO
        self.assertEqual(middle, expected_middle)
        self.assertEqual(multiply(ea, middle), multiply(middle, ea))
        actual = native_columns(width, word)
        # All 32 columns include both clean and dirty values, with literal
        # phase. This is a full operator equality and preserves references.
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(native_columns(width, reduced), actual)
        self.assertEqual(native_columns(width, _inverse(reduced)), adjoint(actual))
        for source in (word, _inverse(word)):
            self.assertEqual([len(band) for band in _native_t_bands(source)],
                             [3, 4, 4, 3])
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in source), 14)
        for source in (reduced, _inverse(reduced)):
            self.assertEqual([len(band) for band in _native_t_bands(source)], [4, 4])
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in source), 8)

    def test_signed_disjoint_pairs_give_a_controlled_swap_middle(self):
        width, size = 5, 32
        q_word = [('CX', 0, 1), ('H', 0)]
        rows = ([9, 11, 12, 8], [13, 15, 10, 14])
        exponents = ([1, -1, 1, -1], [-1, 1, 1, -1])
        coordinates = [_parity_basis(list(range(4)), basis) for basis in rows]
        diagonal = [_odd_parity_band(basis, list(enumerate(powers)))
                    for basis, powers in zip(coordinates, exponents)]
        b_word, d_word = [q_word + band + _inverse(q_word) for band in diagonal]
        x_second = [('H', 1), ('S', 1), ('S', 1), ('H', 1)]
        endpoint_coordinates = _pair_coordinates(2) + [('S', 0)] + x_second
        ea_word = (_inverse(endpoint_coordinates)
                   + [('T', 0), ('TDG', 1), ('T', 4)] + endpoint_coordinates)
        k_word = [('CX', 4, 0), ('H', 4)]
        word = (k_word + ea_word + b_word + d_word
                + _inverse(ea_word) + _inverse(k_word))
        # The diagonal lists still multiply to the original dirty CCZ.
        # X2 changes its control to 1+x2, giving another literal two-band word.
        reduced = k_word + x_second + diagonal[0] + diagonal[1] + x_second + _inverse(k_word)
        ea, eb, ed = [native_columns(width, band) for band in (ea_word, b_word, d_word)]
        j = _conjugate(ea, _tau(48, width))
        labels = ([259, 355, 384, 256], [387, 483, 352, 480])
        self.assertFalse(set(labels[0]) & set(labels[1]))
        for basis, axes in zip(coordinates, labels):
            encoder = native_columns(width, _inverse(basis) + _inverse(q_word))
            signed = [scale(_tau(label, width), sign)
                      for label, sign in zip(axes, [1, -1, 1, 1])]
            self.assertEqual(_rank(axes), 4)
            for bit, axis in enumerate(signed):
                self.assertEqual(_conjugate(encoder, _tau(1 << (width + bit), width)), axis)
            self.assertEqual(_conjugate(j, signed[0]), scale(signed[1], -1))
            self.assertEqual(_conjugate(j, signed[1]), scale(signed[0], -1))
            self.assertEqual(scale(multiply(signed[0], signed[1]), -1),
                             scale(_tau(96, width), -1))
            for axis in signed[2:]:
                self.assertEqual(_conjugate(j, axis), axis)
        for band in (ea, eb, ed):
            self.assertEqual(_conjugate(band, _tau(96, width)), _tau(96, width))
        middle = multiply(ed, eb)
        expected_middle = [[CZ for _ in range(size)] for _ in range(size)]
        expected = [[CZ for _ in range(size)] for _ in range(size)]
        for column in range(size):
            x1, x2, z, w = [(column >> bit) & 1 for bit in (0, 1, 2, 4)]
            row = column ^ (3 if z and x1 != x2 else 0)
            expected_middle[row][column] = CO
            expected[column][column] = -CO if (x1 ^ w) * (1 ^ x2) * z else CO
        self.assertEqual(middle, expected_middle)
        self.assertNotEqual(multiply(ea, middle), multiply(middle, ea))
        actual = native_columns(width, word)
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(native_columns(width, reduced), actual)
        self.assertEqual(native_columns(width, _inverse(reduced)), adjoint(actual))
        for source, expected_sizes in ((word, [3, 4, 4, 3]),
                                       (_inverse(word), [3, 4, 4, 3]),
                                       (reduced, [4, 4]), (_inverse(reduced), [4, 4])):
            sizes = [len(band) for band in _native_t_bands(source)]
            self.assertEqual(sizes, expected_sizes)
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in source), sum(expected_sizes))

    def test_unmatched_native_axis_has_nonclifford_probe_and_boolean_return(self):
        width, size = 5, 32
        q_word = [('CX', 0, 1), ('H', 0)]
        first = [('CX', 3, 0), ('CX', 0, 1), ('CX', 0, 2)]
        second = ([('CX', wire, 0) for wire in (1, 2, 3)]
                  + [('CX', 1, 2), ('CX', 2, 1), ('CX', 1, 2)]
                  + [('CX', 3, 1), ('CX', 3, 2),
                     ('CX', 1, 3), ('CX', 2, 3)])
        rows = ([9, 11, 13, 8], [15, 12, 10, 14])
        coordinates = (first, second)
        exponents = ([1, -1, -1, -1], [1, 1, 1, -1])
        diagonal = [_odd_parity_band(basis, list(enumerate(powers)))
                    for basis, powers in zip(coordinates, exponents)]
        b_word, d_word = [q_word + band + _inverse(q_word) for band in diagonal]
        # Cancel the adjacent Q-inverse/Q pair literally in the full word.
        middle_word = q_word + diagonal[0] + diagonal[1] + _inverse(q_word)
        x_second = [('H', 1), ('S', 1), ('S', 1), ('H', 1)]
        endpoint_coordinates = _pair_coordinates(2) + [('S', 0)] + x_second
        ea_word = (_inverse(endpoint_coordinates)
                   + [('T', 0), ('TDG', 1), ('T', 4)] + endpoint_coordinates)
        k_word = [('CX', 4, 0), ('H', 4)]
        word = (k_word + ea_word + middle_word
                + _inverse(ea_word) + _inverse(k_word))
        reduced = (k_word + x_second + diagonal[0] + diagonal[1]
                   + x_second + _inverse(k_word))
        ea, eb, ed = [native_columns(width, band) for band in (ea_word, b_word, d_word)]
        labels = ([259, 355, 387, 256], [483, 384, 352, 480])
        signs = ([1, -1, 1, 1], [-1, 1, 1, 1])
        signed = [[scale(_tau(label, width), sign)
                   for label, sign in zip(axes, axis_signs)]
                  for axes, axis_signs in zip(labels, signs)]
        self.assertFalse(set(labels[0]) & set(labels[1]))
        self.assertEqual([len(basis) for basis in coordinates], [3, 10])
        for basis, masks, axes, actual_axes in zip(coordinates, rows, labels, signed):
            self.assertEqual(_rank(axes), 4)
            parity_encoder = native_columns(width, _inverse(basis))
            encoder = native_columns(width, _inverse(basis) + _inverse(q_word))
            for bit, (mask, axis) in enumerate(zip(masks, actual_axes)):
                z_axis = _tau(1 << (width + bit), width)
                self.assertEqual(_conjugate(parity_encoder, z_axis),
                                 _tau(mask << width, width))
                self.assertEqual(_conjugate(encoder, z_axis), axis)

        j = _conjugate(ea, _tau(48, width))
        self.assertEqual(_conjugate(j, signed[0][0]), scale(signed[0][1], -1))
        self.assertEqual(_conjugate(j, signed[0][1]), scale(signed[0][0], -1))
        self.assertEqual(scale(multiply(signed[0][0], signed[0][1]), -1),
                         scale(_tau(96, width), -1))
        self.assertEqual(_conjugate(j, signed[0][2]), scale(signed[1][0], -1))
        self.assertEqual(_conjugate(j, signed[1][0]), scale(signed[0][2], -1))
        for axis in [signed[0][3], *signed[1][1:]]:
            self.assertEqual(_conjugate(j, axis), axis)
        for band in (ea, eb, ed):
            self.assertEqual(_conjugate(band, _tau(96, width)), _tau(96, width))
        probe = _conjugate(eb, j)
        self.assertEqual(_conjugate(ed, probe), j)
        # This physical Pauli acquires two entries in every column, so the
        # intermediate reflection is genuinely non-Clifford, not just unequal J.
        image = _conjugate(probe, _tau(34, width))  # Z1 X2.
        self.assertEqual([sum(row[column] != CZ for row in image)
                          for column in range(size)], [2] * size)

        expected_middle = [[CZ for _ in range(size)] for _ in range(size)]
        expected = [[CZ for _ in range(size)] for _ in range(size)]
        for column in range(size):
            x1, x2, z, w = [(column >> bit) & 1 for bit in (0, 1, 2, 4)]
            row = column ^ (3 if z and x1 != x2 else 0)
            expected_middle[row][column] = CO
            expected[column][column] = -CO if (x1 ^ w) * (1 ^ x2) * z else CO
        self.assertEqual(multiply(ed, eb), expected_middle)
        self.assertEqual(native_columns(width, middle_word), expected_middle)
        self.assertEqual(sum(gate[0] not in ('T', 'TDG') for gate in middle_word), 30)
        actual = native_columns(width, word)
        # All 32 physical columns retain both clean blocks and both dirty
        # values. Literal full-operator equality also covers arbitrary references.
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(native_columns(width, reduced), actual)
        self.assertEqual(native_columns(width, _inverse(reduced)), adjoint(actual))
        for source, sizes in ((word, [3, 4, 4, 3]),
                              (_inverse(word), [3, 4, 4, 3]),
                              (reduced, [4, 4]), (_inverse(reduced), [4, 4])):
            self.assertEqual([len(band) for band in _native_t_bands(source)], sizes)
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in source), sum(sizes))


if __name__ == '__main__':
    unittest.main()
