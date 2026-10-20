"""Exact small guards for a fixed original physical-pair completion.

The first fixture retains both clean sectors and actual inverse. The
second protects its non-C3 partner. The third keeps a core-coupled
input axis and both pair-parity sectors; none enumerates completions.
The fourth keeps a longer output quotient relation and pays its shared
clean flip on every sector column, without asserting a Boolean target.
"""
from collections import Counter
import unittest

from tests.exact_matrices import (
    CO, CZ, CI, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_original_middle_pair_insertions import _ghz, _pair_coordinates
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _conjugate, _tau


def _t_phase(pauli):
    omega = native_columns(1, [('T', 0)])[1][1]
    return add(scale(eye(len(pauli)), (CO + omega) / 2),
               scale(pauli, (CO - omega) / 2))


def _product(*factors):
    result = eye(len(factors[0]))
    for factor in factors:
        result = multiply(result, factor)
    return result


class MinimalPhysicalPairCompletionTests(unittest.TestCase):
    def test_original_pair_identity_and_retained_unused_sector(self):
        width = 3
        identity = eye(8)
        k_word = _ghz(2)
        ub_word = ([('CX', 2, 0), ('CX', 2, 1)]
                   + _pair_coordinates(2) + [('H', 2)])
        uc_word = [('CX', 0, 2), ('CX', 1, 2), ('H', 2), ('S', 2), ('S', 2)]
        eb_word = _inverse(ub_word) + [('T', 0), ('T', 1), ('TDG', 2)] + ub_word
        ec_word = _inverse(uc_word) + [('T', 0), ('T', 1), ('T', 2)] + uc_word
        word = (k_word + [('T', 0), ('T', 1)] + eb_word
                + _inverse(ec_word) + _inverse(k_word))

        k = native_columns(width, k_word)
        a = native_columns(width, [('T', 0), ('T', 1)])
        ub = native_columns(width, ub_word)
        uc = native_columns(width, uc_word)
        p = _tau(3, width)
        q = _tau(3 | (3 << width), width)
        r = multiply(p, q)  # Literal R_p=-Z0 Z1, not unsigned Z0 Z1.
        b = _tau(4, width)
        pb, qb, rb = (multiply(axis, b) for axis in (p, q, r))
        clean = _tau(1 << width, width)
        for index, axis in enumerate((pb, qb, b)):
            self.assertEqual(_conjugate(ub, _tau(1 << (width + index), width)), axis)
        for index, axis in enumerate((clean, _tau(2 << width, width), rb)):
            self.assertEqual(_conjugate(uc, _tau(1 << (width + index), width)), axis)
        self.assertEqual(_rank([7, 7 | (3 << width), 4]), 3)
        self.assertEqual(_conjugate(k, clean), p)
        self.assertEqual(r, scale(_tau(3 << width, width), -1))

        j = _conjugate(a, p)
        self.assertEqual(_conjugate(j, pb), qb)
        self.assertEqual(_conjugate(j, qb), pb)
        self.assertEqual(_conjugate(j, b), b)
        eb = native_columns(width, eb_word)
        ec = native_columns(width, ec_word)
        self.assertEqual(_conjugate(eb, j), j)
        self.assertEqual(_conjugate(ec, p), j)

        actual = native_columns(width, word)
        # All eight columns: occupied identity, unused c=a=1 sector i X2.
        for column in range(8):
            expected = [CZ] * 8
            if column & 3 == 3:
                expected[column ^ 4] = CI
            else:
                expected[column] = CO
            self.assertEqual([actual[row][column] for row in range(8)], expected)
        self.assertEqual(_conjugate(actual, clean), clean)
        self.assertEqual(_product(adjoint(actual), actual), identity)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))

        pair_residual = _product(adjoint(k), adjoint(a),
                                 _t_phase(pb), _t_phase(qb), a, k)
        logical_r = _conjugate(adjoint(k), r)
        self.assertEqual(logical_r, scale(_tau(2 << width, width), -1))
        expected_residual = _product(_t_phase(b), _t_phase(multiply(logical_r, b)))
        for column in range(0, 8, 2):
            self.assertEqual([pair_residual[row][column] for row in range(8)],
                             [expected_residual[row][column] for row in range(8)])

        counts = Counter(gate[0] for gate in word)
        self.assertEqual(counts, Counter({'CX': 14, 'H': 12, 'T': 4,
                                          'TDG': 4, 'S': 4, 'SDG': 4}))
        self.assertEqual(Counter(gate[0] for gate in _inverse(word)), counts)
        self.assertEqual([len(band) for band in _native_t_bands(word)], [2, 3, 3])
        self.assertEqual(counts['H'] + counts['CX'] + counts['S'] + 3 * counts['SDG'], 42)

    def test_one_moving_axis_retains_a_non_c3_partner(self):
        # Logical a=0,t=1: H=-Z_a X_t, B=X_t and W=T_Zt.
        width = 2
        coordinates = [('CX', 0, 1), ('H', 1), ('S', 1), ('S', 1)]
        native_coordinates = native_columns(width, coordinates)
        h = scale(_tau(2 | (1 << width), width), -1)
        self.assertEqual(_conjugate(native_coordinates, _tau(2 << width, width)), h)
        partner_word = _inverse(coordinates) + [('T', 1)] + coordinates
        word = [('T', 1)] + partner_word
        actual = native_columns(width, word)
        self.assertEqual(actual, multiply(_t_phase(h), native_columns(width, [('T', 1)])))
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(multiply(adjoint(actual), actual), eye(4))

        # A single exact Pauli witness suffices; no Clifford enumeration.
        # If actual were C3, this first image would be Clifford. Its image
        # of X_a instead has two nonzeros per row, impossible for a Pauli.
        first_image = _conjugate(actual, _tau(2, width))
        second_image = _conjugate(first_image, _tau(1, width))
        self.assertEqual([sum(value != CZ for value in row) for row in second_image],
                         [2, 2, 2, 2])


    def test_core_coupled_input_and_two_rank_two_parity_sectors(self):
        # c=0 is clean; a=1,t=2,s=3 are arbitrary. The new input
        # axis is R_p Z_t, and the extra middle singleton is R_p Z_s.
        # This is a sector-reduction fixture, not a Boolean-degree witness.
        width = 4
        k_word = _ghz(2)
        x_t = [('H', 2), ('S', 2), ('S', 2), ('H', 2)]
        ua_word = [('CX', 0, 2), ('CX', 1, 2)] + x_t
        ub_word = ([('CX', 0, 3), ('CX', 1, 3),
                    ('CX', 2, 0), ('CX', 2, 1)]
                   + _pair_coordinates(2) + [('H', 2)])
        ea_word = _inverse(ua_word) + [('T', 0), ('T', 1), ('T', 2)] + ua_word
        eb_word = (_inverse(ub_word)
                   + [('T', 0), ('T', 1), ('TDG', 2), ('T', 3)] + ub_word)
        word = (k_word + ea_word + eb_word
                + [('TDG', 0), ('TDG', 1)] + _inverse(k_word))

        p_axis = _tau(3, width)
        q_axis = _tau(3 | (3 << width), width)
        rp = multiply(p_axis, q_axis)
        b = _tau(4, width)
        zt, zs = _tau(4 << width, width), _tau(8 << width, width)
        ua = native_columns(width, ua_word)
        ub = native_columns(width, ub_word)
        self.assertEqual(_conjugate(ua, zt), multiply(rp, zt))
        middle_axes = (multiply(p_axis, b), multiply(q_axis, b),
                       b, multiply(rp, zs))
        for index, axis in enumerate(middle_axes):
            self.assertEqual(_conjugate(ub, _tau(1 << (width + index), width)), axis)
        self.assertEqual(_rank([1 << width, 2 << width, 7 << width]), 3)
        self.assertEqual(_rank([7, 7 | (3 << width), 4, 11 << width]), 4)

        actual = native_columns(width, word)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(_product(adjoint(actual), actual), eye(16))
        self.assertEqual(_conjugate(actual, _tau(1 << width, width)),
                         _tau(1 << width, width))
        self.assertEqual([len(band) for band in _native_t_bands(word)], [3, 4, 2])

        k = native_columns(width, k_word)
        a = native_columns(width, [('T', 0), ('T', 1)])
        pair = _product(adjoint(k), adjoint(a), _t_phase(middle_axes[0]),
                        _t_phase(middle_axes[1]), a, k)
        omega = native_columns(1, [('T', 0)])[1][1]
        for central_bit in (0, 1):
            positions = [2 * central_bit + 4 * t + 8 * s
                         for s in (0, 1) for t in (0, 1)]
            block = [[actual[row][column] for column in positions] for row in positions]
            pair_block = [[pair[row][column] for column in positions] for row in positions]
            sign = 2 * central_bit - 1  # R_l=-Z_a.
            expected = _product(_t_phase(scale(_tau(1, 2), sign)),
                                _t_phase(scale(_tau(2 << 2, 2), sign)),
                                _t_phase(scale(_tau(1 << 2, 2), sign)))
            self.assertEqual(block, expected)
            expected_pair = (scale(eye(4), omega) if central_bit == 0
                             else native_columns(2, [('H', 0), ('S', 0), ('H', 0)]))
            self.assertEqual(pair_block, expected_pair)

    def test_one_sided_rank_two_with_four_axis_output_relation(self):
        # Five physical wires, active ranks (2,5), one original pair,
        # and three negative singletons. The C basis changes so that
        # the pair difference uses four actual C axes, not two.
        width = 5
        k_a_word = _ghz(2)
        cz_words = [gate for j in (2, 3, 4)
                    for gate in [('H', j), ('CX', 0, j), ('H', j)]]
        k_c_word = k_a_word + cz_words
        ub_word = _pair_coordinates(2) + [('CX', 0, j) for j in (2, 3, 4)]
        # Orthogonal transvection I+vv^T, v=(1,0,1,1,1), followed
        # by the canonical axes Z0,Z1,Z0Z2,Z0Z3,Z0Z4.
        orthogonal = ([('CX', j, 0) for j in (2, 3, 4)]
                      + [('CX', 0, j) for j in (2, 3, 4)]
                      + [('CX', j, 0) for j in (2, 3, 4)])
        uc_word = orthogonal + [('CX', 0, j) for j in (2, 3, 4)]
        ea_word = [('T', 0), ('T', 1)]
        eb_word = _inverse(ub_word) + [('T', j) for j in range(width)] + ub_word
        ec_word = _inverse(uc_word) + [('T', j) for j in range(width)] + uc_word
        word = k_a_word + ea_word + eb_word + _inverse(ec_word) + _inverse(k_c_word)

        a = native_columns(width, ea_word)
        k_a = native_columns(width, k_a_word)
        k_c = native_columns(width, k_c_word)
        eb, ec = native_columns(width, eb_word), native_columns(width, ec_word)
        ub, uc = native_columns(width, ub_word), native_columns(width, uc_word)
        p_a = _tau(3, width)
        p_c = _tau(3 | (28 << width), width)
        clean = _tau(1 << width, width)
        self.assertEqual(_conjugate(k_a, clean), p_a)
        self.assertEqual(_conjugate(k_c, clean), p_c)
        b_labels = [31, 31 | (3 << width), 5 << width, 9 << width, 17 << width]
        c_labels = [29 << width, 2 << width, 25 << width, 21 << width, 13 << width]
        for index, (b_label, c_label) in enumerate(zip(b_labels, c_labels)):
            z = _tau(1 << (width + index), width)
            self.assertEqual(_conjugate(ub, z), _tau(b_label, width))
            self.assertEqual(_conjugate(uc, z), _tau(c_label, width))
        self.assertEqual(_rank(b_labels), 5)
        self.assertEqual(_rank(c_labels), 5)
        j = _conjugate(a, p_a)
        self.assertEqual(_conjugate(j, _tau(b_labels[0], width)),
                         _tau(b_labels[1], width))
        for label in b_labels[2:]:
            self.assertEqual(_conjugate(j, _tau(label, width)),
                             scale(_tau(label, width), -1))
        self.assertEqual(_conjugate(eb, j), _conjugate(ec, p_c))

        actual = native_columns(width, word)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(_product(adjoint(actual), actual), eye(32))
        self.assertEqual(_conjugate(actual, clean), clean)
        self.assertEqual([len(band) for band in _native_t_bands(word)], [2, 5, 5])

        # R=-Z0Z1 becomes -Z1. Every entry of both physical parity
        # sectors is retained, including both values of the clean wire.
        parity_coordinates = native_columns(width, [('CX', 0, 1)])
        c_parity = _conjugate(adjoint(parity_coordinates), ec)
        j_parity = _conjugate(adjoint(parity_coordinates), j)
        projected_c = [15, 1, 13, 11, 7]
        self.assertEqual(_rank(projected_c), 4)
        self.assertEqual(projected_c[1] ^ projected_c[2]
                         ^ projected_c[3] ^ projected_c[4], 0)
        sector_width = 4
        singleton_axes = [_tau(mask << sector_width, sector_width) for mask in (3, 5, 9)]
        e_minus = _product(*[_t_phase(axis) for axis in singleton_axes])
        k_c_sector_word = ([('H', 0)]
                           + [gate for target in (1, 2, 3)
                              for gate in [('H', target), ('CX', 0, target), ('H', target)]])
        k_c_sector = native_columns(sector_width, k_c_sector_word)
        dressing = native_columns(sector_width, [('H', 1), ('CX', 0, 1), ('H', 1)])
        sector_clean = _tau(1 << sector_width, sector_width)
        h_axis = singleton_axes[0]
        for central_bit in (0, 1):
            positions = [(value & 1) + 2 * central_bit + 4 * (value >> 1)
                         for value in range(16)]
            ec_sector = [[c_parity[row][column] for column in positions] for row in positions]
            j_sector = [[j_parity[row][column] for column in positions] for row in positions]
            k_s_word = [('H', 0)] + ([('S', 0)] if central_bit == 0 else [])
            k_s = native_columns(sector_width, k_s_word)
            self.assertEqual(_conjugate(k_s, sector_clean), j_sector)
            expected_c = _product(*[
                _t_phase(scale(_tau(mask << sector_width, sector_width),
                               -1 if index == 1 and central_bit else 1))
                for index, mask in enumerate(projected_c)
            ])
            self.assertEqual(ec_sector, expected_c)
            common = _conjugate(e_minus, j_sector)
            self.assertEqual(common, _conjugate(_product(ec_sector, k_c_sector), sector_clean))
            self.assertEqual(_conjugate(h_axis, common), scale(common, -1))
            self.assertEqual(_conjugate(h_axis, ec_sector), ec_sector)
            self.assertEqual(_conjugate(h_axis, e_minus), e_minus)

            transition = _product(adjoint(k_c_sector), adjoint(ec_sector), e_minus, k_s)
            extended = _conjugate(adjoint(dressing), transition)
            self.assertEqual(_conjugate(extended, sector_clean), sector_clean)
            self.assertEqual(_conjugate(extended, _tau(1, sector_width)),
                             _tau(1, sector_width))
            for row in range(16):
                for column in range(16):
                    expected = (extended[row & ~1][column & ~1]
                                if row % 2 == column % 2 else CZ)
                    self.assertEqual(extended[row][column], expected)


if __name__ == '__main__':
    unittest.main()
