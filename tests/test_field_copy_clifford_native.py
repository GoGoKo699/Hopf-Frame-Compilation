"""Paid finite noncommuting field copy and its complete two-bit indicator.

The exact block checks include every arbitrary carrier/output input and
actual inverse. This four-row pilot is not a growing-address compiler.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CO, CZ, Gate, add, adjoint, apply_word, eye, multiply, native_columns, scale,
)
from tests.test_bent_query_baseline import _compact_ccz
from tests.test_dirty_indicator_depth import _exact_phase_action
from tests.test_field_copy_clifford_pilot import _pauli, _seed_and_orbit, _seed_words
from tests.test_field_shear_intertwiner import _native_x, _p2_native_decoder, _pencil_columns, _shear
from tests.test_joint_lookup_pauli_source import _bands, _inverse
from tests.test_t_depth import _word


_SIGN_PARTNERS = ((1, 3), (0, 2), (3, 1), (2, 0))


def _coordinate_word():
    # K(H0 H1)K^-1=R and K^-1 P K=-X3. Every macro is
    # Hermitian and its displayed expansion retains literal phases.
    return ([('CX', 0, 1), ('CX', 1, 0), ('CX', 0, 1),
             ('CX', 2, 3), ('CX', 3, 2), ('CX', 2, 3), ('CX', 3, 1),
             ('H', 1), ('CX', 0, 1), ('H', 1),
             ('H', 2), ('CX', 0, 2), ('H', 2)]
            + [('S', 1)] * 2 + [('S', 2)] * 2)


def _fixed_query(address):
    return [('CX', source, 2 + target)
            for source, column in enumerate(_pencil_columns(2, address))
            for target in range(2) if column >> target & 1]


def _y_rotation_pair():
    # F=SH on both targets maps Z to Y. U=F T F^-1 maps Z0Z1
    # to H0H1, with the native scalar retained and its actual inverse used.
    f = [('H', 0), ('S', 0), ('H', 1), ('S', 1)]
    return _inverse(f) + [('T', 0), ('T', 1)] + f


def _controlled_hh(output=4):
    # The arbitrary carrier wire 2 is a returned mask for the CCZ.
    u = _y_rotation_pair()
    return (_inverse(u) + [('CX', 0, 1)]
            + _word(_compact_ccz(output, 3, 1, 2))
            + [('CX', 0, 1)] + u)


def _sector_factor(address, positive, output=4):
    k = _coordinate_word() + _fixed_query(address)
    # K_x^-1 P K_x=-X3; after H3, positive P means bit 3=1.
    flip = [] if positive else [('H', 3), ('S', 3), ('S', 3), ('H', 3)]
    return (_inverse(k) + [('H', 3)] + flip + _controlled_hh(output)
            + flip + [('H', 3)] + k)


def _copy_word():
    wall = [('H', 4 + j) for j in range(4)]
    word = []
    for j, (positive, negative) in enumerate(_SIGN_PARTNERS):
        word += _sector_factor(positive, True, 4 + j)
        word += _sector_factor(negative, False, 4 + j)
    return wall + word + wall


def _three_cube_query():
    # x=(0,1), u=(2,3), v=(4,5), private arbitrary masks=(6,7,8).
    # f=x0*u0*(v0+v1)+x1*u1*v0+v1*(x0+x1)*(u0+u1).
    terms = (((0,), (2,), (4, 5)), ((1,), (3,), (4,)),
             ((5,), (0, 1), (2, 3)))
    helpers = (6, 7, 8)
    core = []
    old = 0
    for step in range(8):
        vertex = step ^ (step >> 1)
        difference = old ^ vertex
        if difference:
            coordinate = difference.bit_length() - 1
            core += [('CX', wire, helper)
                     for term, helper in zip(terms, helpers)
                     for wire in term[coordinate]]
        core += [('TDG' if vertex.bit_count() & 1 else 'T', helper)
                 for helper in helpers]
        old = vertex
    core += [('CX', wire, helper)
             for term, helper in zip(terms, helpers) for wire in term[2]]
    wall = [('H', 4), ('H', 5)]
    return wall + core + wall, core


def _indicator_word():
    # Physical address 0,1; dirty carrier 2..5; output 6..9;
    # private returned query masks 10..12.
    query, _ = _three_cube_query()
    query = [(g[0], *(wire + 4 if wire >= 6 else wire for wire in g[1:]))
             for g in query]
    seed = [(g[0], *(wire + 2 for wire in g[1:])) for g in _seed_words()[0]]
    copy = [(g[0], *(wire + 2 for wire in g[1:])) for g in _copy_word()]
    source = _inverse(query) + seed + query
    return _inverse(source) + _inverse(copy) + source + copy


def _decoder_indicator_word():
    # The improved query is also available to the existing direct decoder.
    # Address 0,1; physical output 2..5; returned dirty masks 6..8.
    decoder = [(g[0], *(wire + 2 for wire in g[1:]))
               for g in _p2_native_decoder()]
    query, _ = _three_cube_query()
    source = _inverse(decoder) + query + decoder
    return _inverse(source) + _native_x(2) + source


class FieldCopyCliffordNativeTests(unittest.TestCase):
    def test_native_control_core_and_all_sector_factors(self):
        r, orbit = _seed_and_orbit()
        p = _pauli(4, 1)
        core = _controlled_hh()
        core_matrix = native_columns(5, core)
        hh = native_columns(4, [('H', 0), ('H', 1)])
        expected = eye(32)
        for col in range(32):
            if col & 16 and col & 8:
                for row in range(32):
                    expected[row][col] = hh[row & 15][col & 15] if row & 16 else CZ
        self.assertEqual(core_matrix, expected)
        self.assertEqual(native_columns(5, _inverse(core)), adjoint(expected))
        k = native_columns(4, _coordinate_word())
        self.assertEqual(multiply(multiply(k, hh), adjoint(k)), r)
        for address in range(4):
            kx = native_columns(4, _coordinate_word() + _fixed_query(address))
            self.assertEqual(multiply(multiply(adjoint(kx), p), kx),
                             scale(_pauli(8, 0), -1))
            for positive in (True, False):
                word = _sector_factor(address, positive)
                actual = native_columns(5, word)
                self.assertEqual(native_columns(5, _inverse(word)), adjoint(actual))
                projector = scale(add(eye(16), scale(p, 1 if positive else -1)),
                                  Fraction(1, 2))
                wanted = add(add(eye(16), scale(projector, -1)),
                             multiply(projector, orbit[address]))
                for col in range(32):
                    for row in range(32):
                        target = ((CO if row == col else CZ) if not col & 16
                                  else wanted[row & 15][col & 15] if row & 16 else CZ)
                        self.assertEqual(actual[row][col], target)

    def test_complete_copy_all_output_x_sectors_and_actual_inverse(self):
        _, orbit = _seed_and_orbit()
        blocks = {}
        for address in range(4):
            for positive in (True, False):
                actual = native_columns(5, _sector_factor(address, positive))
                blocks[address, positive] = [row[16:] for row in actual[16:]]
        not_involution = False
        # The emitted output-H walls convert these exact Z-control blocks
        # to all sixteen output-X sectors. Off-sector coherences are also
        # covered by equality of complete block-diagonal operators.
        for mask in range(16):
            action = eye(16)
            reverse = eye(16)
            for j, (a, b) in enumerate(_SIGN_PARTNERS):
                if mask >> j & 1:
                    for key in ((a, True), (b, False)):
                        action = multiply(blocks[key], action)
                        reverse = multiply(reverse, adjoint(blocks[key]))
            self.assertEqual(reverse, adjoint(action))
            self.assertEqual(multiply(action, reverse), eye(16))
            not_involution |= action != reverse
            for j, reflection in enumerate(orbit):
                self.assertEqual(multiply(multiply(action, reflection), reverse),
                                 scale(reflection, -1 if mask >> j & 1 else 1))
        self.assertTrue(not_involution)

    def test_three_cube_query_all_principal_and_dirty_columns(self):
        query, core = _three_cube_query()
        for word in (core, _inverse(core)):
            for basis in range(512):
                x, u, v = basis & 3, (basis >> 2) & 3, (basis >> 4) & 3
                cols = _pencil_columns(2, x)
                product = (cols[0] if u & 1 else 0) ^ (cols[1] if u & 2 else 0)
                phase = 4 * ((product & v).bit_count() & 1)
                self.assertEqual(_exact_phase_action(basis, word), (basis, phase))
        for word in (query, _inverse(query)):
            gates = [Gate(g[0], tuple(g[1:])) for g in word]
            for basis in range(512):
                address, carrier = basis & 3, (basis >> 2) & 15
                expected = (basis & ~60) | (_shear(2, address, carrier) << 2)
                self.assertEqual(apply_word(gates, {basis: CO}), {expected: CO})

    def test_finite_native_ledgers_and_common_allocation(self):
        query, _ = _three_cube_query()
        copy = _copy_word()
        indicator = _indicator_word()
        decoder_indicator = _decoder_indicator_word()
        expected = {
            'query': (24, 8, 58, 9),
            'copy': (96, 32, 792, 8),
            'indicator': (288, 96, 1868, 13),
            'decoder_indicator': (300, 160, 892, 9),
        }
        for name, word in (('query', query), ('copy', copy), ('indicator', indicator),
                           ('decoder_indicator', decoder_indicator)):
            for actual in (word, _inverse(word)):
                counts = Counter(g[0] for g in actual)
                t, depth, gates, width = expected[name]
                self.assertEqual(counts['T'] + counts['TDG'], t)
                self.assertEqual(len(_bands(actual)), depth)
                self.assertEqual(len(actual), gates)
                self.assertEqual({wire for g in actual for wire in g[1:]}, set(range(width)))
                for band in _bands(actual):
                    self.assertEqual(len({g[1] for g in band}), len(band))
        # Four arbitrary carrier wires and three returned query masks are
        # private to the final six-wire indicator: seven dirty wires total.
        self.assertEqual(13 - 2 - 4, 7)

    def test_assembled_indicator_wiring_and_return_in_both_directions(self):
        # Component identities prove every column. These independent full
        # native-word executions additionally audit the common allocation,
        # output bit order, chronology and actual inverse wiring.
        word = _indicator_word()
        for actual in (word, _inverse(word)):
            gates = [Gate(g[0], tuple(g[1:])) for g in actual]
            for address in range(4):
                column = (address | ((3 + address) << 2)
                          | (10 << 6) | ((7 - address) << 10))
                self.assertEqual(apply_word(gates, {column: CO}),
                                 {column ^ (1 << (6 + address)): CO})
        decoder = _decoder_indicator_word()
        for actual in (decoder, _inverse(decoder)):
            gates = [Gate(g[0], tuple(g[1:])) for g in actual]
            for address in range(4):
                column = address | (10 << 2) | ((7 - address) << 6)
                self.assertEqual(apply_word(gates, {column: CO}),
                                 {column ^ (1 << (2 + address)): CO})


if __name__ == '__main__':
    unittest.main()
