"""A paid four-address flat seed, six-T bank decoder and complete indicator.

Every native comparison retains scalar phases, arbitrary output and dirty
inputs, and the actual inverse. This finite pilot is not an all-size rule.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CO, CZ, Q2, Gate, add, adjoint, apply_word, eye, multiply, native_columns, scale,
)
from tests.test_field_copy_clifford_native import _fixed_query, _three_cube_query
from tests.test_field_copy_clifford_pilot import _pauli
from tests.test_field_shear_intertwiner import _native_x
from tests.test_joint_lookup_pauli_source import _bands, _inverse
from tests.test_one_band_copy_feasibility import _conjugate, _lift


SEED = {86: Q2(Fraction(1, 2)), 80: Q2(Fraction(1, 2)),
        113: Q2(Fraction(-1, 2)), 115: Q2(Fraction(-1, 2))}
ORBIT = [_conjugate(SEED, _fixed_query(j), 4) for j in range(4)]
AXES = ((33, 41), (2, 14), (35, 39))
TARGET = ({113: Q2(-1)}, {64: Q2(1)}, {89: Q2(-1)}, {96: Q2(1)})


def preparation_word():
    # V*=T(X0 X1) T(Y0 Z1) T(Z0), followed by the fixed field embedding.
    u = [('H', 1), ('CX', 1, 0), ('S', 0), ('S', 0), ('S', 0),
         ('H', 0), ('H', 1)]
    v = [('T', 0)] + _inverse(u) + [('T', 0), ('T', 1)] + u
    f = [('S', 1), ('H', 0), ('H', 1), ('CX', 0, 2), ('S', 2),
         ('S', 0), ('CX', 1, 0), ('H', 2), ('S', 0), ('CX', 2, 0),
         ('H', 2), ('S', 2), ('H', 1)]
    return v + f


def seed_word():
    v = preparation_word()
    return _inverse(v) + _native_x(0) + v


def diagonalizers():
    # K maps the listed physical Z axes to each pair in AXES, with + signs.
    k1 = [('CX', 0, 3), ('H', 0), ('H', 3),
          ('H', 1), ('CX', 0, 1), ('H', 1)]
    k2 = [('CX', 1, 2), ('H', 1), ('H', 2), ('CX', 2, 3)]
    k3 = [('CX', 1, 2), ('H', 1), ('S', 1), ('CX', 1, 0), ('H', 2)]
    return ((k1, (0, 3)), (k2, (1, 2)), (k3, (1, 2)))


def decoder_word():
    word = []
    for k, wires in diagonalizers():
        word += _inverse(k) + [('T', j) for j in wires] + k
    return word


def bank_clifford():
    # Sends the four signed TARGET Paulis to +Z0,+Z1,+Z2,+Z3.
    return [('S', 0), ('S', 0), ('S', 0), ('H', 0),
            ('CX', 1, 0), ('CX', 2, 0),
            ('H', 0), ('S', 0), ('S', 0), ('H', 0),
            ('CX', 1, 2), ('CX', 2, 1), ('CX', 1, 2), ('H', 3),
            ('CX', 0, 2), ('CX', 3, 2), ('CX', 0, 3),
            ('CX', 1, 3), ('CX', 2, 3)]


def full_decoder():
    # E is a complete Pauli-bank decoder; its residual query phase is allowed.
    return decoder_word() + bank_clifford()


def full_intertwiner():
    # The extra diagonal Clifford fixes the complementary Pauli orbit too.
    correction = [('S', 1)] * 3 + [('S', 3)]
    correction += [('H', 3), ('CX', 0, 3), ('CX', 2, 3), ('H', 3)]
    return full_decoder() + correction


def copy_word():
    # In the decoded coordinates H_j CX(j;4+j) H_j copies Z_j into X_Yj.
    e = full_decoder()
    canonical = []
    for j in range(4):
        canonical += [('H', j), ('CX', j, 4 + j), ('H', j)]
    return e + canonical + _inverse(e)


def source_word():
    # Address 0,1; physical output 2..5; arbitrary returned query masks 6..8.
    e = [(g[0], *(j + 2 for j in g[1:])) for g in full_decoder()]
    query, _ = _three_cube_query()
    return _inverse(e) + query + e


def indicator_word():
    source = source_word()
    wall = [('H', j) for j in range(2, 6)]
    return wall + _inverse(source) + [('S', 2), ('S', 2)] + source + wall


def copy_indicator_word():
    # Address 0,1; carrier 2..5; output 6..9; returned masks 10..12.
    query, _ = _three_cube_query()
    query = [(g[0], *(j + 4 if j >= 6 else j for j in g[1:])) for g in query]
    seed = [(g[0], *(j + 2 for j in g[1:])) for g in seed_word()]
    copy = [(g[0], *(j + 2 for j in g[1:])) for g in copy_word()]
    source = _inverse(query) + seed + query
    return _inverse(source) + _inverse(copy) + source + copy


def ledger(word):
    counts = Counter(g[0] for g in word)
    return (counts['T'] + counts['TDG'], len(_bands(word)), len(word),
            max(j for g in word for j in g[1:]) + 1)


def _copy_blocks():
    # Independent Fourier assembly from all sixteen output-X sectors.
    # The carrier action in sector s is E^dagger X^s E, including its scalar.
    e = native_columns(4, full_decoder())
    ed = adjoint(e)
    sectors = [multiply(ed, [e[row ^ s] for row in range(16)])
               for s in range(16)]
    return [[[sum((sectors[s][row][col]
                   * (-1 if (s & difference).bit_count() & 1 else 1)
                   for s in range(16)), CZ) * Fraction(1, 16)
              for col in range(16)] for row in range(16)]
            for difference in range(16)]


class FlatFieldBankNativeTests(unittest.TestCase):
    def test_paid_flat_seed_and_literal_inverse(self):
        v = preparation_word()
        self.assertEqual(_conjugate({1: Q2(1)}, v, 4), SEED)
        expected = [[CZ] * 16 for _ in range(16)]
        for label, coefficient in SEED.items():
            self.assertEqual(coefficient.b, 0)
            expected = add(expected, scale(_pauli(label & 15, label >> 4),
                                           coefficient.a))
        self.assertEqual(native_columns(4, seed_word()), expected)
        self.assertEqual(native_columns(4, _inverse(seed_word())), adjoint(expected))
        self.assertEqual(multiply(expected, expected), eye(16))
        self.assertEqual(native_columns(4, _inverse(v)),
                         adjoint(native_columns(4, v)))

    def test_three_physical_bands_and_complete_bank_decoding(self):
        for (k, wires), axes in zip(diagonalizers(), AXES):
            for wire, axis in zip(wires, axes):
                self.assertEqual(_conjugate({1 << (4 + wire): Q2(1)}, k, 4),
                                 {axis: Q2(1)})
        d, e = decoder_word(), full_decoder()
        for j, (reflection, target) in enumerate(zip(ORBIT, TARGET)):
            self.assertEqual(_conjugate(reflection, d, 4), target)
            self.assertEqual(_conjugate(target, _inverse(d), 4), reflection)
            self.assertEqual(_conjugate(reflection, e, 4),
                             {1 << (4 + j): Q2(1)})
        for word in (d, e):
            self.assertEqual(native_columns(4, _inverse(word)),
                             adjoint(native_columns(4, word)))

    def test_complete_copy_all_columns_and_actual_inverse(self):
        word = copy_word()
        for j, reflection in enumerate(ORBIT):
            lifted = {_lift(k, 4, 8): v for k, v in reflection.items()}
            target = {k ^ (1 << (4 + j)): v for k, v in lifted.items()}
            self.assertEqual(_conjugate(lifted, word, 8), target)
            self.assertEqual(_conjugate(target, _inverse(word), 8), lifted)
        blocks = _copy_blocks()
        for actual in (word, _inverse(word)):
            gates = [Gate(g[0], tuple(g[1:])) for g in actual]
            for column in range(256):
                carrier, output = column & 15, column >> 4
                expected = {row | (other << 4): blocks[other ^ output][row][carrier]
                            for other in range(16) for row in range(16)
                            if blocks[other ^ output][row][carrier] != CZ}
                self.assertEqual(apply_word(gates, {column: CO}), expected)

    def test_direct_indicator_every_input_and_actual_inverse(self):
        word = indicator_word()
        for actual in (word, _inverse(word)):
            gates = [Gate(g[0], tuple(g[1:])) for g in actual]
            for column in range(512):
                expected = column ^ (1 << (2 + (column & 3)))
                self.assertEqual(apply_word(gates, {column: CO}), {expected: CO})

    def test_native_ledgers_and_common_allocations(self):
        rows = ((preparation_word(), (3, 2, 28, 3), (3, 2, 38, 3)),
                (seed_word(), (6, 4, 70, 3), (6, 4, 70, 3)),
                (decoder_word(), (6, 3, 38, 4), (6, 3, 38, 4)),
                (bank_clifford(), (0, 0, 19, 4), (0, 0, 17, 4)),
                (full_decoder(), (6, 3, 57, 4), (6, 3, 55, 4)),
                (full_intertwiner(), (6, 3, 65, 4), (6, 3, 63, 4)),
                (copy_word(), (12, 6, 124, 8), (12, 6, 124, 8)),
                (source_word(), (36, 14, 170, 9), (36, 14, 170, 9)),
                (indicator_word(), (72, 28, 350, 9), (72, 28, 350, 9)),
                (copy_indicator_word(), (132, 52, 620, 13), (132, 52, 620, 13)))
        for word, forward, inverse in rows:
            self.assertEqual(ledger(word), forward)
            self.assertEqual(ledger(_inverse(word)), inverse)
            for actual in (word, _inverse(word)):
                self.assertLessEqual({g[0] for g in actual}, {'H', 'S', 'CX', 'T', 'TDG'})
                for band in _bands(actual):
                    self.assertEqual(len({g[1] for g in band}), len(band))
        self.assertEqual(ledger(indicator_word())[3] - 2 - 4, 3)
        self.assertEqual(ledger(copy_indicator_word())[3] - 2 - 4, 7)

    def test_encoded_echo_native_wiring_and_return(self):
        word = copy_indicator_word()
        for actual in (word, _inverse(word)):
            gates = [Gate(g[0], tuple(g[1:])) for g in actual]
            for address in range(4):
                column = address | ((3 + address) << 2) | (10 << 6) | ((7 - address) << 10)
                self.assertEqual(apply_word(gates, {column: CO}),
                                 {column ^ (1 << (6 + address)): CO})


if __name__ == '__main__':
    unittest.main()
