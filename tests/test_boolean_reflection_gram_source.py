"""Exact guards for paid Boolean spectral reshaping and its flat-mixer gap.

The fixtures test literal phases and arbitrary borrowed columns.  They
do not supply the unpaid nonflat reflection factorization or a scalable
complete indicator.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, native_columns,
    scale, zeros,
)
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _ccz, _counts
from tests.test_signal_interleaved_kernel_source import (
    _ct_word, _one_clean_block, _root_powers,
)
from tests.test_t_wall_xor_gram_source import _p1_mixer


HALF = CQ2(Fraction(1, 2))
QUARTER = CQ2(Fraction(1, 4))


def _cs(a, b):
    return [('T', a), ('T', b), ('CX', a, b), ('TDG', b), ('CX', a, b)]


def _ct(a, b, helpers):
    mapping = dict(enumerate((a, b, *helpers)))
    return [tuple([gate[0]] + [mapping[wire] for wire in gate[1:]])
            for gate in _ct_word()]


def _qft8_word():
    helpers = (3, 4, 5, 6)
    return ([('H', 2)] + _cs(1, 2) + _ct(0, 2, helpers)
            + [('H', 1)] + _cs(0, 1) + [('H', 0)]
            + [('CX', 0, 2), ('CX', 2, 0), ('CX', 0, 2)])


def _fourier8():
    powers = _root_powers()
    normalization = CQ2(Q2(0, Fraction(1, 4)))
    return [[normalization * powers[(row * column) % 8]
             for column in range(8)] for row in range(8)]


def _quarter_phase_echo():
    # data bit 0 is f, signal 1 is a, and 2,3 are arbitrary u,v.
    return ([('CX', 0, 2)] + _inverse(_cs(1, 2))
            + [('CX', 0, 2)] + _cs(1, 2)
            + [('CX', 0, 3)] + _ccz((1, 2, 3))
            + [('CX', 0, 3)] + _ccz((1, 2, 3)))


def _quarter_phase(reflection):
    plus = scale(add(eye(len(reflection)), reflection), HALF)
    minus = scale(add(eye(len(reflection)), scale(reflection, -CO)), HALF)
    e = add(plus, scale(minus, -CI))
    v = add(scale(plus, CI), scale(minus, _root_powers()[1]))
    return e, v


def _boolean_parts(address, f, g):
    raw_f, raw_g = zeros(16), zeros(16)
    for value in range(16):
        z, t = value & 3, value >> 2
        shifted = t ^ (z & 2)
        h = ((shifted ^ (address * (z >> 1))).bit_count()
             - shifted.bit_count()) % 8
        raw_f[value][value] = -CO if f(h) else CO
        raw_g[value][value] = -CO if g(h) else CO
    mixer = native_columns(4, _p1_mixer())
    r = multiply(multiply(adjoint(mixer), raw_f), mixer)
    s = multiply(native_columns(4, [('Z', 0), ('Z', 2)]), raw_g)
    return r, s


class BooleanReflectionGramSourceTests(unittest.TestCase):
    def test_literal_qft8_has_correct_orientation_and_returns_all_helpers(self):
        word = _qft8_word()
        self.assertEqual(_counts(word), (76, 98))
        actual = native_columns(7, word)
        expected = _fourier8()
        for column in range(128):
            for row in range(128):
                value = expected[row & 7][column & 7] if row >> 3 == column >> 3 else CZ
                self.assertEqual(actual[row][column], value)
        self.assertEqual(native_columns(7, _inverse(word)), adjoint(actual))

    def test_exact_phase_evaluation_and_truth_lookup_return_clean_selectors(self):
        fourier = _fourier8()
        hadamards = native_columns(3, [('H', bit) for bit in range(3)])
        powers = _root_powers()
        for h in range(8):
            phase = zeros(8)
            for selector in range(8):
                phase[selector][selector] = powers[(h * selector) % 8]
            w = multiply(multiply(adjoint(fourier), phase), hadamards)
            self.assertEqual([row[0] for row in w],
                             [CO if row == h else CZ for row in range(8)])
            extended = zeros(16)
            for output in (0, 1):
                for row in range(8):
                    for column in range(8):
                        extended[row + 8 * output][column + 8 * output] = w[row][column]
            for table in (lambda k: k >> 2, lambda k: k & 1,
                          lambda k: int(k in (2, 3, 7))):
                lookup = zeros(16)
                for column in range(16):
                    lookup[column ^ (8 * table(column & 7))][column] = CO
                source = multiply(multiply(adjoint(extended), lookup), extended)
                for output in (0, 1):
                    column = 8 * output
                    target = 8 * (output ^ table(h))
                    self.assertEqual([row[column] for row in source],
                                     [CO if row == target else CZ for row in range(16)])

    def test_quarter_phase_echo_retains_literal_scalar_and_arbitrary_dirty_bits(self):
        word = _quarter_phase_echo()
        self.assertEqual(_counts(word), (20, 28))  # includes four XOR-query CNOTs
        actual = native_columns(4, word)
        for column in range(16):
            phase = -CI if column & 3 == 3 else CO
            for row in range(16):
                self.assertEqual(actual[row][column], phase if row == column else CZ)
        self.assertEqual(native_columns(4, _inverse(word)), adjoint(actual))

    def test_general_quarter_phase_normal_form_keeps_final_unitary_conjugation(self):
        r = native_columns(1, [('H', 0)])
        f = native_columns(1, [('T', 0)])
        e, v = _quarter_phase(r)
        accepted = _one_clean_block(e, f)
        product = multiply(f, r)
        expected = scale(multiply(add(eye(2), product), v),
                         CQ2(Q2(0, Fraction(1, 4))))
        self.assertEqual(accepted, expected)
        real = scale(add(product, adjoint(product)), HALF)
        expected_gram = scale(multiply(multiply(adjoint(v), add(eye(2), real)), v),
                              QUARTER)
        self.assertEqual(multiply(adjoint(accepted), accepted), expected_gram)

    def test_two_reflection_gram_and_complete_p1_opposite_charge_witnesses(self):
        tables = ((lambda h: h >> 2, lambda h: h >> 2, CQ2(Fraction(3, 8))),
                  (lambda h: h & 1, lambda h: int(h == 7), QUARTER))
        for f, g, expected_witness in tables:
            for address in (0, 1):
                r, s = _boolean_parts(address, f, g)
                self.assertEqual(multiply(r, r), eye(16))
                self.assertEqual(multiply(s, s), eye(16))
                e, _ = _quarter_phase(r)
                accepted = _one_clean_block(e, s)
                gram = multiply(adjoint(accepted), accepted)
                anticommutator = scale(add(multiply(s, r), multiply(r, s)), HALF)
                expected = scale(add(eye(16), anticommutator), QUARTER)
                self.assertEqual(gram, expected)
                self.assertEqual(sum(gram[index][index] for index in range(16)),
                                 CQ2(4))
                if address:
                    # z=0,t=0 and z=0,t=e1 have Q charges +1 and -1.
                    self.assertEqual(gram[0][0], expected_witness)
                    self.assertEqual(gram[8][8], expected_witness)
                else:
                    for row in range(16):
                        for column in range(16):
                            equal = not ((column ^ (column >> 2)) & 1)
                            value = HALF if row == column and equal else CZ
                            self.assertEqual(gram[row][column], value)

    def test_commuting_reflection_product_gives_the_exact_half_projector(self):
        r = native_columns(2, [('Z', 0)])
        s = native_columns(2, [('Z', 1)])
        e, _ = _quarter_phase(r)
        accepted = _one_clean_block(e, s)
        target = scale(add(eye(4), multiply(s, r)), QUARTER)
        self.assertEqual(multiply(adjoint(accepted), accepted), target)


if __name__ == '__main__':
    unittest.main()
