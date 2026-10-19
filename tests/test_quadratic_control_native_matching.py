"""Exact guards for quadratic controls; no hard-query counterexample is claimed."""
from collections import Counter
import unittest

from tests.exact_matrices import CO, CZ, adjoint, eye, multiply, native_columns
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _clifford_columns, _conjugate, _tau


def _controlled_x(width, target, control):
    result = [[CZ for _ in range(1 << width)] for _ in range(1 << width)]
    for value in range(1 << width):
        result[value ^ ((1 << target) if control(value) else 0)][value] = CO
    return result


def _family(blocks, with_clean):
    offset = int(with_clean)
    h = offset
    c_mask = 1 if with_clean else 0
    changes, phases, lists = [[], [], []], [[], [], []], [[], [], []]
    terms = []
    for j in range(blocks):
        x, y = offset + 1 + 2 * j, offset + 2 + 2 * j
        t = 1 if j < blocks // 2 else -1
        xm, ym, hm = 1 << x, 1 << y, 1 << h
        aa = [(c_mask ^ xm, t), (c_mask ^ hm ^ ym, -t)]
        bb = [(c_mask ^ ym, t), (c_mask ^ xm ^ ym, -t)]
        cc = [(c_mask ^ hm ^ xm, t), (c_mask ^ hm ^ xm ^ ym, -t)]
        for index, axes in enumerate((aa, bb, cc)):
            lists[index].extend(axis for axis, _ in axes)
        terms.extend(aa + bb + [(axis, -sign) for axis, sign in cc])
        if with_clean:
            changes[0].extend([('CX', 0, x), ('CX', 0, y)])
        changes[0].append(('CX', h, y))
        changes[1].append(('CX', y, x))
        if with_clean:
            changes[1].append(('CX', 0, y))
        changes[2].append(('CX', x, y))
        if with_clean:
            changes[2].append(('CX', 0, x))
        changes[2].append(('CX', h, x))
        for index, targets in enumerate(((x, y), (y, x), (x, y))):
            phases[index].extend([('T' if t == 1 else 'TDG', targets[0]),
                                  ('TDG' if t == 1 else 'T', targets[1])])
    bands = [_inverse(change) + phase + change
             for change, phase in zip(changes, phases)]
    endpoint = [('H', 0)] if with_clean else []
    word = endpoint + bands[0] + bands[1] + _inverse(bands[2]) + endpoint
    return word, lists, terms


class QuadraticControlNativeMatchingTests(unittest.TestCase):
    def test_literal_affine_absorption_and_output_band_conjugation(self):
        # The test is the full operator repacking identity, not a diagonal target.
        # H=X0, b=z1*z2+z1+1; its affine factor is actual X0 CNOT1,0.
        affine_word = [('H', 0), ('S', 0), ('S', 0), ('H', 0), ('CX', 1, 0)]
        cl = native_columns(3, affine_word)
        self.assertEqual(native_columns(3, _inverse(affine_word)), adjoint(cl))
        _clifford_columns(cl, 3)
        fq = _controlled_x(3, 0, lambda v: ((v >> 1) & 1) & ((v >> 2) & 1))
        f = _controlled_x(3, 0, lambda v: (((v >> 1) & 1) & ((v >> 2) & 1))
                          ^ ((v >> 1) & 1) ^ 1)
        self.assertEqual(f, multiply(cl, fq))
        ec = native_columns(3, [('T', 0)])
        eb = native_columns(3, [('H', 0), ('T', 0), ('H', 0)])
        ea = native_columns(3, [('T', 2)])
        self.assertEqual(multiply(cl, eb), multiply(eb, cl))
        ec_prime = _conjugate(adjoint(cl), ec)
        self.assertNotEqual(ec_prime, ec)
        original = multiply(adjoint(ec), multiply(f, multiply(eb, ea)))
        repacked = multiply(cl, multiply(adjoint(ec_prime), multiply(fq, multiply(eb, ea))))
        self.assertEqual(repacked, original)
        self.assertEqual(multiply(adjoint(original), original), eye(8))

    def test_rank_two_support_boundary_is_an_actual_native_band(self):
        # C maps Z1,Z2 to X0 Z1,X0 Z2. Only the two T targets form the band.
        change = [('CX', 0, 1), ('CX', 0, 2), ('H', 0)]
        word = _inverse(change) + [('T', 1), ('T', 2)] + change
        e = native_columns(3, word)
        self.assertEqual(native_columns(3, _inverse(word)), adjoint(e))
        h = _tau(1 << 3, 3)
        self.assertNotEqual(multiply(e, h), multiply(h, e))
        controlled_h = eye(8)
        for value in range(8):
            b = ((value >> 1) & 1) & (1 ^ ((value >> 2) & 1))
            controlled_h[value][value] = -CO if b and value & 1 else CO
        self.assertEqual(_conjugate(adjoint(e), controlled_h), controlled_h)
        c = _conjugate(adjoint(e), h)
        _clifford_columns(c, 3)
        self.assertTrue(any(c[row][column] != CZ
                            for row in range(8) for column in range(8) if row != column))
        omega = native_columns(1, [('T', 0)])[1][1]
        for value in (2, 3):
            for row in range(8):
                self.assertEqual(e[row][value], omega if row == value else CZ)

    def test_actual_family_all_columns_ranks_and_native_ledgers(self):
        for blocks in (2, 4):
            for with_clean in (False, True):
                word, lists, terms = _family(blocks, with_clean)
                width = 2 * blocks + 1 + int(with_clean)
                self.assertEqual([_rank(axes) for axes in lists], [2 * blocks] * 3)
                self.assertEqual([len(band) for band in _native_t_bands(word)],
                                 [2 * blocks] * 3)
                for actual_word in (word, _inverse(word)):
                    counts = Counter(gate[0] for gate in actual_word)
                    self.assertEqual(counts['T'] + counts['TDG'], 6 * blocks)
                    self.assertEqual(counts['H'] + counts['CX'],
                                     16 * blocks + 2 if with_clean else 8 * blocks)
                for value in range(1 << width):
                    data = value >> int(with_clean)
                    bit = (data & 1) & (sum(((data >> (1 + 2*j)) & 1)
                                             * ((data >> (2 + 2*j)) & 1)
                                             for j in range(blocks)) % 2)
                    exponent = sum(sign * ((axis & value).bit_count() % 2)
                                   for axis, sign in terms) % 8
                    self.assertEqual(exponent, 4 * bit)
                # Only the smallest member uses a dense exact matrix.
                if blocks == 2:
                    actual = native_columns(width, word)
                    self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
                    for column in range(1 << width):
                        exponent = sum(sign * ((axis & column).bit_count() % 2)
                                       for axis, sign in terms) % 8
                        for row in range(1 << width):
                            self.assertEqual(actual[row][column],
                                             (-CO if exponent else CO)
                                             if row == column else CZ)
            # In actual shared-center coordinates, b=sum(zj*yj+yj).
            polar_rows = [1 << (j ^ 1) for j in range(2 * blocks)]
            self.assertEqual(_rank(polar_rows), 2 * blocks)


if __name__ == '__main__':
    unittest.main()
