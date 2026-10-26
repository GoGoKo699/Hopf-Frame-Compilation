"""Reduced pair-bank sources and the strictly weaker carrier-linear probe.

The exact fixtures pay one native factor and check full arbitrary-input
identities. They do not supply a shallow all-pair compiler.
"""
from collections import Counter
from itertools import combinations
import unittest

from tests.exact_matrices import CO
from tests.test_coherent_feedback_sources import _shared_ccz
from tests.test_common_probe_source import _native
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _pairs(p):
    return tuple((1 << i) | (1 << j) for i, j in combinations(range(p), 2))


def _lower(c, subset, p):
    return sum(((c >> (t | subset)) & 1) << t
               for t in range(1 << p) if not t & subset)


def _f(s, c, p):
    for bit, pair in enumerate(_pairs(p)):
        if s >> bit & 1:
            c ^= _lower(c, pair, p)
    return c


def _clique(x, p):
    return sum(int(x & pair == pair) << bit for bit, pair in enumerate(_pairs(p)))


def _j(x, c, p):
    result = c
    for i in range(p):
        if x >> i & 1:
            result ^= _lower(c, 1 << i, p)
    return result


def _l(x, c, p):
    return c ^ (sum(((x >> i) & (c >> (1 << i)) & 1) for i in range(p)) & 1)


def _linear_completion(p):
    root = p
    return ([('H', root)]
            + _shared_ccz(root, [(i, p + (1 << i)) for i in range(p)])
            + [('H', root)])


def _rows(columns):
    return tuple(sum(((column >> row) & 1) << j for j, column in enumerate(columns))
                 for row in range(len(columns)))


def _apply(rows, c):
    return sum(((row & c).bit_count() & 1) << bit for bit, row in enumerate(rows))


def _invert(rows):
    width = len(rows)
    left, right = list(rows), [1 << i for i in range(width)]
    for bit in range(width):
        pivot = next(i for i in range(bit, width) if left[i] >> bit & 1)
        left[bit], left[pivot] = left[pivot], left[bit]
        right[bit], right[pivot] = right[pivot], right[bit]
        for i in range(width):
            if i != bit and left[i] >> bit & 1:
                left[i] ^= left[bit]
                right[i] ^= right[bit]
    assert left == [1 << i for i in range(width)]
    return tuple(right)


def _interpolate(a, b, p):
    values = [a[_clique(x, p)] for x in range(1 << p)]
    for x in range(1 << p):
        for i in range(p):
            if x >> i & 1:
                values[x] ^= b[i][_clique(x, p)]
    for i in range(p):
        for x in range(1 << p):
            if x >> i & 1:
                values[x] ^= values[x ^ (1 << i)]
    return tuple(values)


def _weak_u(s):
    # Pair order 12,13,14,23,24,34; coefficient labels are subset masks.
    s12, s13 = s & 1, (s >> 1) & 1
    return ((s12 * s13) << 0) | (s13 << 3) | (s12 << 5) | (1 << 7)


def _weak_h(s):
    columns = [_f(s, 1 << column, 4) for column in range(16)]
    columns[-1] ^= _weak_u(s)
    return _rows(columns)


def _factor(p, pair_index):
    q, pair = len(_pairs(p)), _pairs(p)[pair_index]
    matching = [(q + (t | pair), q + t)
                for t in range(1 << p) if not t & pair]
    wall = [('H', target) for _, target in matching]
    return wall + _shared_ccz(pair_index, matching) + list(reversed(wall))


class ReducedPairGaugeTests(unittest.TestCase):
    def test_reduced_echo_is_complete_for_every_pair_word(self):
        # Full carrier inputs at p=2,3; all columns suffice by linearity at p=4.
        for p in (2, 3, 4):
            q, m = len(_pairs(p)), 1 << p
            carriers = range(1 << m) if p < 4 else [1 << bit for bit in range(m)]
            for s in range(1 << q):
                for k in range(1 << q):
                    for c in carriers:
                        actual = _f(s ^ k, _f(s, c, p), p)
                        self.assertEqual(actual, _f(k, c, p))
                        self.assertEqual(_f(s, _f(s ^ k, actual, p), p), c)
            for x in range(m):
                for c in carriers:
                    actual = _j(x, _f(_clique(x, p), c, p), p)
                    direct = c
                    for i in range(p):
                        if x >> i & 1:
                            direct ^= _lower(direct, 1 << i, p)
                    self.assertEqual(actual, direct)

    def test_literal_shared_factor_and_actual_native_inverse(self):
        # p=3 has two disjoint carrier pairs sharing one arbitrary control.
        for p in (2, 3):
            q, m = len(_pairs(p)), 1 << p
            word = _factor(p, 0)
            active = sorted({wire for gate in word for wire in gate[1:]})
            mapping = {wire: i for i, wire in enumerate(active)}
            compact = [(gate[0], *(mapping[wire] for wire in gate[1:])) for gate in word]
            inverse = _inverse(compact)
            for basis in range(1 << len(active)):
                expanded = sum(((basis >> i) & 1) << wire for i, wire in enumerate(active))
                s, c = expanded & ((1 << q) - 1), expanded >> q
                output = s | (_f(s & 1, c, p) << q)
                target = sum(((output >> wire) & 1) << i for i, wire in enumerate(active))
                self.assertEqual(_native(compact, {basis: CO}), {target: CO})
                self.assertEqual(_native(inverse, {target: CO}), {basis: CO})
            r = m // 4
            for native in (compact, inverse):
                counts = Counter(gate[0] for gate in native)
                t = counts['T'] + counts['TDG']
                self.assertEqual(t, 6 * r + (r % 2))
                self.assertLessEqual(len(native) - t, 12 * r + 3)
                self.assertLessEqual(len(_bands(native)), 4)

        # Larger batches check the reduced root S power and both native ledgers.
        for p in (4, 5, 6):
            r = (1 << p) // 4
            for native in (_factor(p, 0), _inverse(_factor(p, 0))):
                counts = Counter(gate[0] for gate in native)
                t = counts['T'] + counts['TDG']
                self.assertEqual(t, 6 * r + (r % 2))
                self.assertLessEqual(len(native) - t, 12 * r + 3)
                self.assertLessEqual(len(_bands(native)), 4)

    def test_literal_root_completion_includes_all_inputs_and_actual_inverse(self):
        for p in (2, 3):
            word = _linear_completion(p)
            active = sorted({wire for gate in word for wire in gate[1:]})
            mapping = {wire: i for i, wire in enumerate(active)}
            compact = [(gate[0], *(mapping[wire] for wire in gate[1:])) for gate in word]
            inverse = _inverse(compact)
            for basis in range(1 << len(active)):
                expanded = sum(((basis >> i) & 1) << wire for i, wire in enumerate(active))
                x, c = expanded & ((1 << p) - 1), expanded >> p
                output = x | (_l(x, c, p) << p)
                target = sum(((output >> wire) & 1) << i for i, wire in enumerate(active))
                self.assertEqual(_native(compact, {basis: CO}), {target: CO})
                self.assertEqual(_native(inverse, {target: CO}), {basis: CO})
                self.assertEqual(_l(x, c, p) & 1, _j(x, c, p) & 1)
        for p in range(2, 10):
            for native in (_linear_completion(p), _inverse(_linear_completion(p))):
                counts = Counter(gate[0] for gate in native)
                t = counts['T'] + counts['TDG']
                self.assertEqual(t, 6 * p + (p % 2))
                self.assertLessEqual(len(native) - t, 10 * p + 5)
                self.assertLessEqual(len(_bands(native)), 4)

    def test_interpolation_exhausts_the_formal_probe_constraint(self):
        p, q, m = 3, 3, 8
        # Deliberately arbitrary row functions, with no invertibility promise.
        a = [((s * 29) ^ (s << 4) ^ 17) & 255 for s in range(1 << q)]
        b = [[((s + 3 * i) * 53 ^ (s << i) ^ (1 << i)) & 255
              for s in range(1 << q)] for i in range(p)]
        family = []
        for s in range(1 << q):
            family.append(_interpolate([a[s ^ k] for k in range(1 << q)],
                                       [[row[s ^ k] for k in range(1 << q)] for row in b], p))
        for s, rows in enumerate(family):
            self.assertEqual(rows[0], a[s])
            for i in range(p):
                self.assertEqual(rows[1 << i], b[i][s])
            for x in range(m):
                other = family[s ^ _clique(x, p)]
                left = other[0]
                for i in range(p):
                    if x >> i & 1:
                        left ^= other[1 << i]
                right = 0
                for subset in range(m):
                    if subset & x == subset:
                        right ^= rows[subset]
                self.assertEqual(left, right)

    def test_invertible_weak_gauge_has_correct_probe_but_not_full_transport(self):
        p, q, m = 4, 6, 16
        family = [_weak_h(s) for s in range(1 << q)]
        inverses = [_invert(rows) for rows in family]
        for s in range(1 << q):
            for c in (1 << bit for bit in range(m)):
                self.assertEqual(_apply(inverses[s], _apply(family[s], c)), c)
            for x in range(m):
                k = _clique(x, p)
                for c in (1 << bit for bit in range(m)):
                    # This is V_H dagger, including its actual inverse H_s^-1.
                    moved = _l(x, _apply(family[s ^ k], _apply(inverses[s], c)), p)
                    target_sign = sum((c >> subset) & 1 for subset in range(m)
                                      if subset & x == subset) & 1
                    self.assertEqual(moved & 1, target_sign)
                    returned = _apply(family[s], _apply(inverses[s ^ k], _l(x, moved, p)))
                    self.assertEqual(returned, c)
        self.assertNotEqual(_weak_u(1), _f(1, _weak_u(0), p))
        self.assertNotEqual(_apply(family[1], 1 << 15),
                            _f(1, _apply(family[0], 1 << 15), p))

    def test_low_degree_rows_force_a_zero_top_row_before_invertibility(self):
        # Each linear pair function and each singleton-weighted one has degree <=3.
        # At p=4 this checks an entire basis for every possible degree-one entry.
        p, q = 4, 6
        functions = [[1] * (1 << q)] + [
            [(s >> j) & 1 for s in range(1 << q)] for j in range(q)]
        zero = [0] * (1 << q)
        for location in range(p + 1):
            for function in functions:
                a, b = zero, [zero] * p
                if location == 0:
                    a = function
                else:
                    b = list(b)
                    b[location - 1] = function
                for s in range(1 << q):
                    rows = _interpolate([a[s ^ k] for k in range(1 << q)],
                                        [[row[s ^ k] for k in range(1 << q)] for row in b], p)
                    self.assertEqual(rows[-1], 0)
        # The canonical quadratic low row reaches the missing quartic term.
        self.assertEqual(_f(0b111111, 1 << 15, p) & 1, 1)


if __name__ == '__main__':
    unittest.main()
