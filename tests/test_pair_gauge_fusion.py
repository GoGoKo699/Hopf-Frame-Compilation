"""Complete small-space checks for pair-gauge probes and paid native pieces.

Full encoded-Hadamard failures and weaker root-probe contracts are tested
separately. None of these fixtures supplies an all-pair shallow compiler.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import unittest

from tests.exact_matrices import CI, CO
from tests.test_coherent_feedback_sources import _shared_ccz
from tests.test_common_probe_source import ROOT_HALF, _native
from tests.test_complete_observable_bases import (
    _hadamard_numerator, _hrt_unit, _permutation, _quadratic_inverse, _ring_product,
)
from tests.test_dirty_cubic_phase_compiler import _gray_phase_batch, _dirty_cubic_phase
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _pairs(r, first=0):
    return tuple((1 << i) | (1 << j) for i, j in combinations(range(first, r), 2))


def _mask(family):
    return sum(1 << s for s in family)


def _b(y, family):
    for s in family:
        i = s & -s
        j = s ^ i
        y ^= (((y >> i) & (y >> j)) & 1) << s
    return y


def _c(y, r, family):
    return _ring_product(_hrt_unit(y & _mask(family), r), y, r)


def _t(y, r, family, inverse=False):
    z = y & _mask(family)
    if inverse:
        return _ring_product(_hrt_unit(z, r), 1 ^ y, r) ^ 1 ^ z
    return _ring_product(_hrt_unit(z, r), 1 ^ y ^ z, r) ^ 1


def _encoded_kernel(y, r, family, naive=False, corrected=False):
    if naive:
        key = _b(_c(y, r, family), family)
        forward = lambda z: _c(_b(z, family), r, family)
    elif corrected:
        key = _b(_t(y, r, family, inverse=True), family)
        forward = lambda z: _t(_b(z, family), r, family)
    else:
        key = _quadratic_inverse(y, r)
        forward = lambda z: _hrt_unit(z, r) ^ 1
    result = {}
    for bits in range(1 << len(family)):
        toggle = sum(((bits >> i) & 1) << s for i, s in enumerate(family))
        sign = (-1) ** ((key & _mask(family) & (key ^ toggle)).bit_count() & 1)
        result[forward(key ^ toggle)] = sign
    return result  # common denominator sqrt(2**len(family))


def _clique(x, p):
    return sum(1 << ((1 << (i + 1)) | (1 << (j + 1)))
               for i, j in combinations(range(p), 2) if x >> i & 1 and x >> j & 1)


def _echo(y, x, p, linear=False):
    r, family = p + 1, _pairs(p + 1, first=1)
    def gauge(z):
        if linear:
            return z ^ _ring_product(z & _mask(family), z, r)
        return _c(z, r, family)
    return gauge(gauge(y) ^ _clique(x, p))


def _carrier_slots(r):
    return sum(1 << s for s in range(1 << r) if s & 1)


def _j(y, x, p):
    r = p + 1
    linear = sum(1 << (1 << (i + 1)) for i in range(p) if x >> i & 1)
    return y ^ _ring_product(linear, y & _carrier_slots(r), r)


def _middle(p):
    pairs = tuple(combinations(range(p), 2))
    q = len(pairs)
    h, gray = tuple(range(2 * p + q, 2 * p + 2 * q)), tuple(range(2 * p + 2 * q, 2 * p + 3 * q))
    outer_h = [('H', wire) for wire in h]
    u = outer_h + _gray_phase_batch([(i, j, h[k]) for k, (i, j) in enumerate(pairs)],
                                    gray, [1] * q, 3) + list(reversed(outer_h))
    c = _gray_phase_batch([(h[k], p + i, p + j) for k, (i, j) in enumerate(pairs)],
                          gray, [1] * q, 3)
    for k in range(q):
        c += [('H', 2 * p + k), ('CX', h[k], 2 * p + k), ('H', 2 * p + k)]
    return u + c + _inverse(u) + _inverse(c), 2 * p + 3 * q


def _middle_sign(basis, p):
    exponent = 0
    for k, (i, j) in enumerate(combinations(range(p), 2)):
        inverse_pair = ((basis >> (2 * p + k)) & 1) ^ (
            ((basis >> (p + i)) & (basis >> (p + j))) & 1)
        exponent ^= ((basis >> i) & (basis >> j) & 1) * inverse_pair
    return (-1) ** exponent


def _star(r, center, leaves):
    w = (1 << r) - 1
    family = tuple((1 << center) | (1 << j) for j in leaves)
    triples = [(s - 1, t - 1, (s | t) - 1)
               for s in family for t in range(1, 1 << r) if not s & t]
    targets = sorted({triple[2] for triple in triples})
    wall = [('H', wire) for wire in targets]
    phase, width = _dirty_cubic_phase([0] * w, [], triples)
    c = wall + phase + list(reversed(wall))
    pair_wall = [('H', s - 1) for s in family]
    b = pair_wall + _shared_ccz((1 << center) - 1,
                                [((1 << j) - 1, s - 1) for j, s in zip(leaves, family)])
    b += list(reversed(pair_wall))
    return c + b + pair_wall + b + c, width, family, len(triples), len(targets)


class PairGaugeFusionTests(unittest.TestCase):
    def test_complete_weak_echo_and_root_completion_with_arbitrary_dirty_bank(self):
        for p in (2, 3):
            r = p + 1
            inputs = range(128) if p == 2 else (0, 1, 2, 3, 17, 37, 255, 1023, 8191, 16383, 32767)
            family = _pairs(r, first=1)
            for bits in inputs:
                y = bits << 1
                s = y & _mask(family)
                for x in range(1 << p):
                    k = _clique(x, p)
                    actual = _echo(y, x, p)
                    expected = (_ring_product(_hrt_unit(k, r), y, r)
                                ^ _ring_product(_hrt_unit(s ^ k, r), k, r))
                    self.assertEqual(actual, expected)
                    self.assertEqual(_echo(actual, x, p), y)
                    carrier = y & _carrier_slots(r)
                    self.assertEqual(actual & _carrier_slots(r),
                                     _ring_product(_hrt_unit(k, r), carrier, r))
                    # Direct zeta product is independent of the central echo.
                    direct = carrier
                    for i in range(p):
                        if x >> i & 1:
                            direct = _ring_product(1 | (1 << (1 << (i + 1))), direct, r)
                    completed = _j(actual, x, p)
                    self.assertEqual(completed & _carrier_slots(r), direct)
                    self.assertEqual(_echo(_j(completed, x, p), x, p), y)

    def test_correct_partial_gauge_and_distinct_full_hadamard_failures(self):
        r, family = 3, _pairs(3)
        for bits in range(128):
            y = bits << 1
            self.assertEqual(_t(_t(y, r, family), r, family, inverse=True), y)
            self.assertEqual(_encoded_kernel(y, r, family, corrected=True),
                             _encoded_kernel(y, r, family))
        r = 4
        for family, expected_error in (((3, 12), Fraction(1, 2)), (_pairs(4), Fraction(7, 8))):
            exact = _encoded_kernel(0, r, family)
            naive = _encoded_kernel(0, r, family, naive=True)
            self.assertEqual(exact, _encoded_kernel(0, r, family, corrected=True))
            squared_numerator = sum((exact.get(row, 0) - naive.get(row, 0)) ** 2
                                    for row in exact.keys() | naive.keys())
            self.assertEqual(Fraction(squared_numerator, 1 << len(family)), expected_error)

    def test_intersecting_pair_fusion_includes_nonstar_triangles(self):
        for r, family, inputs in (
                (3, (3, 5, 6), range(128)),
                (4, (3, 5, 9), (0, 1, 7, 21, 42, 255, 1023, 8191, 32767)),
                (4, (3, 5, 6), (0, 1, 7, 21, 42, 255, 1023, 8191, 32767))):
            self.assertTrue(all(a & b for a, b in combinations(family, 2)))
            for bits in inputs:
                self.assertEqual(_encoded_kernel(bits << 1, r, family, naive=True),
                                 _encoded_kernel(bits << 1, r, family))

    def test_star_native_schedule_returns_arbitrary_helpers_and_literal_phases(self):
        r, leaves = 3, (1, 2)
        word, width, family, edges, targets = _star(r, 0, leaves)
        w = (1 << r) - 1
        self.assertEqual((edges, targets, width - w), (2, 1, 2))
        coefficient = ROOT_HALF * ROOT_HALF
        reverse = _inverse(word)
        for carrier, helpers in product((0, 1, 7, 21, 42, 64, 85, 127), range(4)):
            basis = carrier | (helpers << w)
            expected = {(row >> 1) | (helpers << w): coefficient * sign
                        for row, sign in _encoded_kernel(carrier << 1, r, family).items()}
            actual = _native(word, {basis: CO})
            self.assertEqual(actual, expected)
            self.assertEqual(_native(reverse, actual), {basis: CO})
        # A coherent reference can correlate different data/helper inputs.
        left, right = 21, 85 | (3 << w)
        state = {left: ROOT_HALF, right: CI * ROOT_HALF}
        actual = _native(word, state)
        self.assertEqual(_native(reverse, actual), state)
        for r in (3, 4, 6):
            leaves = tuple(range(1, r))
            word, width, family, edges, targets = _star(r, 0, leaves)
            k, w = len(leaves), (1 << r) - 1
            self.assertEqual(edges, k * ((1 << (r - 2)) - 1))
            self.assertEqual(targets, (1 << (r - 1)) - (1 << (r - 1 - k)) - k)
            self.assertEqual(width - w, edges)
            for native in (word, _inverse(word)):
                counts = Counter(gate[0] for gate in native)
                tcount = counts['T'] + counts['TDG']
                self.assertLessEqual(tcount, 16 * edges + 12 * k + 2 * (k & 1))
                self.assertLessEqual(len(native) - tcount, 16 * edges + 4 * targets + 25 * k + 6)
                self.assertLessEqual(len(_bands(native)), 24)
                self.assertTrue(all(len({gate[1] for gate in band}) == len(band)
                                    for band in _bands(native)))

    def test_paid_middle_all_small_inputs_and_actual_native_inverse(self):
        p = 2
        word, width = _middle(p)
        reverse = _inverse(word)
        for basis in range(1 << width):
            expected = {basis: CO * _middle_sign(basis, p)}
            actual = _native(word, {basis: CO})
            self.assertEqual(actual, expected)
            self.assertEqual(_native(reverse, actual), {basis: CO})
        # Shared addresses/singletons at p=3, with correlated nonzero helpers.
        p = 3
        word, width = _middle(p)
        principal = 2 * p + 3
        for low, helpers in product((0, 7, 63, 173, 511), (0, 1, 21, 63)):
            basis = low | (helpers << principal)
            self.assertEqual(_native(word, {basis: CO}),
                             {basis: CO * _middle_sign(basis, p)})

    def test_paid_middle_simultaneous_ledger_and_private_t_band_targets(self):
        for p in (2, 3, 5):
            q = p * (p - 1) // 2
            word, width = _middle(p)
            self.assertEqual(width - (2 * p + q), 2 * q)
            for native in (word, _inverse(word)):
                counts = Counter(gate[0] for gate in native)
                self.assertEqual(counts['T'] + counts['TDG'], 32 * q)
                self.assertEqual(len(native), 74 * q)
                self.assertEqual(len(_bands(native)), 32)
                for band in _bands(native):
                    self.assertEqual(len({gate[1] for gate in band}), len(band))
                    self.assertTrue(all(gate[1] >= 2 * p + 2 * q for gate in band))

    def test_linear_gauge_truncation_fails_the_actual_root_probe(self):
        p, r = 4, 5
        top = (1 << r) - 1
        for x, dirty in ((15, 0), (12, 1 << 6)):
            y = (1 << 1) | dirty  # arbitrary carrier constant equals one
            wrong = _j(_echo(y, x, p, linear=True), x, p)
            right = _j(_echo(y, x, p), x, p)
            self.assertNotEqual((wrong >> top) & 1, (right >> top) & 1)
            self.assertEqual(abs((-1) ** ((wrong >> top) & 1)
                                 - (-1) ** ((right >> top) & 1)), 2)
            self.assertEqual(_echo(_echo(y, x, p, linear=True), x, p, linear=True), y)

    def test_every_pair_hadamard_mask_fails_the_linear_gauge_pivot(self):
        def bilinear(s, c):
            return sum(((s >> i) & 1) * ((c >> (5 - i)) & 1) for i in range(6)) & 1
        def shear(index):
            s, c, z = index & 63, (index >> 6) & 63, index >> 12
            return s | ((c ^ s) << 6) | ((z ^ bilinear(s, c)) << 12)
        # Check this reduced sector against the full square-zero algebra map.
        family = _pairs(5, first=1)
        for s, c, z in product((0, 1, 3, 63), (0, 2, 21, 63), (0, 1)):
            y = (1 << 1) | (z << 31)
            y |= sum(((s >> i) & 1) << pair for i, pair in enumerate(family))
            y |= sum(((c >> i) & 1) << (pair | 1) for i, pair in enumerate(family))
            actual = y ^ _ring_product(y & _mask(family), y, 5)
            expected = shear(s | (c << 6) | (z << 12))
            self.assertEqual((actual >> 31) & 1, expected >> 12)
            self.assertEqual(sum(((actual >> (pair | 1)) & 1) << i
                                 for i, pair in enumerate(family)), (expected >> 6) & 63)
        for mask in range(64):
            def encoded(state):
                state = _permutation(state, shear)
                for bit in range(6):
                    if mask >> bit & 1:
                        state = _hadamard_numerator(state, bit)
                return _permutation(state, shear)
            def source(state):
                return encoded(_permutation(encoded(state), lambda index: index ^ 63))
            state = source({0: 1})
            state = {index: (-amplitude if index >> 12 else amplitude)
                     for index, amplitude in state.items()}
            state = source(state)
            output = mask | (mask << 6) | ((mask.bit_count() & 1) << 12)
            self.assertEqual(state, {output: 1 << (2 * mask.bit_count())})
            # The required full-clique root phase is -1 on the original input.
            self.assertNotEqual(state, {0: -(1 << (2 * mask.bit_count()))})


if __name__ == '__main__':
    unittest.main()
