"""Exact bounded checks for Clifford restrictions on a one-band clean code.

The arbitrary-width classification is analytic. These fixtures retain
literal phases, actual inverse words, and cheaper complete baselines for
the two failed extension mechanisms. No dirty-return query is inferred
from an arbitrary-data transition.
"""
from collections import Counter, deque
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CQ2, Q2, CO, CZ, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import (
    _clifford_columns, _conjugate, _tau,
)


def _key(matrix):
    return tuple(tuple(row) for row in matrix)


def _one_qubit_cliffords():
    generators = [native_columns(1, [(name, 0)]) for name in ('H', 'S')]
    initial = eye(2)
    found = {_key(initial): initial}
    pending = deque([initial])
    while pending:
        current = pending.popleft()
        for generator in generators:
            candidate = multiply(generator, current)
            key = _key(candidate)
            if key not in found:
                found[key] = candidate
                pending.append(candidate)
    return list(found.values())


def _trace(matrix):
    return sum((matrix[j][j] for j in range(len(matrix))), CZ) / len(matrix)


def _pauli_label(matrix, width):
    """Recognize a signed Hermitian Pauli without enumerating the group."""
    size = 1 << width
    first = [row for row in range(size) if matrix[row][0] != CZ]
    if len(first) != 1:
        return None
    x = first[0]
    phase = matrix[x][0]
    z = 0
    for bit in range(width):
        column = 1 << bit
        value = matrix[column ^ x][column]
        if value == -phase:
            z |= column
        elif value != phase:
            return None
    label = x | (z << width)
    candidate = _tau(label, width)
    return label if matrix == candidate or matrix == scale(candidate, -1) else None


def _rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def _high_rank_words(pairs):
    # Only wire c_1=0 is initialized. All other wires are arbitrary data.
    link = [('CX', 2 * j, 0) for j in range(1, pairs)]
    local_k = [gate for j in range(pairs)
               for gate in [('H', 2 * j), ('CX', 2 * j, 2 * j + 1)]]
    coordinates = [gate for j in range(pairs)
                   for gate in [('CX', 2 * j, 2 * j + 1),
                                ('H', 2 * j),
                                ('CX', 2 * j, 2 * j + 1)]]
    e = [('T', wire) for wire in range(2 * pairs)]
    middle = [gate for j in range(pairs)
              for gate in [('T', 2 * j), ('TDG', 2 * j + 1)]]
    b = coordinates + middle + _inverse(coordinates)
    k = link + local_k
    ua = k + e
    ub = k + e + _inverse(b)
    displayed = k + e + b + _inverse(e) + _inverse(k)
    pair_cx = [('CX', 2 * j, 2 * j + 1) for j in range(pairs)]
    baseline = link + pair_cx + middle + _inverse(pair_cx) + _inverse(link)
    return ua, ub, displayed, baseline


def _nonclifford_extension_words():
    # D=(S_t H_t) CX(t,c), E0=D (T_c† T_t) D†.
    # All S inverses are explicitly charged as three S gates.
    d = [('CX', 1, 0), ('H', 1), ('S', 1)]
    d_inverse = [('S', 1)] * 3 + [('H', 1), ('CX', 1, 0)]
    e0 = d_inverse + [('TDG', 0), ('T', 1)] + d
    e0_inverse = d_inverse + [('T', 0), ('TDG', 1)] + d
    ua = [('H', 0), ('T', 0), ('T', 1)]
    ua_inverse = [('TDG', 1), ('TDG', 0), ('H', 0)]
    ub = e0 + ua
    displayed = ua + ua_inverse + e0_inverse
    displayed_inverse = e0 + ua + ua_inverse
    return ua, ub, displayed, displayed_inverse, e0_inverse, e0


class CleanCodeCliffordRigidityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.one_qubit_group = _one_qubit_cliffords()
        t = native_columns(1, [('T', 0)])
        cls.magic = _conjugate(t, _tau(1, 1))
        cls.high_words = _high_rank_words(2)
        cls.high_ua = native_columns(4, cls.high_words[0])
        cls.high_ub = native_columns(4, cls.high_words[1])
        cls.high_s = _conjugate(cls.high_ua, _tau(16, 4))

    def test_one_qubit_literal_stabilizer_has_exactly_two_elements(self):
        self.assertEqual(len(self.one_qubit_group), 192)
        twice_code = add(eye(2), self.magic)
        fixes = {_key(w) for w in self.one_qubit_group
                 if multiply(w, twice_code) == twice_code}
        self.assertEqual(fixes, {_key(eye(2)), _key(self.magic)})
        # A common phase is retained; multiplying a pointwise stabilizer
        # by a nontrivial eighth root does not remain pointwise identity.
        omega = CQ2(Q2(0, Fraction(1, 2)), Q2(0, Fraction(1, 2)))
        self.assertNotEqual(multiply(scale(self.magic, omega), twice_code),
                            twice_code)

    def test_literal_trace_gap_on_all_one_qubit_cliffords(self):
        literal_squared_floor = Q2(Fraction(3, 2), -1)
        for w in self.one_qubit_group:
            overlap = _trace(w) + _trace(multiply(w, self.magic))
            distance_squared = 2 - 2 * overlap.real
            if distance_squared != 0:
                self.assertGreaterEqual(distance_squared, literal_squared_floor)
        # These bounded examples check exact phase arithmetic. The universal
        # gap follows from the analytic Clifford trace/support argument.

    def test_rank_two_exception_is_a_literal_klein_four_group(self):
        e = native_columns(2, [('T', 0), ('T', 1)])
        s = _conjugate(e, _tau(3, 2))
        swap = native_columns(2, [('CX', 0, 1), ('CX', 1, 0), ('CX', 0, 1)])
        other = multiply(s, swap)
        group = [eye(4), s, swap, other]
        keys = {_key(w) for w in group}
        self.assertEqual(len(keys), 4)
        twice_code = add(eye(4), s)
        for w in group:
            self.assertEqual(multiply(w, twice_code), twice_code)
            self.assertEqual(multiply(w, w), eye(4))
            self.assertEqual(adjoint(w), w)
            _clifford_columns(w, 2)
            for v in group:
                self.assertIn(_key(multiply(w, v)), keys)
        self.assertEqual(add(swap, other), twice_code)
        self.assertEqual(swap, scale(add(add(add(eye(4), _tau(3, 2)),
                                             _tau(15, 2)), _tau(12, 2)), Fraction(1, 2)))
        self.assertEqual(other, scale(add(add(add(eye(4), _tau(11, 2)),
                                             _tau(7, 2)), scale(_tau(12, 2), -1)),
                                      Fraction(1, 2)))

    def test_high_rank_single_clean_image_and_all_logical_images_match(self):
        for bit in range(4):
            pauli = _tau(1 << (4 + bit), 4)
            self.assertEqual(_conjugate(self.high_ua, pauli),
                             _conjugate(self.high_ub, pauli))
        columns = _clifford_columns(self.high_s, 4)
        self.assertEqual(_rank([column ^ (1 << bit)
                                for bit, column in enumerate(columns)]), 4)
        for j in range(2):
            x = (1 << (2 * j)) | (1 << (2 * j + 1))
            self.assertEqual(_conjugate(self.high_s, _tau(x, 4)),
                             _tau(x | (x << 4), 4))
        # The rank-two SWAP exception does not persist after adding a
        # second active pair under the same single clean parity condition.
        swap = native_columns(4, [('CX', 0, 1), ('CX', 1, 0), ('CX', 0, 1)])
        twice_code = add(eye(16), self.high_s)
        self.assertNotEqual(multiply(swap, twice_code), twice_code)

    def test_high_rank_split_has_no_common_pauli_clean_flip(self):
        common = []
        for z in range(16):
            pauli = _tau(z << 4, 4)
            first = _conjugate(adjoint(self.high_ua), pauli)
            second = _conjugate(adjoint(self.high_ub), pauli)
            self.assertIsNotNone(_pauli_label(first, 4))
            expected = all(((z >> (2 * j)) & 3).bit_count() % 2 == 0
                           for j in range(2))
            self.assertEqual(_pauli_label(second, 4) is not None, expected)
            if expected:
                common.append(z)
                self.assertEqual(_conjugate(pauli, self.high_s), self.high_s)
        self.assertEqual(common, [0, 3, 12, 15])
        for bit in range(4):
            self.assertIsNone(_pauli_label(
                _conjugate(adjoint(self.high_ua), _tau(1 << bit, 4)), 4))
        # Analytically, any physical X/Y component expands through E into
        # distinct native-axis terms. Thus all possible common flips were
        # among the Z-type candidates checked above, not just generators.

    def test_high_rank_transition_full_baseline_and_actual_inverse(self):
        displayed, baseline = self.high_words[2:]
        actual = native_columns(4, displayed)
        self.assertEqual(actual, native_columns(4, baseline))
        self.assertEqual(native_columns(4, _inverse(displayed)), adjoint(actual))
        self.assertEqual(native_columns(4, _inverse(baseline)), adjoint(actual))
        for row in range(16):
            for column in range(16):
                if row != column:
                    self.assertEqual(actual[row][column], CZ)
        self.assertNotEqual(actual, adjoint(actual))
        # Every physical column is checked. Only c_1=0 is initialized;
        # all other wires are arbitrary data, not identity-return helpers.

    def test_high_rank_family_parallel_bands_and_fair_ledgers(self):
        for pairs in (1, 2, 5):
            _, _, displayed, baseline = _high_rank_words(pairs)
            for word in (displayed, _inverse(displayed)):
                counts = Counter(gate[0] for gate in word)
                self.assertEqual(counts['T'] + counts['TDG'], 6 * pairs)
                self.assertEqual(counts['CX'] + counts['H'], 12 * pairs - 2)
                self.assertEqual([len(b) for b in _native_t_bands(word)],
                                 [2 * pairs] * 3)
            for word in (baseline, _inverse(baseline)):
                counts = Counter(gate[0] for gate in word)
                self.assertEqual(counts['T'] + counts['TDG'], 2 * pairs)
                self.assertEqual(counts['CX'], 4 * pairs - 2)
                self.assertEqual([len(b) for b in _native_t_bands(word)], [2 * pairs])
            for word in (displayed, baseline):
                self.assertEqual({wire for gate in word for wire in gate[1:]},
                                 set(range(2 * pairs)))
                for band in _native_t_bands(word):
                    self.assertEqual({gate[1] for gate in band}, set(range(2 * pairs)))

    def test_correlated_two_band_image_need_not_be_clifford(self):
        ua_word, ub_word, *_ = _nonclifford_extension_words()
        ua, ub = native_columns(2, ua_word), native_columns(2, ub_word)
        clean, target = _tau(4, 2), _tau(8, 2)
        s = _conjugate(ua, clean)
        a, c = _conjugate(ub, target), _conjugate(ua, target)
        self.assertEqual(s, _conjugate(ub, clean))
        self.assertEqual(c, target)
        self.assertEqual(multiply(a, s), multiply(s, a))
        twice_code = add(eye(4), s)
        self.assertEqual(multiply(a, twice_code), multiply(c, twice_code))
        # R=Z_t is an actual probe commuting with the outer T axes and S.
        # Its nonzero character produces a genuine non-Clifford A:
        # A R A=(C R C)S, with L_S contained in the outer axis span.
        self.assertEqual(_conjugate(a, target), multiply(target, s))
        self.assertIsNone(_pauli_label(multiply(target, s), 2))
        with self.assertRaises(ValueError):
            _clifford_columns(a, 2)
        self.assertEqual([len(b) for b in _native_t_bands(ua_word)], [2])
        self.assertEqual([len(b) for b in _native_t_bands(ub_word)], [2, 2])

    def test_nonclifford_extension_has_a_paid_one_band_identity_code_baseline(self):
        _, _, displayed, inverse, baseline, baseline_inverse = _nonclifford_extension_words()
        actual = native_columns(2, displayed)
        self.assertEqual(actual, native_columns(2, baseline))
        self.assertEqual(native_columns(2, inverse), adjoint(actual))
        self.assertEqual(native_columns(2, baseline_inverse), adjoint(actual))
        for column in (0, 2):
            for row in range(4):
                self.assertEqual(actual[row][column], CO if row == column else CZ)
        self.assertNotEqual(actual, eye(4))
        self.assertNotEqual(actual, adjoint(actual))
        for word, t_count, cliffords in ((displayed, 6, 10), (inverse, 6, 10),
                                         (baseline, 2, 8), (baseline_inverse, 2, 8)):
            counts = Counter(gate[0] for gate in word)
            self.assertEqual(counts['T'] + counts['TDG'], t_count)
            self.assertEqual(counts['H'] + counts['S'] + counts['CX'], cliffords)
        # The displayed composition has three scheduled two-wire bands;
        # two adjacent inverse bands cancel. The fairly simplified target
        # is the actual E0† word, with one band and its full unused block.
        self.assertEqual([len(b) for b in _native_t_bands(baseline)], [2])
        self.assertEqual([len(b) for b in _native_t_bands(baseline_inverse)], [2])


if __name__ == '__main__':
    unittest.main()
