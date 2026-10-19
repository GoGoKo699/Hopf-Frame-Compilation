"""Exact module-barrier evidence and a paid two-representation escape.

The lower bound is proved algebraically. These fixtures retain symplectic
coordinate changes, literal commutator phases, and complete native inputs.
"""
from fractions import Fraction
from random import Random
import unittest

from tests.exact_matrices import CZ, CO, CI, adjoint, native_columns, scale
from tests.test_balanced_algebra_permutation import _inverse, _ring_product
from tests.test_returned_diagonal_echo_degree import _ccz
from tests.test_two_layer_diagonal_rigidity import _primitive_inverse, _reduce_adjacent_s
from tests.test_two_layer_galois_return import (
    _tau, _product, _conjugate, _coefficients, _clifford_columns, _symplectic,
)


def _rank(columns):
    pivots = {}
    for column in columns:
        while column:
            bit = column.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = column
                break
            column ^= pivots[bit]
    return len(pivots)


def _apply(columns, vector):
    result = 0
    for bit, column in enumerate(columns):
        if vector >> bit & 1:
            result ^= column
    return result


def _symplectic_word(label, width, gates):
    for gate in gates:
        if gate[0] == 'H':
            bit = gate[1]
            if ((label >> bit) ^ (label >> (width + bit))) & 1:
                label ^= (1 << bit) | (1 << (width + bit))
        else:
            control, target = gate[1:]
            label ^= ((label >> control) & 1) << target
            label ^= ((label >> (width + target)) & 1) << (width + control)
    return label


def _regular_representation(d, mixed=False):
    size, top = 1 << d, (1 << d) - 1
    width = size // 2
    subsets = [s for s in range(size) if not s & 1]
    labels = {s: 1 << j for j, s in enumerate(subsets)}
    labels.update({top ^ s: 1 << (width + j) for j, s in enumerate(subsets)})
    inverse_labels = {v.bit_length() - 1: s for s, v in labels.items()}
    def encode(value):
        return sum(labels[s] for s in range(size) if value >> s & 1)
    def decode(value):
        return sum(1 << inverse_labels[bit] for bit in range(size) if value >> bit & 1)
    rng = Random(20261019 + d)
    gates = []
    if mixed:
        for _ in range(4 * size):
            if width == 1 or rng.randrange(2):
                gates.append(('H', rng.randrange(width)))
            else:
                control, target = rng.sample(range(width), 2)
                gates.append(('CX', control, target))
    forward = lambda x: _symplectic_word(x, width, gates)
    reverse = lambda x: _symplectic_word(x, width, list(reversed(gates)))
    matrices = []
    for bit in range(d):
        columns = []
        for col in range(size):
            original = decode(reverse(1 << col))
            image = _ring_product(1 << (1 << bit), original)
            columns.append(forward(encode(image)))
        matrices.append(columns)
    unit = forward(encode(1))
    return width, matrices, unit


def _fredkin(control, first, second):
    ccx = ([('H', first)]
           + [g for _, gates in _ccz(control, second, first) for g in gates]
           + [('H', first)])
    return [('CX', first, second)] + ccx + [('CX', first, second)]


def _split_native_query(table):
    # Original addresses 0,1; arbitrary output 2; arbitrary banks 3,4 and 5,6.
    left = _fredkin(0, 3, 4)
    right = _fredkin(1, 5, 6)
    fixed = [('CX', 3 + i, 5 + j) for i in range(2) for j in range(2)
             if table >> (2 * i + j) & 1]
    a = _inverse(left) + [('CX', 2, 3)] + left
    b = _inverse(fixed) + _inverse(right) + [('S', 5)] * 2 + right + fixed
    inverse_b = _reduce_adjacent_s(_primitive_inverse(b))
    word = [('H', 2)] + a + b + _inverse(a) + inverse_b + [('H', 2)]
    return word


def _actual_inverse(word):
    return _reduce_adjacent_s(_primitive_inverse(word))


def _controlled_pauli(control, label, carrier):
    """Positive Hermitian tensor Pauli, with literal positive-S words."""
    word = []
    width = len(carrier)
    for bit, target in enumerate(carrier):
        x, z = (label >> bit) & 1, (label >> (width + bit)) & 1
        if x and z:
            word += [('S', target)] * 3 + [('CX', control, target), ('S', target)]
        elif x:
            word.append(('CX', control, target))
        elif z:
            word += [('H', target), ('CX', control, target), ('H', target)]
    return word


def _native_probe_fixture(p):
    if p == 1:
        # At address 0 use X, at address 1 use Y; retain omega^address.
        carrier, gamma, dual = [1], [1, 3], [3, 1]
        source = ([('S', 0)] + _controlled_pauli(0, 2, carrier)
                  + [('H', 1), ('S', 1), ('S', 1), ('H', 1), ('T', 0)])
    else:
        # At address b0+2b1 use X/Z on carrier b1. Both Fredkins are paid.
        carrier, gamma, dual = [2, 3], [1, 4, 2, 8], [4, 1, 8, 2]
        swap = _fredkin(1, 2, 3)
        base = ([('S', 0)] * 3 + _controlled_pauli(0, 3, [2])
                + [('H', 2), ('S', 2), ('S', 2), ('H', 2)])
        source = swap + base + _actual_inverse(swap) + [('T', 0)]
    outputs = list(range(p + len(carrier), p + len(carrier) + len(gamma)))
    wrapper = [gate for bit, label in zip(outputs, dual)
               for gate in _controlled_pauli(bit, label, carrier)]
    hadamards = [('H', bit) for bit in outputs]
    word = (hadamards + wrapper + source + _actual_inverse(wrapper)
            + _actual_inverse(source) + hadamards)
    return carrier, gamma, dual, source, wrapper, word


def _ledger(word):
    """Count the emitted word and certify its consecutive native T bands."""
    t_count = g_count = depth = 0
    targets = set()
    for gate in word:
        if gate[0] in ('T', 'TDG'):
            t_count += 1
            if not targets or gate[1] in targets:
                depth += 1
                targets = set()
            targets.add(gate[1])
        else:
            g_count += 1
            targets = set()
    return t_count, g_count, depth


class SymplecticModuleQueryBarrierTests(unittest.TestCase):
    def test_free_cyclic_module_and_sharp_rank_in_mixed_coordinates(self):
        for d in range(1, 7):
            for mixed in (False, True):
                width, deviations, unit = _regular_representation(d, mixed)
                size = 2 * width
                for left in deviations:
                    self.assertEqual(_rank(left), 1 << (d - 1))
                    self.assertTrue(all(_apply(left, col) == 0 for col in left))
                    for right in deviations:
                        self.assertEqual([_apply(left, col) for col in right],
                                         [_apply(right, col) for col in left])
                vectors = []
                for subset in range(1 << d):
                    value = unit
                    for bit in range(d):
                        if subset >> bit & 1:
                            value = _apply(deviations[bit], value)
                    vectors.append(value)
                self.assertEqual(_rank(vectors), 1 << d)
                self.assertEqual(_symplectic(unit, vectors[-1], width), 1)
                for address in range(1 << d):
                    value = unit
                    for bit in range(d):
                        if address >> bit & 1:
                            value ^= _apply(deviations[bit], value)
                    self.assertEqual(_symplectic(unit, value ^ unit, width),
                                     int(address == (1 << d) - 1))
                # Every I+N_i is symplectic, including after basis mixing.
                for columns in deviations:
                    image = [column ^ (1 << i) for i, column in enumerate(columns)]
                    for i in range(size):
                        for j in range(size):
                            self.assertEqual(_symplectic(image[i], image[j], width),
                                             _symplectic(1 << i, 1 << j, width))

    def test_clifford_support_and_nonreal_commutator_phase(self):
        fixtures = [[], [('H', 0)], [('S', 0)], [('CX', 0, 1)],
                    [('H', 0), ('S', 1), ('CX', 0, 1), ('H', 1)],
                    [('S', 0), ('CX', 1, 0), ('H', 0), ('S', 1)]]
        for word in fixtures:
            operator = native_columns(2, word)
            columns = _clifford_columns(operator, 2)
            rank = _rank([column ^ (1 << bit) for bit, column in enumerate(columns)])
            coefficients = _coefficients(operator, 2)
            self.assertEqual(sum(value != CZ for value in coefficients), 1 << rank)
        # C=S, v=u=X has K=Y X=-i Z. The actual adjoint must cancel
        # this nonreal scalar before the final commutator is -I.
        c = native_columns(1, [('S', 0)])
        p = _tau(1, 1)
        k = _product(c, p, adjoint(c), adjoint(p))
        self.assertEqual(k, scale(_tau(2, 1), -CI))
        self.assertEqual(_product(k, p, adjoint(k), adjoint(p)),
                         scale(native_columns(1, []), -1))

    def test_native_source_support_counts_conjugation_once(self):
        source = ([('H', 2)]
                  + [g for _, gates in _ccz(0, 1, 2) for g in gates]
                  + [('H', 2)])
        f = native_columns(3, source)
        image = _conjugate(f, _tau(1, 3))
        columns = _clifford_columns(image, 3)
        rank = _rank([col ^ (1 << i) for i, col in enumerate(columns)])
        self.assertEqual(rank, 2)
        self.assertEqual(sum(c != CZ for c in _coefficients(image, 3)), 1 << rank)
        t_count = sum(g[0] in ('T', 'TDG') for g in source)
        self.assertEqual(t_count, 7)
        self.assertLessEqual(rank, t_count)

    def test_two_half_native_query_full_operator_and_actual_inverse(self):
        for table in (1, 6, 8, 15):
            word = _split_native_query(table)
            inverse = _reduce_adjacent_s(_primitive_inverse(word))
            actual = native_columns(7, word)
            for column in range(128):
                row_index = 2 * (column & 1) + ((column >> 1) & 1)
                output = column ^ (((table >> row_index) & 1) << 2)
                for row in range(128):
                    self.assertEqual(actual[row][column], CO if row == output else CZ)
            self.assertEqual(native_columns(7, inverse), adjoint(actual))
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 56)
            self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in word),
                             120 + 4 * table.bit_count())
            self.assertLessEqual(len(inverse), len(word))

    def test_two_families_escape_the_single_representation_width_bound(self):
        # Six-address AND: two eight-bit banks use w=16, below 2^(6-1).
        half, m = 3, 8
        width = 2 * m
        self.assertLess(width, 1 << (2 * half - 1))
        for u in range(m):
            gamma = 1 << u
            for v in range(m):
                delta = (1 << (width + m + v))
                if v == m - 1:
                    delta ^= 1 << (width + m - 1)
                self.assertEqual(_symplectic(gamma, delta, width),
                                 int(u == m - 1 and v == m - 1))
        def translate(label, bank, shift):
            result = 0
            for phase in (0, 1):
                for wire in range(width):
                    destination = wire
                    if bank * m <= wire < (bank + 1) * m:
                        destination = bank * m + ((wire - bank * m) ^ shift)
                    result ^= ((label >> (phase * width + wire)) & 1) << (phase * width + destination)
            return result
        def fixed(label):
            return _symplectic_word(label, width, [('CX', m - 1, 2 * m - 1)])
        left = lambda label: translate(label, 0, 1)
        right = lambda label: fixed(translate(fixed(label), 1, 1))
        probe = 1 << (width + 2 * m - 1)
        self.assertNotEqual(left(right(probe)), right(left(probe)))

    def _check_native_probe_indicator(self, p, expected_source, expected_indicator):
        carrier, gamma, dual, source, wrapper, word = _native_probe_fixture(p)
        width, m = len(carrier), len(gamma)
        self.assertEqual(_rank(gamma), m)
        for j, q in enumerate(dual):
            for i, g in enumerate(gamma):
                self.assertEqual(_symplectic(q, g, width), int(i == j))
        self.assertTrue(any(_symplectic(q, r, width) for q in dual for r in dual))
        omega = native_columns(1, [('T', 0)])[1][1]
        actual_source = native_columns(p + width, source)
        for column in range(1 << (p + width)):
            address, data = column & (m - 1), column >> p
            phase = omega if address & 1 else CO
            pauli = _tau(gamma[address], width)
            for row in range(1 << (p + width)):
                expected = (phase * pauli[row >> p][data]
                            if row & (m - 1) == address else CZ)
                self.assertEqual(actual_source[row][column], expected)
        inverse_source = _actual_inverse(source)
        self.assertEqual(native_columns(p + width, inverse_source), adjoint(actual_source))
        self.assertEqual((_ledger(source), _ledger(inverse_source)), expected_source)
        self.assertLessEqual(_ledger(wrapper)[1], 5 * m * width)
        self.assertLessEqual(_ledger(_actual_inverse(wrapper))[1], 5 * m * width)

        total_width = p + width + m
        actual = native_columns(total_width, word)
        for column in range(1 << total_width):
            address = column & (m - 1)
            output = column ^ (1 << (p + width + address))
            for row in range(1 << total_width):
                self.assertEqual(actual[row][column], CO if row == output else CZ)
        inverse = _actual_inverse(word)
        self.assertEqual(native_columns(total_width, inverse), adjoint(actual))
        self.assertEqual((_ledger(word), _ledger(inverse)), expected_indicator)
        source_t = max(_ledger(source)[0], _ledger(inverse_source)[0])
        source_g = max(_ledger(source)[1], _ledger(inverse_source)[1])
        source_d = max(_ledger(source)[2], _ledger(inverse_source)[2])
        self.assertLessEqual(_ledger(word)[0], 2 * source_t)
        self.assertLessEqual(_ledger(word)[1], 2 * source_g + 10 * m * width + 2 * m)
        self.assertLessEqual(_ledger(word)[2], 2 * source_d)

    def test_native_noncommuting_probe_indicator_retains_address_scalar(self):
        self._check_native_probe_indicator(
            1, ((1, 8, 1), (1, 10, 1)), ((2, 34, 2), (2, 34, 2)))

    def test_native_four_probe_indicator_returns_every_carrier_column(self):
        self._check_native_probe_indicator(
            2, ((15, 40, 11), (15, 38, 11)), ((30, 102, 22), (30, 102, 22)))

    def test_quadratic_indicator_cliffords_fit_every_r91_allocation_branch(self):
        branches = set()
        for a in range(1, 13):
            q_rows = 1 << a
            for q in list(range(1, 65)) + [q_rows - 1, q_rows, q_rows + 1, 3 * q_rows]:
                for t_exp in range(a + 8):
                    t = 1 << t_exp
                    if q >= q_rows and q_rows <= t:
                        branches.add('direct')
                        self.assertLessEqual(q_rows * q_rows, q_rows * q)
                        continue
                    j = min(q_rows, t)
                    while j * j > q_rows * q:
                        j //= 2
                    g = min(q, t)
                    h_limit = int(min(Fraction(q_rows, j), Fraction(t, g)))
                    h = 1 << (h_limit.bit_length() - 1)
                    k = q_rows // (h * j)
                    self.assertEqual(k * h * j, q_rows)
                    self.assertLessEqual(g * h, t)
                    self.assertLessEqual(h * h + j * j, 5 * q_rows * q)
                    self.assertLessEqual(h + j, 2 * t)
                    if t * t >= q_rows * q:
                        branches.add('wide')
                        self.assertLessEqual(k * ((q + g - 1) // g), 4)
                    elif q <= t:
                        branches.add('bounded_outputs')
                        self.assertLessEqual(k * ((q + g - 1) // g) * t * t,
                                             2 * q_rows * q)
                    else:
                        branches.add('many_outputs')
                        self.assertEqual(h, 1)
                        self.assertLessEqual(k * ((q + g - 1) // g) * t * t,
                                             2 * q_rows * q)
        self.assertEqual(branches, {'direct', 'wide', 'bounded_outputs', 'many_outputs'})

    def test_formal_walsh_indicator_return_inverse_and_non_native_phase_gap(self):
        # R denotes a formal phase of +/- pi/m, not a native T gate.
        def run(word, value, modulus):
            exponent = 0
            for gate, first, second in word:
                if gate == 'R':
                    exponent += second * ((value >> first) & 1)
                else:
                    value ^= ((value >> first) & 1) << second
            return value, exponent % modulus

        for p in range(1, 4):
            m = 1 << p
            gray = [a ^ (a >> 1) for a in range(m)]
            word, parity_labels = [], set()
            for index, a in enumerate(gray):
                word += [('R', p + u, -1 if (a & u).bit_count() & 1 else 1)
                         for u in range(m)]
                parity_labels.update((1 << (p + u)) | a for u in range(m))
                edge = a ^ gray[(index + 1) % m]
                self.assertEqual(edge.bit_count(), 1)
                word += [('CX', edge.bit_length() - 1, p + u) for u in range(m)]
            inverse = [(gate, first, -second if gate == 'R' else second)
                       for gate, first, second in reversed(word)]
            self.assertEqual(len(parity_labels), m * m)
            for orientation, sign in ((word, 1), (inverse, -1)):
                self.assertEqual(sum(g[0] == 'R' for g in orientation), m * m)
                self.assertEqual(sum(g[0] == 'CX' for g in orientation), m * m)
                bands = [orientation[start:start + m]
                         for start in range(0 if sign == 1 else m, len(word), 2 * m)]
                self.assertEqual(len(bands), m)
                self.assertTrue(all(len({g[1] for g in band}) == m
                                    and all(g[0] == 'R' for g in band) for band in bands))
                for column in range(1 << (p + m)):
                    address = column & (m - 1)
                    expected_phase = sign * m * ((column >> (p + address)) & 1)
                    self.assertEqual(run(orientation, column, 2 * m),
                                     (column, expected_phase % (2 * m)))

        # Every native affine parity has vanishing fourth difference mod 8.
        # The four-way Boolean sign (the p=3 restricted indicator) does not.
        for mask in range(16):
            for bias in range(2):
                difference = sum((-1) ** (4 - x.bit_count())
                                 * (((mask & x).bit_count() ^ bias) & 1)
                                 for x in range(16))
                self.assertEqual(difference % 8, 0)
        sign_difference = sum((-1) ** (4 - x.bit_count()) * (4 if x == 15 else 0)
                              for x in range(16))
        self.assertEqual(sign_difference % 8, 4)


if __name__ == '__main__':
    unittest.main()
