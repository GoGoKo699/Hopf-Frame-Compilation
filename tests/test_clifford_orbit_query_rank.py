"""Exact finite receipts for the phase lift and multiplicative rank argument.

The all-width theorem is analytic. These fixtures retain all normal-Pauli
phases and use actual query permutations, not affine symplectic data alone.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, Q2, adjoint, eye, multiply, native_columns, scale,
)


def _xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def _bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def _anf(values):
    coefficients = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for mask in range(len(values)):
            if mask >> bit & 1:
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return coefficients


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


@lru_cache(None)
def _feature_masks(width):
    z_size = 2 * width + 1
    return ((0,) + tuple(1 << i for i in range(z_size))
            + tuple((1 << i) | (1 << j)
                    for i, j in combinations(range(z_size), 2))
            + (1 << z_size,))


def _features(width, state):
    return sum(int(state & mask == mask) << i
               for i, mask in enumerate(_feature_masks(width)))


def _lift(width, action):
    coefficients = _anf([_features(width, action(state))
                         for state in range(1 << (2 * width + 2))])
    masks = _feature_masks(width)
    assert all(not value or mask in masks
               for mask, value in enumerate(coefficients))
    return tuple(coefficients[mask] for mask in masks)


def _apply(columns, vector):
    return _xor(columns[i] for i in _bits(vector))


def _compose(left, right):
    return tuple(_apply(left, column) for column in right)


def _pack(columns):
    return sum(column << (len(columns) * i)
               for i, column in enumerate(columns))


def _degree(values):
    return max((mask.bit_count() for mask, value in enumerate(_anf(values))
                if value), default=0)


@lru_cache(None)
def _normal(width, state):
    mask = (1 << width) - 1
    a, b, phase = state & mask, (state >> width) & mask, state >> (2 * width)
    root = (CO, CI, -CO, -CI)[phase]
    return tuple(tuple(root * (-1 if (column & b).bit_count() & 1 else 1)
                       if row == column ^ a else CZ
                       for column in range(1 << width))
                 for row in range(1 << width))


def _key(matrix):
    return tuple(tuple(row) for row in matrix)


def _tableau(width, unitary):
    normal = {_normal(width, state): state
              for state in range(1 << (2 * width + 2))}
    return tuple(normal[_key(multiply(multiply(unitary, _normal(width, 1 << j)),
                                     adjoint(unitary)))]
                 for j in range(2 * width))


def _normal_action(width, tableau, state):
    mask = (1 << width) - 1
    label, phase = 0, state >> (2 * width)
    for j in _bits(state & ((1 << (2 * width)) - 1)):
        image = tableau[j]
        phase += (image >> (2 * width)) + 2 * (
            ((label >> width) & (image & mask)).bit_count() & 1)
        label ^= image & ((1 << (2 * width)) - 1)
    return label | ((phase % 4) << (2 * width))


def _unitary_lift(width, unitary):
    tableau = _tableau(width, unitary)
    return _lift(width, lambda state: _normal_action(width, tableau, state))


def _scalar(matrix):
    return all(entry == (matrix[0][0] if i == j else CZ)
               for i, row in enumerate(matrix) for j, entry in enumerate(row))


@lru_cache(None)
def _one_qubit_cliffords():
    gates = [native_columns(1, [(name, 0)]) for name in ('H', 'S')]
    lifts = [_unitary_lift(1, gate) for gate in gates]
    identity = _unitary_lift(1, eye(2))
    found, queue = {identity: eye(2)}, [identity]
    for current in queue:
        for gate, lift in zip(gates, lifts):
            candidate = _compose(lift, current)
            if candidate not in found:
                found[candidate] = multiply(gate, found[current])
                queue.append(candidate)
    return found


def _shear_basis(x, basis):
    """Two source wires, one target; M=(x0,x1), offset=x2."""
    shift = (((x & 1) & (basis & 1))
             ^ (((x >> 1) & 1) & ((basis >> 1) & 1))
             ^ ((x >> 2) & 1))
    return basis ^ (shift << 2)


def _shear_action(x, state):
    a, b, phase = state & 7, (state >> 3) & 7, state >> 6
    av = ((x & 1) & (a & 1)) ^ (((x >> 1) & 1) & ((a >> 1) & 1))
    bv = (b >> 2) & 1
    a ^= av << 2
    b ^= (x & 3) if bv else 0
    phase = (phase + 2 * bv * ((x >> 2) & 1)) % 4
    return a | (b << 3) | (phase << 6)


def _algebra_multiply(left, right):
    """Coefficient bitsets of the free commutative square-zero algebra."""
    return _xor(1 << (i | j) for i in _bits(left) for j in _bits(right)
                if not i & j)


def _group(generators):
    values = {1}
    for generator in generators:
        values |= {_algebra_multiply(value, generator) for value in values}
    return values


class CliffordOrbitQueryRankTests(unittest.TestCase):
    def test_phase_lift_matches_every_small_normal_pauli_and_retains_carry(self):
        words = [[('H', 0)], [('S', 0)], [('CX', 0, 1)],
                 [('H', 0), ('S', 1), ('CX', 0, 1), ('H', 1), ('S', 0)]]
        width = 2
        for word in words:
            unitary = native_columns(width, word)
            tableau = _tableau(width, unitary)
            lift = _unitary_lift(width, unitary)
            for state in range(1 << (2 * width + 2)):
                image = _normal_action(width, tableau, state)
                actual = multiply(multiply(unitary, _normal(width, state)),
                                  adjoint(unitary))
                self.assertEqual(_key(actual), _normal(width, image))
                self.assertEqual(_apply(lift, _features(width, state)),
                                 _features(width, image))
        s_lift = _unitary_lift(width, native_columns(width, [('S', 0)]))
        carry = (1 << (2 * width)) | 1  # t0 times the input X label.
        carry_column = s_lift[_feature_masks(width).index(carry)]
        self.assertTrue(carry_column >> (len(s_lift) - 1) & 1)
        evaluations = [_features(width, state) for state in range(64)]
        self.assertEqual(_rank(evaluations), len(s_lift))

    def test_composition_and_faithfulness_are_matrix_statements(self):
        left = native_columns(2, [('H', 0), ('S', 1), ('CX', 0, 1)])
        right = native_columns(2, [('H', 1), ('S', 0)])
        self.assertEqual(_unitary_lift(2, multiply(left, right)),
                         _compose(_unitary_lift(2, left), _unitary_lift(2, right)))
        self.assertEqual(_unitary_lift(2, scale(left, CI)), _unitary_lift(2, left))
        cliffords = _one_qubit_cliffords()
        self.assertEqual(len(cliffords), 24)
        for a, ua in cliffords.items():
            for b, ub in cliffords.items():
                self.assertEqual(a == b, _scalar(multiply(adjoint(ua), ub)))

    def test_literal_affine_query_has_degree_two_with_every_input_phase(self):
        query_lifts = []
        for x in range(8):
            lift = _lift(3, lambda state: _shear_action(x, state))
            query_lifts.append(lift)
            for state in range(256):
                image = _shear_action(x, state)
                a, b, phase = state & 7, (state >> 3) & 7, state >> 6
                ia, ib, ip = image & 7, (image >> 3) & 7, image >> 6
                for basis in range(8):
                    before = _shear_basis(x, basis)
                    actual_basis = _shear_basis(x, before ^ a)
                    actual_phase = (phase + 2 * ((b & before).bit_count() & 1)) % 4
                    self.assertEqual(actual_basis, basis ^ ia)
                    self.assertEqual(actual_phase,
                                     (ip + 2 * ((ib & basis).bit_count() & 1)) % 4)
            self.assertEqual(_compose(lift, lift),
                             tuple(1 << i for i in range(len(lift))))
        self.assertEqual(_degree([_pack(lift) for lift in query_lifts]), 2)
        wall = _unitary_lift(3, native_columns(3, [('H', 0), ('S', 2)]))
        inverse_wall = _unitary_lift(3, native_columns(3, [('S', 2)] * 3 + [('H', 0)]))
        r = _unitary_lift(3, native_columns(3, [('H', 1)]))
        orbit = [_compose(_compose(_compose(_compose(wall, query), r), query),
                          inverse_wall) for query in query_lifts]
        self.assertLessEqual(_degree([_pack(value) for value in orbit]), 4)
        # Independent address variables make the two-query degree bound nonvacuous.
        two_query_words = []
        for x in range(64):
            first, second = query_lifts[x & 7], query_lifts[x >> 3]
            forward = _compose(second, _compose(wall, first))
            inverse = _compose(first, _compose(inverse_wall, second))
            self.assertEqual(_compose(forward, inverse),
                             tuple(1 << i for i in range(len(forward))))
            two_query_words.append(_pack(forward))
        self.assertEqual(_degree(two_query_words), 4)

    def test_affine_symplectic_entries_do_not_bound_arbitrary_sign_lifts(self):
        values = []
        for x in range(8):
            nonlinear = int(x == 7)
            action = lambda state: state ^ (2 * nonlinear * (state & 1) << 2)
            values.append(_pack(_lift(1, action)))
            for label in range(4):
                self.assertEqual(action(label) & 3, label)
        self.assertEqual(_degree(values), 3)

    def test_square_zero_factorization_and_monomial_collection(self):
        for count in range(9):
            factors = [1 << (1 << i) for i in range(count)]
            actual = 1
            for element in factors:
                actual = _algebra_multiply(actual, 1 ^ element)
            for left, right in combinations(factors, 2):
                actual = _algebra_multiply(actual, 1 ^ _algebra_multiply(left, right))
            self.assertEqual(actual, 1 ^ _xor(factors))
        # Degree-two address family with nontrivial, dependent algebra coefficients.
        coefficients = {0: 2 ^ 8, 1: 4 ^ 32, 2: 16 ^ 128,
                        4: 256 ^ 1024, 3: 2 ^ 4, 5: 16 ^ 256, 6: 8 ^ 32}
        grouped = {}
        for address_mask, element in coefficients.items():
            grouped[address_mask] = _algebra_multiply(grouped.get(address_mask, 1),
                                                       1 ^ element)
            self.assertEqual(_algebra_multiply(element, element), 0)
        for (a, left), (b, right) in combinations(coefficients.items(), 2):
            grouped[a | b] = _algebra_multiply(grouped.get(a | b, 1),
                                               1 ^ _algebra_multiply(left, right))
        values = []
        for x in range(8):
            actual = 1 ^ _xor(element for mask, element in coefficients.items()
                              if mask & x == mask)
            product = 1
            for mask, generator in grouped.items():
                if mask & x == mask:
                    product = _algebra_multiply(product, generator)
            self.assertEqual(actual, product)
            values.append(actual)
        self.assertLessEqual(len(_group(values)), 1 << len(grouped))

    def test_group_rank_is_not_raw_coefficient_rank(self):
        a, b, c = (1 << (1 << i) for i in range(3))
        generators = [1 ^ a, 1 ^ b, 1 ^ c, 1 ^ a ^ b ^ c]
        self.assertEqual(_xor(generators), 0)
        self.assertEqual(_rank(generators), 3)
        self.assertEqual(len(_group(generators)), 16)
        product = 1
        for value in generators:
            product = _algebra_multiply(product, value)
        self.assertNotEqual(product, 1)

    def test_nonscalar_clifford_trace_gap_and_projective_independence(self):
        cliffords = _one_qubit_cliffords()
        for matrix in cliffords.values():
            if not _scalar(matrix):
                trace = (matrix[0][0] + matrix[1][1]) * Q2(Fraction(1, 2))
                self.assertLessEqual(trace.real * trace.real + trace.imag * trace.imag,
                                     Q2(Fraction(1, 2)))
        # An anticommuting pair has commuting, independent adjoint lifts.
        lifts = [_unitary_lift(1, _normal(1, label)) for label in (1, 2)]
        identity = tuple(1 << i for i in range(len(lifts[0])))
        self.assertEqual(_compose(lifts[0], lifts[1]), _compose(lifts[1], lifts[0]))
        self.assertEqual(len({identity, lifts[0], lifts[1],
                              _compose(lifts[0], lifts[1])}), 4)
        x, z = _normal(1, 1), _normal(1, 2)
        self.assertEqual(multiply(multiply(multiply(x, z), x), z), scale(eye(2), -1))
        for p in range(1, 7):
            rows = []
            for x in range(1 << p):
                row = 0
                for y in range(1 << p):
                    value = 0
                    for i in range(p):
                        term = ((x ^ y) >> i) & 1
                        for j in range(i + 1, p):
                            term &= 1 ^ (((x >> j) & 1) & ((y >> j) & 1))
                        value ^= term
                    row |= value << y
                rows.append(row)
            self.assertEqual(_rank(rows), 1 << p)


if __name__ == '__main__':
    unittest.main()
