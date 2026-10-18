"""Exact two-band quadratic XORs with one arbitrary returned dirty bit.

These fixtures emit every CNOT, H, S, T and T-dagger.  The quadratic
coordinate change is invertible, needs no initialized workspace, and is
undone before the affine correction.  No shared-control batching claim
is inferred from this one-output construction.
"""
from collections import Counter
from itertools import combinations
from math import comb
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_bent_query_baseline import _compact_toffoli
from tests.test_dirty_indicator_depth import _parity_basis
from tests.test_t_depth import _word
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import (
    _primitive_inverse, _reduce_adjacent_s,
)


def _actual_inverse(word):
    return _reduce_adjacent_s(_primitive_inverse(word))


def _odd_parity_band(coordinates, coefficients):
    """One physical band, with each odd coefficient's even part paid by S."""
    even, odd = [], []
    for wire, coefficient in coefficients:
        assert coefficient % 2
        power, name = min((((coefficient - base) % 8) // 2, name)
                          for base, name in ((1, 'T'), (7, 'TDG')))
        even.extend([('S', wire)] * power)
        odd.append((name, wire))
    assert len({wire for _, wire in odd}) == len(odd)
    return coordinates + even + odd + list(reversed(coordinates))


def _symmetric_quadratic_query(address, target, dirty):
    """Literal XOR of e2 on positive even arity; zero clean, one dirty."""
    n = len(address)
    assert n > 0 and n % 2 == 0
    assert len(set(address + [target, dirty])) == n + 2
    first = [('CX', dirty, wire) for wire in address + [target]]
    second = ([('CX', target, dirty)]
              + [('CX', dirty, wire) for wire in address]
              + [('CX', wire, target) for wire in address]
              + [('CX', dirty, target)]
              + [('CX', wire, dirty) for wire in address])
    a = _odd_parity_band(first, [(wire, 1) for wire in address]
                         + [(target, n - 1), (dirty, 1 - n)])
    b = _odd_parity_band(second, [(wire, -1) for wire in address]
                         + [(target, -1), (dirty, 1)])
    return [('H', target)] + a + b + [('H', target)]


def _polar_basis(rows):
    """Columns for a symplectic-pair basis followed by the radical."""
    def pairing(a, b):
        return sum((rows[i] & b).bit_count()
                   for i in range(len(rows)) if a >> i & 1) % 2

    remaining = [1 << i for i in range(len(rows))]
    pairs = []
    while True:
        found = next(((i, j) for i in range(len(remaining))
                      for j in range(i + 1, len(remaining))
                      if pairing(remaining[i], remaining[j])), None)
        if found is None:
            return pairs + remaining, len(pairs)
        i, j = found
        a, b = remaining[i], remaining[j]
        pairs.extend((a, b))
        remaining = [w ^ (a if pairing(w, b) else 0)
                     ^ (b if pairing(w, a) else 0)
                     for k, w in enumerate(remaining) if k not in (i, j)]


def _rows_from_columns(columns):
    return [sum(((column >> i) & 1) << j
                for j, column in enumerate(columns))
            for i in range(len(columns))]


def _inverse_rows(rows):
    n = len(rows)
    augmented = [row | (1 << (n + i)) for i, row in enumerate(rows)]
    for column in range(n):
        pivot = next(i for i in range(column, n)
                     if augmented[i] >> column & 1)
        augmented[pivot], augmented[column] = augmented[column], augmented[pivot]
        for i in range(n):
            if i != column and augmented[i] >> column & 1:
                augmented[i] ^= augmented[column]
    return [row >> n for row in augmented]


def _quadratic_coordinates(n, pairs):
    polar = [0] * n
    for i, j in pairs:
        polar[i] ^= 1 << j
        polar[j] ^= 1 << i
    basis, rank = _polar_basis(polar)
    inverse = _inverse_rows(_rows_from_columns(basis))
    if not rank:
        return [1 << i for i in range(n)], 0
    standard_e2 = [((1 << rank) - 1) ^ (1 << i) for i in range(rank)]
    target_columns, target_rank = _polar_basis(standard_e2)
    assert target_rank == rank
    target_rows = _rows_from_columns(target_columns)
    target_rows += [1 << i for i in range(rank, n)]
    coordinates = []
    for row in target_rows:
        out = 0
        for i in range(n):
            if row >> i & 1:
                out ^= inverse[i]
        coordinates.append(out)
    return coordinates, rank


def _apply_rows(rows, x):
    return sum(((row & x).bit_count() % 2) << i
               for i, row in enumerate(rows))


def _quadratic_value(x, pairs, linear=0, constant=0):
    return (sum(((x >> i) & 1) * ((x >> j) & 1) for i, j in pairs)
            + (x & linear).bit_count() + constant) % 2


def _quadratic_query(n, pairs, linear=0, constant=0):
    """Generic rank-congruence emitter, using address 0..n-1 and z=n,d=n+1."""
    rows, rank = _quadratic_coordinates(n, pairs)
    coordinates = _parity_basis(list(range(n)), rows)
    correction = linear
    for i in range(n):
        y = _apply_rows(rows, 1 << i) & ((1 << rank) - 1)
        correction ^= (comb(y.bit_count(), 2) % 2) << i
    affine = [('CX', i, n) for i in range(n) if correction >> i & 1]
    if constant:
        # Literal X=HSS H, so the primitive H/S/CX ledger is explicit.
        affine += [('H', n), ('S', n), ('S', n), ('H', n)]
    core = (_symmetric_quadratic_query(list(range(rank)), n, n + 1)
            if rank else [])
    return coordinates + core + list(reversed(coordinates)) + affine, rank


def _exact_diagonal_action(basis, word):
    phase = 0
    for name, *wires in word:
        if name == 'CX':
            if basis >> wires[0] & 1:
                basis ^= 1 << wires[1]
        else:
            assert name in ('S', 'T', 'TDG')
            phase += {'S': 2, 'T': 1, 'TDG': -1}[name] * (basis >> wires[0] & 1)
    return basis, phase % 8


class QuadraticNativeQueryTests(unittest.TestCase):
    def test_even_arity_phase_identity_both_dirty_values_and_actual_inverse(self):
        # Covers every coefficient residue, both signs, and growing arity.
        for n in (2, 4, 6, 8, 10):
            word = _symmetric_quadratic_query(list(range(n)), n, n + 1)
            for native in (word, _actual_inverse(word)):
                for basis in range(1 << (n + 2)):
                    address = basis & ((1 << n) - 1)
                    expected = 4 * ((basis >> n) & 1) * (comb(address.bit_count(), 2) % 2)
                    self.assertEqual(_exact_diagonal_action(basis, native[1:-1]),
                                     (basis, expected))

    def test_two_physical_bands_and_forward_inverse_primitive_ledgers(self):
        for n in (2, 4, 6, 8, 10, 12):
            word = _symmetric_quadratic_query(list(range(n)), n, n + 1)
            shifted = n % 8 in (4, 6)
            for native, extra in ((word, 3 * shifted),
                                  (_actual_inverse(word), 5 * shifted)):
                counts = Counter(g[0] for g in native)
                self.assertEqual(counts['T'] + counts['TDG'], 2 * n + 4)
                self.assertEqual(counts['CX'], 8 * n + 6)
                self.assertEqual(counts['H'], 2)
                self.assertEqual(counts['S'], extra)
                bands = _native_t_bands(native)
                self.assertEqual(len(bands), 2)
                for band in bands:
                    self.assertEqual(sorted(g[1] for g in band), list(range(n + 2)))

    def test_polar_congruence_all_quadratic_functions_through_four_addresses(self):
        for n in range(5):
            possible = list(combinations(range(n), 2))
            for mask in range(1 << len(possible)):
                pairs = [pair for i, pair in enumerate(possible) if mask >> i & 1]
                rows, rank = _quadratic_coordinates(n, pairs)
                self.assertEqual(rank % 2, 0)
                self.assertEqual(len({_apply_rows(rows, x) for x in range(1 << n)}), 1 << n)
                residual = sum((comb((_apply_rows(rows, 1 << i) & ((1 << rank) - 1)).bit_count(), 2) % 2) << i
                               for i in range(n))
                for linear in range(1 << n):
                    for constant in (0, 1):
                        for x in range(1 << n):
                            y = _apply_rows(rows, x) & ((1 << rank) - 1)
                            actual = (comb(y.bit_count(), 2) + ((linear ^ residual) & x).bit_count() + constant) % 2
                            self.assertEqual(actual, _quadratic_value(x, pairs, linear, constant))

    def test_all_small_functions_full_columns_and_literal_native_adjoint(self):
        cases = []
        for n in range(3):
            possible = list(combinations(range(n), 2))
            for mask in range(1 << len(possible)):
                pairs = [pair for i, pair in enumerate(possible) if mask >> i & 1]
                cases.extend((n, pairs, linear, constant)
                             for linear in range(1 << n) for constant in (0, 1))
        cases += [(4, [(0, 1), (1, 2)], 9, 1),
                  (4, [(0, 1), (2, 3), (0, 2)], 7, 1),
                  (6, list(combinations(range(6), 2)), 37, 1)]
        for n, pairs, linear, constant in cases:
            word, rank = _quadratic_query(n, pairs, linear, constant)
            size = 1 << (n + 2)
            expected = [[CZ for _ in range(size)] for _ in range(size)]
            for column in range(size):
                row = column ^ (_quadratic_value(column & ((1 << n) - 1), pairs, linear, constant) << n)
                expected[row][column] = CO
            actual = native_columns(n + 2, word)
            self.assertEqual(actual, expected)
            self.assertEqual(native_columns(n + 2, _actual_inverse(word)), adjoint(actual))
            self.assertEqual(len(_native_t_bands(word)), 2 if rank else 0)
            self.assertLessEqual(max((wire for gate in word for wire in gate[1:]), default=0), n + 1)

    def test_symmetric_cubic_fair_baseline_complete_occupied_action(self):
        m = 7
        target, pivot, dirty, clean = 0, m, m + 1, m + 2
        coordinates = [('CX', wire, pivot) for wire in range(1, m)]
        query = _symmetric_quadratic_query(list(range(1, m)), clean, dirty)
        middle = _word(_compact_toffoli(pivot, clean, target, dirty))
        # The query is a literal full-space permutation and its own actual
        # inverse.  This reuse is justified for both dirty values, unlike an
        # unspecified clean-one extension of an initialized phase source.
        word = coordinates + query + middle + query + list(reversed(coordinates))
        counts = Counter(g[0] for g in word)
        self.assertEqual(counts, Counter({'CX': 136, 'H': 6, 'S': 6, 'T': 22, 'TDG': 18}))
        self.assertEqual(len(_native_t_bands(word)), 6)
        inverse = _actual_inverse(word)
        self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in inverse), 152)
        self.assertEqual(len(_native_t_bands(inverse)), 6)
        actual = native_columns(m + 3, word, range(1 << (m + 2)))
        for column in range(1 << (m + 2)):
            address = (column >> 1) & ((1 << m) - 1)
            output = column ^ (comb(address.bit_count(), 3) % 2)
            for row in range(1 << (m + 3)):
                self.assertEqual(actual[row][column], CO if row == output else CZ)
        # Retain selected clean-one columns as well: the physical extension
        # is the explicit phase-free XOR e3(x)+c*e1(x), also self-inverse.
        columns = [0, 1, 127, 128, 255, 256, 511, 512, 513, 767, 768, 1023]
        adjoint_columns = native_columns(m + 3, inverse, columns)
        for position, column in enumerate(columns):
            address = (column >> 1) & ((1 << m) - 1)
            bit = (comb(address.bit_count(), 3)
                   + ((column >> clean) & 1) * address.bit_count()) % 2
            for row in range(1 << (m + 3)):
                self.assertEqual(adjoint_columns[row][position],
                                 CO if row == column ^ bit else CZ)


if __name__ == '__main__':
    unittest.main()
