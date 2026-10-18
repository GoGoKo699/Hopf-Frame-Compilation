"""Exact native three-layer constructions and a split-specific scope test.

The cubic XOR uses one clean qubit and one arbitrary dirty qubit. Its
physical unused block is retained when checking the actual inverse.
The clean-flip example is deliberately not an optimal synthesis claim.
"""
from collections import Counter
from math import comb
import unittest

from tests.exact_matrices import CI, CO, CZ, adjoint, native_columns, scale
from tests.test_bent_query_baseline import _compact_toffoli
from tests.test_t_depth import _word
from tests.test_two_layer_diagonal_rigidity import (
    _counterexample, _inverse, _primitive_inverse, _quadratic,
    _reduce_adjacent_s,
)
from tests.test_two_layer_galois_return import _conjugate, _tau


def _quartic_correction():
    # Wire order: arbitrary t,x1,x2,x3,x4 and clean c.  M replaces
    # each of the four x coordinates by the parity of the other three.
    coordinates, _, _ = _counterexample()
    complement = ([('CX', j, 4) for j in range(1, 4)]
                  + [('CX', 4, j) for j in range(1, 4)]
                  + [('CX', j, 4) for j in range(1, 4)])
    basis = complement + coordinates
    clifford = [('S', j) for j in range(5)]
    clifford += [gate for j in range(1, 5)
                 for gate in [('H', j), ('CX', 0, j), ('H', j)]]
    word = (basis + [('T', j) for j in range(6)]
            + list(reversed(basis)) + clifford)
    return basis, word


def _three_layer_and3():
    _, _, previous = _counterexample()
    _, correction = _quartic_correction()
    return [('H', 4)] + previous + correction + [('H', 4)]


def _quartic_family(k):
    """Three physical bands for (-1)^e4, retaining the literal S orientation."""
    assert k > 0 and k % 4 == 0
    clean = k + 1
    coordinates = ([('CX', 0, j) for j in range(1, k + 1)]
                   + [('CX', j, clean) for j in range(1, k + 1)]
                   + [('CX', 0, clean)])
    complement = ([('CX', j, k) for j in range(1, k)]
                  + [('CX', k, j) for j in range(1, k)]
                  + [('CX', j, k) for j in range(1, k)])
    parity_map = complement + coordinates
    orientation = 'T' if k % 8 == 4 else 'TDG'
    band = (coordinates + [(orientation, j) for j in range(k + 2)]
            + list(reversed(coordinates)))
    first = ([('H', 0)] + band + [('H', 0)] + _inverse(band)
             + [('S', 0)] * 3 + [('H', 0), ('S', 0), ('H', 0)])
    power = 1 if k % 8 == 4 else 3
    correction = (parity_map + [('T', j) for j in range(k + 2)]
                  + list(reversed(parity_map))
                  + [gate for j in range(k + 1)
                     for gate in [('S', j)] * power]
                  + [gate for j in range(1, k + 1)
                     for gate in [('H', j), ('CX', 0, j), ('H', j)]])
    return first + correction


def _native_t_bands(word):
    bands = []
    active = []
    for gate in word:
        if gate[0] in ('T', 'TDG'):
            active.append(gate)
        elif active:
            bands.append(active)
            active = []
    if active:
        bands.append(active)
    return bands


def _basis_rows(width, word):
    rows = [1 << bit for bit in range(width)]
    for name, source, target in word:
        assert name == 'CX'
        rows[target] ^= rows[source]
    return rows


def _clean_flip_fixture():
    # Operator K=CX H_c, C=CX H_c CX, E=T_c T_d,
    # B'=C(T_c T_d†)C†; chronological words follow.
    k = [('H', 0), ('CX', 0, 1)]
    c = [('CX', 0, 1), ('H', 0), ('CX', 0, 1)]
    e = [('T', 0), ('T', 1)]
    b = c + [('T', 0), ('TDG', 1)] + c
    ua = k + e
    ub = k + e + _inverse(b)
    word = k + e + b + _inverse(e) + _inverse(k)
    return ua, ub, word


def _relative_toffoli(a, b, target):
    """Literal (-i)^(ab) X_target^(ab), four T gates in two bands."""
    return [
        ('H', target), ('CX', target, b), ('T', target), ('TDG', b),
        ('CX', target, b), ('CX', target, a), ('CX', a, b),
        ('T', b), ('TDG', a), ('CX', a, b), ('CX', target, a),
        ('H', target),
    ]


class ThreeLayerNativeMechanismTests(unittest.TestCase):
    def test_correction_independent_rows_and_exact_occupied_phase(self):
        basis, word = _quartic_correction()
        rows = _basis_rows(6, basis)
        self.assertEqual(rows, [1, 29, 27, 23, 15, 63])
        # An actual invertible CNOT word certifies six physical independent
        # axes. The sole relation after c=0 is not treated as physical.
        self.assertEqual(_basis_rows(6, basis + list(reversed(basis))),
                         [1 << bit for bit in range(6)])
        self.assertEqual(len(basis), 18)
        actual = native_columns(6, word, range(32))
        phases = [CO, CI, -CO, -CI]
        for column in range(32):
            xs = [(column >> j) & 1 for j in range(1, 5)]
            pairs = sum(xs[j] * xs[k] for j in range(4)
                        for k in range(j + 1, 4))
            exponent = (-sum(xs) + pairs
                        + 2 * (column & 1) * _quadratic(column)) % 4
            for row in range(64):
                self.assertEqual(actual[row][column],
                                 phases[exponent] if row == column else CZ)

    def test_cubic_query_all_occupied_columns_and_full_actual_inverse(self):
        word = _three_layer_and3()
        actual = native_columns(6, word)
        for column in range(32):
            output = column ^ 16 if column & 14 == 14 else column
            for row in range(64):
                self.assertEqual(actual[row][column], CO if row == output else CZ)
        # Both values of the arbitrary helper t occur in these columns.
        # Full equality on this input subspace includes its coherences with
        # arbitrary external references, with no postselection or leakage.
        inverse = _reduce_adjacent_s(_primitive_inverse(word))
        self.assertEqual(native_columns(6, inverse), adjoint(actual))
        # The occupied XOR is self-inverse, but the unspecified physical
        # clean-one extension is not. Reusing the forward word is unsound.
        self.assertNotEqual(actual, adjoint(actual))

    def test_three_physical_bands_and_complete_resource_ledger(self):
        word = _three_layer_and3()
        self.assertEqual(Counter(gate[0] for gate in word),
                         Counter({'CX': 76, 'H': 14, 'S': 9,
                                  'T': 12, 'TDG': 6}))
        inverse = _reduce_adjacent_s(_primitive_inverse(word))
        self.assertEqual(Counter(gate[0] for gate in inverse),
                         Counter({'CX': 76, 'H': 14, 'S': 19,
                                  'T': 6, 'TDG': 12}))
        for native in (word, inverse):
            bands = _native_t_bands(native)
            self.assertEqual(len(bands), 3)
            self.assertEqual([len(band) for band in bands], [6, 6, 6])
            for band in bands:
                self.assertEqual({gate[1] for gate in band}, set(range(6)))
            self.assertEqual({wire for gate in native for wire in gate[1:]},
                             set(range(6)))

    def test_fair_relative_toffoli_baseline_pays_all_six_layers(self):
        compute = _relative_toffoli(1, 2, 5)
        middle = _word(_compact_toffoli(5, 3, 4, 0))
        word = compute + middle + _inverse(compute)
        actual = native_columns(6, word)
        for column in range(32):
            output = column ^ 16 if column & 14 == 14 else column
            for row in range(64):
                self.assertEqual(actual[row][column], CO if row == output else CZ)
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))
        for native in (word, _inverse(word)):
            self.assertEqual(Counter(g[0] for g in native),
                             Counter({'CX': 28, 'H': 6, 'T': 8, 'TDG': 8}))
            self.assertEqual([len(band) for band in _native_t_bands(native)],
                             [2, 2, 4, 4, 2, 2])
        # This independently paid baseline has lower count and Clifford
        # cost than the three-layer construction. Six is only its displayed
        # schedule's depth, not a lower bound on resynthesizing this word.

    def test_eight_variable_family_negative_orientation_and_actual_inverse(self):
        self.assertEqual(_quartic_family(4), _three_layer_and3()[1:-1])
        k = 8
        word = _quartic_family(k)
        inverse = _reduce_adjacent_s(_primitive_inverse(word))
        self.assertEqual(Counter(gate[0] for gate in word),
                         Counter({'CX': 152, 'H': 20, 'S': 31,
                                  'T': 20, 'TDG': 10}))
        self.assertEqual(sum(gate[0] not in ('T', 'TDG') for gate in word), 203)
        self.assertEqual(sum(gate[0] not in ('T', 'TDG') for gate in inverse), 185)
        for native in (word, inverse):
            bands = _native_t_bands(native)
            self.assertEqual([len(band) for band in bands], [10, 10, 10])
            for band in bands:
                self.assertEqual({gate[1] for gate in band}, set(range(10)))
        self.assertEqual([gate[0] for gate in _native_t_bands(word)[0]],
                         ['TDG'] * 10)
        occupied = 1 << (k + 1)
        expected = [[CZ for _ in range(occupied)] for _ in range(2 * occupied)]
        for column in range(occupied):
            weight = (column >> 1).bit_count()
            expected[column][column] = -CO if comb(weight, 4) % 2 else CO
        # Every x input and both arbitrary dirty t values are included;
        # the complete clean-one output block is checked for zero leakage.
        # The inverse is the actual reversed native word with only literal
        # S-power simplifications, not an assumed repeat of the forward word.
        self.assertEqual(native_columns(10, word, range(occupied)), expected)
        self.assertEqual(native_columns(10, inverse, range(occupied)), expected)

    def test_three_layer_split_has_a_cheaper_complete_baseline(self):
        _, _, word = _clean_flip_fixture()
        actual = native_columns(2, word)
        baseline = [('CX', 0, 1), ('T', 0), ('TDG', 1), ('CX', 0, 1)]
        self.assertEqual(actual, native_columns(2, baseline))
        self.assertEqual(native_columns(2, _inverse(word)), adjoint(actual))
        self.assertEqual(Counter(g[0] for g in word),
                         Counter({'CX': 6, 'H': 4, 'T': 3, 'TDG': 3}))
        self.assertEqual(len(_native_t_bands(word)), 3)
        self.assertEqual(len(_native_t_bands(baseline)), 1)
        # Clean wire is c=0; the occupied action is simply T† on d.
        occupied = native_columns(2, [('TDG', 1)], (0, 2))
        self.assertEqual([[row[0], row[2]] for row in actual], occupied)

    def test_given_split_has_no_common_pauli_clean_flip(self):
        ua_word, ub_word, _ = _clean_flip_fixture()
        ua, ub = native_columns(2, ua_word), native_columns(2, ub_word)
        paulis = [_tau(label, 2) for label in range(16)]
        midpoint = _conjugate(ua, _tau(4, 2))
        self.assertEqual(midpoint, _conjugate(ub, _tau(4, 2)))
        first_candidates = []
        both_candidates = []
        for label, pauli in enumerate(paulis):
            if _conjugate(pauli, midpoint) != scale(midpoint, -1):
                continue
            first = _conjugate(adjoint(ua), pauli)
            second = _conjugate(adjoint(ub), pauli)
            is_pauli = lambda matrix: any(matrix == p or matrix == scale(p, -1)
                                          for p in paulis)
            if is_pauli(first):
                first_candidates.append(label)
                if is_pauli(second):
                    both_candidates.append(label)
        self.assertEqual(first_candidates, [4, 8])
        self.assertEqual(both_candidates, [])
        # This invalidates an inference for this displayed split only;
        # the preceding test gives its strictly simpler one-layer word.


if __name__ == '__main__':
    unittest.main()
