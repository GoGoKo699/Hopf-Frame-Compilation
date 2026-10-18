"""Exact full-input HRT extraction and fast dirty subset-convolution fixtures.

No fixture initializes a coefficient, mask or linear-program temporary.
Bounded checks support, rather than replace, the general proof.
"""
from itertools import product
from math import comb
from random import Random
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_balanced_algebra_permutation import (
    _hrt_f, _hrt_k, _inverse, _ring_product, _tree_nodes,
)
from tests.test_returned_diagonal_echo_degree import _ccz
from tests.test_two_layer_diagonal_rigidity import _primitive_inverse, _reduce_adjacent_s


def _reference(ref, x, nodes):
    if ref is None:
        return 0
    return ((x >> ref[1]) & 1) if ref[0] == 'x' else nodes[ref[1]]


def _ranked_program(t, kind):
    """Ranked subset/superset zeta as a straight-line XOR program."""
    m = 1 << t
    current = [[('x', subset) if subset.bit_count() == rank else None
                for subset in range(m)] for rank in range(t + 1)]
    nodes = []
    for bit in range(t):
        for rank in range(t + 1):
            for subset in range(m):
                active = bool(subset & (1 << bit))
                if active != (kind != 'W'):
                    continue
                other = subset ^ (1 << bit)
                nodes.append((current[rank][subset], current[rank][other]))
                current[rank][subset] = ('h', len(nodes) - 1)
    outputs = [current[i if kind == 'U' else j if kind == 'V' else i + j][subset]
               for subset in range(m)
               for i in range(t + 1) for j in range(t + 1 - i)]
    return nodes, outputs


def _evaluate_program(program, x):
    nodes, outputs = program
    values = []
    for left, right in nodes:
        values.append(_reference(left, x, values) ^ _reference(right, x, values))
    return [_reference(ref, x, values) for ref in outputs]


def _linear_shear(program, inputs, outputs, temporary):
    nodes, references = program
    def source(ref):
        if ref is None:
            return None
        return (inputs if ref[0] == 'x' else temporary)[ref[1]]
    compute, internal = [], []
    for index, refs in enumerate(nodes):
        for ref in refs:
            wire = source(ref)
            if wire is not None:
                compute.append(('CX', wire, temporary[index]))
                if ref[0] == 'h':
                    internal.append(('CX', wire, temporary[index]))
    copy, internal_copy = [], []
    for target, ref in zip(outputs, references):
        wire = source(ref)
        if wire is not None:
            copy.append(('CX', wire, target))
            if ref[0] == 'h':
                internal_copy.append(('CX', wire, target))
    return (compute + copy + _inverse(compute)
            + internal + internal_copy + _inverse(internal))


def _phase_column(word, column):
    state, phase = column, 0
    for gate in word:
        if gate[0] == 'CX':
            state ^= ((state >> gate[1]) & 1) << gate[2]
        elif gate[0] in ('T', 'TDG'):
            phase += (1 if gate[0] == 'T' else -1) * ((state >> gate[1]) & 1)
        else:
            raise ValueError(gate)
    return state, phase % 8


def _fast_phase(t):
    m = 1 << t
    rank = comb(t + 2, 2) * m
    programs = [_ranked_program(t, kind) for kind in ('U', 'V', 'W')]
    size = len(programs[0][0])
    masks = [list(range(3 * m + bank * rank, 3 * m + (bank + 1) * rank))
             for bank in range(3)]
    temporary = list(range(3 * m + 3 * rank, 3 * m + 3 * rank + size))
    shifts = [_linear_shear(programs[bank], list(range(bank * m, (bank + 1) * m)),
                            masks[bank], temporary) for bank in range(3)]
    # Each CCZ has the same schedule shape. Parallelize only disjoint triples.
    schedules = [_ccz(masks[0][i], masks[1][i], masks[2][i]) for i in range(rank)]
    phase_word = [gate for position in range(len(schedules[0]))
                  for schedule in schedules for gate in schedule[position][1]]
    word, orientation = [], [False, False, False]
    for bank in (0, 1, 0, 2, 0, 1, 0, 2):
        word += phase_word
        word += _inverse(shifts[bank]) if orientation[bank] else shifts[bank]
        orientation[bank] = not orientation[bank]
    assert orientation == [False, False, False]
    return word, 3 * m + 3 * rank + size


def _source_conjugation(y, address, r, inverse=False):
    # J = F X^iota(x) F†; actual reverse has the same ideal action.
    shift = sum(((address >> i) & 1) << (1 << i) for i in range(r - 1))
    return _hrt_f(_hrt_k(y, r) ^ shift, r)


def _extract_phase(y, address, table_coefficients, output, r, inverse=False):
    top = (1 << r) - 1
    v = sum(((table_coefficients >> subset) & 1) << (top ^ subset)
            for subset in range(1, 1 << (r - 1)))
    # C in chronological order: J, Z, J†, Z. Its actual inverse reverses it.
    c_word = ['J', 'Z', 'Jinv', 'Z']
    word = c_word + ['L'] + list(reversed(c_word)) + ['L']
    if inverse:
        word = list(reversed(word))
    phase = 0
    for gate in word:
        if gate in ('J', 'Jinv'):
            y = _source_conjugation(y, address, r, inverse=(gate == 'Jinv'))
        elif gate == 'Z':
            phase ^= (y >> top) & 1
        else:
            y ^= output * v
    return y, phase


class AlgebraQueryExtractionTests(unittest.TestCase):
    def test_full_input_table_extraction_and_actual_reverse(self):
        rng = Random(20261019)
        for r in (2, 3, 4):
            a = r - 1
            cases = (range(1 << ((1 << r) - 1)) if r <= 3 else
                     [rng.randrange(1 << ((1 << r) - 1)) for _ in range(32)])
            tables = range(1 << (1 << a)) if r <= 3 else [rng.randrange(256) for _ in range(16)]
            for raw, address, table, output in product(cases, range(1 << a), tables, (0, 1)):
                y = raw << 1
                target = output * (sum((table >> subset) & 1
                                       for subset in range(1, 1 << a)
                                       if subset & address == subset) % 2)
                self.assertEqual(_extract_phase(y, address, table, output, r), (y, target))
                self.assertEqual(_extract_phase(y, address, table, output, r, True), (y, target))

    def test_native_table_query_all_physical_columns_and_ledger(self):
        # r=2: HRT is exactly one Toffoli on the three arbitrary
        # coefficient wires. Address 0 and output 1 are separate.
        source = ([('H', 4)] + [g for _, gates in _ccz(2, 3, 4) for g in gates]
                  + [('H', 4)])
        jword = _inverse(source) + [('CX', 0, 2)] + source
        zword = [('S', 4), ('S', 4)]
        cword = jword + zword + _inverse(jword) + zword
        cinverse = _reduce_adjacent_s(_primitive_inverse(cword))
        for table in range(4):
            beta, coefficient = table & 1, (table >> 1) & 1
            shift = [('CX', 1, 3)] if coefficient else []
            word = ([('H', 1)] + cword + shift + cinverse + _inverse(shift)
                    + [('H', 1)])
            if beta:
                word += [('H', 1), ('S', 1), ('S', 1), ('H', 1)]
            inverse = _reduce_adjacent_s(_primitive_inverse(word))
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 8 * 7)
            self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in word),
                             8 * 12 + 4 + 2 * coefficient + 10 + 4 * beta)
            actual = native_columns(5, word)
            for column in range(32):
                expected = column ^ ((beta ^ (coefficient * (column & 1))) << 1)
                for row in range(32):
                    self.assertEqual(actual[row][column], CO if row == expected else CZ)
            self.assertEqual(native_columns(5, inverse), adjoint(actual))
            self.assertLessEqual(len(inverse), len(word))

    def test_dirty_linear_program_all_columns_and_direct_references(self):
        for kind in ('U', 'V', 'W'):
            program = _ranked_program(1, kind)
            nodes, refs = program
            self.assertEqual(len(nodes), 2)
            self.assertEqual(len(refs), 6)
            word = _linear_shear(program, [0, 1], list(range(2, 8)), [8, 9])
            self.assertLessEqual(len(word), 8 * len(nodes) + 2 * len(refs))
            for column in range(1 << 10):
                expected_mask = sum(bit << (index + 2)
                                    for index, bit in enumerate(_evaluate_program(program, column & 3)))
                self.assertEqual(_phase_column(word, column), (column ^ expected_mask, 0))
                self.assertEqual(_phase_column(word + _inverse(word), column), (column, 0))
        # No node exists at t=0. Its original-input output must be copied
        # only in the first half of the dirty-linear echo.
        for kind in ('U', 'V', 'W'):
            self.assertEqual(_linear_shear(_ranked_program(0, kind), [0], [1], []),
                             [('CX', 0, 1)])

    def test_ranked_factorization_exact_binary_products(self):
        rng = Random(20261019)
        for t in range(5):
            m = 1 << t
            programs = [_ranked_program(t, kind) for kind in ('U', 'V', 'W')]
            rank = comb(t + 2, 2) * m
            size = 0 if t == 0 else t * (t + 1) * (1 << (t - 1))
            for nodes, refs in programs:
                self.assertEqual(len(nodes), size)
                self.assertEqual(len(refs), rank)
            cases = product(range(1 << m), repeat=3) if t <= 1 else (
                (rng.randrange(1 << m), rng.randrange(1 << m), rng.randrange(1 << m))
                for _ in range(128))
            for a, b, z in cases:
                u, v, w = [_evaluate_program(p, value)
                           for p, value in zip(programs, (a, b, z))]
                self.assertEqual(sum(x * y * q for x, y, q in zip(u, v, w)) % 2,
                                 (_ring_product(a, b) & z).bit_count() % 2)

    def test_full_native_small_product_and_actual_adjoint(self):
        phase, width = _fast_phase(0)
        self.assertEqual(width, 6)
        word = [('H', 2)] + phase + [('H', 2)]
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 56)
        self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in word), 90)
        actual = native_columns(width, word)
        for column in range(1 << width):
            target = column ^ (((column & 1) * ((column >> 1) & 1)) << 2)
            for row in range(1 << width):
                self.assertEqual(actual[row][column], CO if row == target else CZ)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))

    def test_larger_native_masks_temporary_return_and_ledger(self):
        rng = Random(20261019)
        for t in (1, 2):
            word, width = _fast_phase(t)
            m = 1 << t
            rank = comb(t + 2, 2) * m
            size = t * (t + 1) * (1 << (t - 1))
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 56 * rank)
            self.assertLessEqual(sum(g[0] == 'CX' for g in word), 96 * rank + 64 * size)
            self.assertEqual(width - 3 * m, 3 * rank + size)
            bands, active = [], []
            for gate in word:
                if gate[0] in ('T', 'TDG'):
                    active.append(gate[1])
                elif active:
                    bands.append(active)
                    active = []
            if active:
                bands.append(active)
            self.assertEqual(len(bands), 40)
            self.assertTrue(all(len(set(band)) == len(band) for band in bands))
            for _ in range(128):
                column = rng.randrange(1 << width)
                a = column & ((1 << m) - 1)
                b = (column >> m) & ((1 << m) - 1)
                z = (column >> (2 * m)) & ((1 << m) - 1)
                phase = 4 * ((_ring_product(a, b) & z).bit_count() % 2)
                self.assertEqual(_phase_column(word, column), (column, phase))
                self.assertEqual(_phase_column(_inverse(word), column), (column, phase))

    def test_fast_balanced_source_simultaneous_resources(self):
        for r in range(2, 65):
            levels = {}
            for lo, left_length, level, n, m, _ in _tree_nodes(r):
                t = r - lo - left_length - 1
                self.assertEqual(m, 1 << t)
                rank = comb(t + 2, 2) * m
                size = 0 if t == 0 else t * (t + 1) * (1 << (t - 1))
                self.assertLessEqual(size, rank)
                actual_t = (4 * n + 2) * 56 * rank
                actual_g = ((4 * n + 2) * (96 * rank + 64 * size + 2 * m)
                            + 2 * n * m + 4 * m - 2)
                work = m + (n + 1) * (3 * rank + size)
                self.assertLessEqual(actual_t, 336 * n * rank)
                self.assertLessEqual(actual_g, 978 * n * rank)
                self.assertLessEqual(work, 9 * n * rank)
                levels[level] = levels.get(level, 0) + n * rank
            bound = comb(r, 2) * (1 << r)
            height = (r - 1).bit_length()
            self.assertEqual(len(levels), height)
            self.assertLess(max(levels.values()), bound)
            self.assertLess(sum(levels.values()), height * bound)

    def test_full_degree_table_pairing_has_no_annihilator(self):
        # Full-degree coefficient is one; lower coefficients vary freely.
        for a in range(1, 4):
            m, top = 1 << a, (1 << a) - 1
            for lower in range(1 << top):
                functional = lower | (1 << top)
                for u in range(1, 1 << m):
                    support = [subset for subset in range(m) if u >> subset & 1]
                    chosen = min(support, key=int.bit_count)
                    paired = _ring_product(u, 1 << (top ^ chosen))
                    self.assertEqual((paired & functional).bit_count() % 2, 1)


if __name__ == '__main__':
    unittest.main()
