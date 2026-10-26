"""Exact receipts for the balanced-Clifford source and its paid phase core.

The small operators test complete action and literal phases. The analytic
normal-form argument, not these fixtures, proves the all-width statement.
"""
from itertools import product
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, CQ2, Q2, adjoint, multiply, native_columns, scale,
)
from tests.test_common_probe_source import ROOT_HALF
from tests.test_joint_lookup_pauli_source import _inverse


OMEGA = (CO + CI) * ROOT_HALF
ROOTS = [CO]
for _ in range(7):
    ROOTS.append(ROOTS[-1] * OMEGA)
ROOTS = tuple(ROOTS)


def _quadratic(linear, pairs, value):
    return (sum(coefficient * ((value >> i) & 1)
                for i, coefficient in enumerate(linear))
            + 2 * sum(((value >> i) & (value >> j)) & 1
                      for i, j in pairs)) % 4


def _quadratic_basis(linear, pairs):
    """Symmetric bilinear elimination; vectors are in the old coordinates."""
    q = lambda value: _quadratic(linear, pairs, value)
    b = lambda u, v: ((q(u ^ v) - q(u) - q(v)) % 4) // 2
    remaining = [1 << i for i in range(len(linear))]
    odd, hyperbolic, radical = [], [], []
    while remaining:
        pivot = next((v for v in remaining if b(v, v)), None)
        if pivot is not None:
            odd.append(pivot)
            remaining = [v ^ (pivot if b(pivot, v) else 0)
                         for v in remaining if v != pivot]
            continue
        pair = next(((u, v) for i, u in enumerate(remaining)
                     for v in remaining[i + 1:] if b(u, v)), None)
        if pair is None:
            radical.extend(remaining)
            break
        u, v = pair
        hyperbolic.append(pair)
        remaining = [w ^ (u if b(v, w) else 0) ^ (v if b(u, w) else 0)
                     for w in remaining if w not in pair]
    basis = odd + [v for pair in hyperbolic for v in pair] + radical
    return basis, len(odd), len(hyperbolic)


def _phase_core(k, signs, pair_count, even_wires=(), spectators=0):
    """Chronological native schedule with each physical T band explicit."""
    singletons = tuple(range(1, 1 + len(signs)))
    pairs = tuple((1 + len(signs) + 2 * j, 2 + len(signs) + 2 * j)
                  for j in range(pair_count))
    width = 1 + len(signs) + 2 * pair_count + spectators
    residue = (k + sum(signs) + pair_count) % 8
    schedule = [('C', [('S', 0)] * (residue // 2))]
    for wire in even_wires:
        schedule.append(('C', [('H', wire), ('CX', 0, wire), ('H', wire)]))
    first = [('T' if sign == 1 else 'TDG', wire)
             for sign, wire in zip(signs, singletons)]
    first += [('T', wire) for pair in pairs for wire in pair]
    if residue % 2:
        first.append(('T', 0))
    schedule.append(('T', first))
    parity = [('CX', 0, wire) for wire in singletons]
    parity += [('CX', 0, wire) for pair in pairs for wire in pair]
    second = [('TDG' if sign == 1 else 'T', wire)
              for sign, wire in zip(signs, singletons)]
    second += [('TDG', wire) for pair in pairs for wire in pair]
    schedule += [('C', parity), ('T', second), ('C', list(reversed(parity)))]
    parity = [('CX', u, v) for u, v in pairs]
    schedule += [('C', parity), ('T', [('TDG', v) for _, v in pairs]),
                 ('C', list(reversed(parity)))]
    parity = [gate for u, v in pairs for gate in [('CX', u, v), ('CX', 0, v)]]
    schedule += [('C', parity), ('T', [('T', v) for _, v in pairs]),
                 ('C', list(reversed(parity)))]
    return schedule, width, residue


def _word(schedule):
    return [gate for _, gates in schedule for gate in gates]


def _target_phase(k, signs, pair_count, even_wires, basis):
    if not basis & 1:
        return 0
    value = k + 2 * sum(sign * ((basis >> (1 + j)) & 1)
                        for j, sign in enumerate(signs))
    value += 4 * sum(((basis >> (1 + len(signs) + 2 * j))
                     & (basis >> (2 + len(signs) + 2 * j))) & 1
                    for j in range(pair_count))
    value += 4 * sum((basis >> wire) & 1 for wire in even_wires)
    return value % 8


def _diagonal(width, phases):
    return [[ROOTS[phases[column] % 8] if row == column else CZ
             for column in range(1 << width)] for row in range(1 << width)]


def _galois(matrix):
    """Fix i and reverse sqrt(2), retaining exact rational components."""
    return [[CQ2(Q2(entry.real.a, -entry.real.b),
                 Q2(entry.imag.a, -entry.imag.b))
             for entry in row] for row in matrix]


def _pauli(width, x, z):
    phase = ROOTS[2 * (x & z).bit_count() % 8]
    return [[phase * (-1 if (column & z).bit_count() % 2 else 1)
             if row == column ^ x else CZ
             for column in range(1 << width)] for row in range(1 << width)]


def _displacement_rank(matrix, width):
    """Recover all generator images exactly; reject a non-Clifford matrix."""
    paulis = [_pauli(width, label & ((1 << width) - 1), label >> width)
              for label in range(1 << (2 * width))]
    columns = []
    for bit in range(2 * width):
        image = multiply(multiply(matrix, paulis[1 << bit]), adjoint(matrix))
        label = next((i for i, candidate in enumerate(paulis)
                      if image == candidate or image == scale(candidate, -1)), None)
        if label is None:
            raise AssertionError('non-Clifford Pauli image')
        columns.append(label ^ (1 << bit))
    pivots = {}
    for column in columns:
        while column:
            pivot = column.bit_length() - 1
            if pivot in pivots:
                column ^= pivots[pivot]
            else:
                pivots[pivot] = column
                break
    return len(pivots)


class CliffordInvolutionSourceTests(unittest.TestCase):
    def test_quadratic_compression_includes_odd_hyperbolic_and_radical_parts(self):
        fixtures = [([0, 0, 0], []),
                    ([2, 0, 2, 0], [(0, 1)]),
                    ([1, 3, 0, 2, 0], [(2, 3)]),
                    ([1, 0, 2, 3, 2, 0],
                     [(0, 1), (0, 3), (1, 2), (2, 4), (3, 5), (4, 5)]),
                    ([0, 2, 0, 2, 0, 2],
                     [(0, 1), (0, 2), (1, 3), (2, 4), (3, 5), (4, 5)])]
        for linear, pairs in fixtures:
            basis, odd, hyperbolic = _quadratic_basis(linear, pairs)
            self.assertEqual(len(basis), len(linear))
            values = set()
            for value in range(1 << len(linear)):
                old = 0
                for i, vector in enumerate(basis):
                    if value >> i & 1:
                        old ^= vector
                values.add(old)
                expected = sum(_quadratic(linear, pairs, vector) * ((value >> i) & 1)
                               for i, vector in enumerate(basis))
                expected += 2 * sum(((value >> (odd + 2 * j))
                                     & (value >> (odd + 2 * j + 1))) & 1
                                    for j in range(hyperbolic))
                self.assertEqual(_quadratic(linear, pairs, old), expected % 4)
            self.assertEqual(len(values), 1 << len(linear))
            self.assertTrue(all(_quadratic(linear, pairs, v) % 2
                                for v in basis[:odd]))
            self.assertTrue(all(not _quadratic(linear, pairs, v) % 2
                                for v in basis[odd:]))

    def test_four_band_schedule_literal_phase_and_all_input_columns(self):
        fixtures = [(0, (), 0, (), 0), (7, (), 0, (), 1),
                    (3, (1, -1), 0, (3,), 1),
                    (0, (), 2, (1, 4), 0),
                    (5, (-1,), 1, (2, 4), 1)]
        for k, signs, pairs, even, spectators in fixtures:
            schedule, width, residue = _phase_core(k, signs, pairs, even, spectators)
            word = _word(schedule)
            for kind, gates in schedule:
                if kind == 'T':
                    targets = [gate[1] for gate in gates]
                    self.assertEqual(len(targets), len(set(targets)))
            t_count = sum(gate[0] in ('T', 'TDG') for gate in word)
            self.assertEqual(t_count, 2 * len(signs) + 6 * pairs + residue % 2)
            self.assertLessEqual(t_count, 3 * width - 2)
            self.assertLessEqual(sum(bool(gates) for kind, gates in schedule if kind == 'T'), 4)
            self.assertEqual(len(word) - t_count,
                             2 * len(signs) + 10 * pairs + 3 * len(even) + residue // 2)
            phases = [_target_phase(k, signs, pairs, even, basis)
                      for basis in range(1 << width)]
            expected = _diagonal(width, phases)
            self.assertEqual(native_columns(width, word), expected)
            self.assertEqual(native_columns(width, _inverse(word)), adjoint(expected))

    def test_eighth_root_offsets_and_controlled_s_signs_are_not_projective(self):
        for k, sign in product(range(8), (1, -1)):
            schedule, width, _ = _phase_core(k, (sign,), 0)
            phases = [_target_phase(k, (sign,), 0, (), value) for value in range(4)]
            source = native_columns(width, _word(schedule))
            self.assertEqual(source, _diagonal(width, phases))
            # The derivative phase is g(0,z) on t=0 and its conjugate on t=1.
            x = native_columns(width, [('H', 0), ('S', 0), ('S', 0), ('H', 0)])
            image = multiply(multiply(source, x), adjoint(source))
            expected = [[ROOTS[((1 - 2 * (column & 1))
                                * (k + 2 * sign * ((column >> 1) & 1))) % 8]
                         if row == (column ^ 1) else CZ
                         for column in range(4)] for row in range(4)]
            self.assertEqual(image, expected)
            self.assertEqual(image, adjoint(image))
            self.assertEqual(multiply(image, image), _diagonal(width, [0] * 4))

    def test_independent_nonmonomial_cliffords_and_actual_source_inverses(self):
        # A Bell change diagonalizes SWAP, while H on the root turns Z into X.
        # U Z_root U† is the independently specified Z_root * SWAP(1,2).
        core, width, _ = _phase_core(0, (), 1)
        swap_source = [('H', 0)] + _word(core) + [('H', 0), ('H', 1), ('CX', 1, 2)]
        swap_target = [('S', 0), ('S', 0), ('CX', 1, 2), ('CX', 2, 1), ('CX', 1, 2)]
        # H S H maps the (X+Y)/sqrt(2) axis to H=(X+Z)/sqrt(2).
        # Root T combines literally with the shared CCZ root phase.
        core, _, _ = _phase_core(1, (), 1)
        hadamard_source = [('H', 0)] + _word(core) + [('H', 0), ('S', 0), ('H', 0)]
        hadamard_target = [('H', 0), ('H', 2), ('CX', 1, 2), ('H', 2)]
        probe = native_columns(width, [('S', 0), ('S', 0)])
        identity = _diagonal(width, [0] * (1 << width))
        for source_word, target_word in ((swap_source, swap_target),
                                         (hadamard_source, hadamard_target)):
            source = native_columns(width, source_word)
            inverse = native_columns(width, _inverse(source_word))
            target = native_columns(width, target_word)
            self.assertEqual(inverse, adjoint(source))
            self.assertEqual(multiply(inverse, source), identity)
            self.assertEqual(multiply(multiply(source, probe), inverse), target)
            self.assertEqual(target, adjoint(target))
            self.assertEqual(sum((target[i][i] for i in range(1 << width)), CZ), CZ)

    def test_hermitian_clifford_galois_sign_matches_displacement_rank(self):
        fixtures = [(1, [('H', 0)], 1),
                    (2, [('H', 0), ('H', 1)], 2),
                    (3, [('H', 0), ('H', 1), ('H', 2)], 3),
                    (3, [('Z', 0)], 0),
                    (3, [('Z', 0), ('CX', 1, 2), ('CX', 2, 1), ('CX', 1, 2)], 2)]
        for width, word, rank in fixtures:
            matrix = native_columns(width, word)
            self.assertEqual(matrix, adjoint(matrix))
            self.assertEqual(_displacement_rank(matrix, width), rank)
            self.assertEqual(_galois(matrix), scale(matrix, (-1) ** rank))
            support = []
            for x, z in product(range(1 << width), repeat=2):
                product_matrix = multiply(_pauli(width, x, z), matrix)
                coefficient = sum((product_matrix[i][i] for i in range(1 << width)), CZ)
                coefficient /= 1 << width
                if coefficient != CZ:
                    self.assertEqual(coefficient.imag, 0)
                    self.assertEqual(coefficient.abs2(), Q2(1) / (1 << rank))
                    support.append(coefficient)
            self.assertEqual(len(support), 1 << rank)

        # An even-rank nonmonomial Clifford can have no nonidentity shared Z.
        pair = native_columns(2, [('H', 0), ('H', 1)])
        z_paulis = [_pauli(2, 0, z) for z in range(4)]
        shared = [z for z in range(4)
                  if any(multiply(multiply(pair, z_paulis[z]), pair) == scale(p, sign)
                         for p in z_paulis for sign in (1, -1))]
        self.assertEqual(shared, [0])
        self.assertTrue(all(entry != CZ for row in pair for entry in row))

    def test_outer_galois_defect_is_the_literal_combined_clifford(self):
        width = 3
        c2 = native_columns(width, [('H', 0), ('CX', 0, 1), ('S', 2), ('H', 1)])
        c3 = native_columns(width, [('CX', 1, 2), ('H', 2), ('S', 0)])
        identity = _diagonal(width, [0] * 8)
        for odd in (False, True):
            phases2, phases3 = [], []
            for value in range(8):
                x, y, z = [(value >> bit) & 1 for bit in range(3)]
                phases2.append((1 if odd else 2) * x + (3 if odd else 6) * y
                               + 2 * x * z + 4 * x * y * z)
                phases3.append((7 if odd else 4) * z + 2 * y * z + 4 * x * y * z)
            d2, d3 = _diagonal(width, phases2), _diagonal(width, phases3)
            z2, z3 = _pauli(width, 0, 3 if odd else 0), _pauli(width, 0, 4 if odd else 0)
            self.assertEqual(_galois(d2), multiply(z2, d2))
            self.assertEqual(_galois(d3), multiply(z3, d3))
            a = multiply(multiply(multiply(c3, d3), c2), d2)
            e = multiply(multiply(c2, z2), adjoint(c2))
            e = multiply(multiply(d3, e), adjoint(d3))
            e = multiply(multiply(multiply(c3, z3), e), adjoint(c3))
            # The three H gates contribute the literal common minus sign.
            self.assertEqual(multiply(_galois(a), adjoint(a)), scale(e, -1))
            _displacement_rank(e, width)  # Checks every Pauli generator image.
            if not odd:
                self.assertEqual(e, identity)

    def test_odd_linear_phase_allows_odd_rank_with_required_anticommutation(self):
        h = native_columns(1, [('H', 0)])
        t = native_columns(1, [('T', 0)])
        x, z = _pauli(1, 1, 0), _pauli(1, 0, 1)
        a = multiply(h, t)
        k = multiply(multiply(adjoint(t), x), t)
        self.assertEqual(k, adjoint(k))
        self.assertEqual(_displacement_rank(k, 1), 1)
        self.assertEqual(_galois(k), scale(k, -1))
        self.assertEqual(multiply(multiply(a, k), adjoint(a)), z)
        e = multiply(multiply(h, z), h)
        self.assertEqual(e, x)
        self.assertEqual(multiply(_galois(a), adjoint(a)), scale(e, -1))
        self.assertEqual(multiply(multiply(e, z), adjoint(e)), scale(z, -1))


if __name__ == '__main__':
    unittest.main()
