"""Exact small fixtures for the central-polarization oracle interface.

The rank statements count character queries with character-independent
interlayers, not native T layers. All matrix comparisons retain literal phases.
"""
from fractions import Fraction
from itertools import combinations
import unittest

from tests.exact_matrices import (
    CO, CQ2, CZ, Q2, Gate, adjoint, apply_word, eye, multiply, native_columns,
)


def _pairs(p):
    return tuple(combinations(range(p), 2))


def _wedge(p, u, v):
    return sum(((((u >> i) & 1) * ((v >> j) & 1))
                ^ (((u >> j) & 1) * ((v >> i) & 1))) << k
               for k, (i, j) in enumerate(_pairs(p)))


def _clique(p, x):
    return sum((((x >> i) & 1) * ((x >> j) & 1)) << k
               for k, (i, j) in enumerate(_pairs(p)))


def _rows(p, mask):
    rows = [0] * p
    for k, (i, j) in enumerate(_pairs(p)):
        if (mask >> k) & 1:
            rows[i] ^= 1 << j
            rows[j] ^= 1 << i
    return rows


def _rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


def _lower(p, c, mask):
    """Apply N_mask as a sum of coordinate lowerings to an arbitrary bank."""
    result = 0
    for source in range(1 << p):
        if not (c >> source) & 1:
            continue
        active = source & mask
        while active:
            bit = active & -active
            result ^= 1 << (source ^ bit)
            active ^= bit
    return result


def _central(p, c, mask):
    """Literal product of pair shears; no Pfaffian formula is used here."""
    for k, (i, j) in enumerate(_pairs(p)):
        if (mask >> k) & 1:
            c ^= _lower(p, _lower(p, c, 1 << i), 1 << j)
    return c


def _zeta(p, c, x):
    for i in range(p):
        if (x >> i) & 1:
            c ^= _lower(p, c, 1 << i)
    return c


def _pfaffian(rows, vertices):
    if not vertices:
        return 1
    i, rest = vertices[0], vertices[1:]
    value = 0
    for index, j in enumerate(rest):
        if (rows[i] >> j) & 1:
            value ^= _pfaffian(rows, rest[:index] + rest[index + 1:])
    return value


def _imatmul(a, b):
    return [[sum(x * y for x, y in zip(row, column))
             for column in zip(*b)] for row in a]


def _sign(character, mask):
    return (-1) ** ((character & mask).bit_count() & 1)


def _oracle_diagonal(character, masks):
    return [[_sign(character, mask) if i == j else 0
             for j in range(len(masks))] for i, mask in enumerate(masks)]


def _reverse(word):
    result = []
    for name, *qubits in reversed(word):
        if name == 'S':
            result.extend([('S', *qubits)] * 3)
        else:
            result.append(({'T': 'TDG', 'TDG': 'T'}.get(name, name), *qubits))
    return result


class CentralPolarizationTests(unittest.TestCase):
    def test_literal_cocycle_and_distinct_quadratic_refinements(self):
        # Equality on every coefficient basis vector proves these binary maps
        # equal on arbitrary bank inputs; x,y are exhaustively varied as well.
        for p in (2, 3, 4):
            for x in range(1 << p):
                for y in range(1 << p):
                    for subset in range(1 << p):
                        c = 1 << subset
                        f_cocycle = c
                        j_cocycle = c
                        for z in (y, x, x ^ y):
                            f_cocycle = _central(p, f_cocycle, _clique(p, z))
                            j_cocycle ^= _lower(p, j_cocycle, z)
                        polarization = c ^ _lower(p, _lower(p, c, y), x)
                        self.assertEqual(f_cocycle, polarization)
                        self.assertEqual(j_cocycle, polarization)
                        fc = _central(p, c, _clique(p, x))
                        self.assertEqual(fc ^ _lower(p, fc, x), _zeta(p, c, x))
            # The same cocycle does not identify its quadratic refinement.
            c = 1 << 1
            self.assertEqual(_central(p, c, _clique(p, 1)), c)
            self.assertNotEqual(c ^ _lower(p, c, 1), c)

    def test_triangle_is_a_complete_fixed_space_root_witness(self):
        for p in (3, 4, 5):
            # Repeat on several embedded triangles, not only the low three bits.
            for labels in combinations(range(p), 3):
                x = sum(1 << i for i in labels)
                c = sum(1 << ((1 << i) | (1 << j))
                        for i, j in combinations(labels, 2))
                for root in (0, 1):
                    bank = c ^ root
                    self.assertEqual(_lower(p, bank, x), 0)
                    for y in range(1 << p):
                        self.assertEqual(_lower(p, _lower(p, bank, y), x), 0)
                    transported = _central(p, bank, _clique(p, x))
                    self.assertEqual(transported, bank ^ 1)
                    lookup = sum((bank >> s) & 1 for s in range(1 << p)
                                 if s & x == s) & 1
                    self.assertEqual(lookup, root ^ 1)
                    self.assertEqual(abs((-1) ** root - (-1) ** lookup), 2)

    def test_all_central_characters_occur_in_a_free_bank_orbit(self):
        p = 4
        count = 1 << len(_pairs(p))
        seed = 1 << ((1 << p) - 1)
        orbit = [_central(p, seed, k) for k in range(count)]
        self.assertEqual(len(set(orbit)), count)
        for k, bank in enumerate(orbit):
            recovered = sum(((bank >> (((1 << p) - 1) ^ (1 << i) ^ (1 << j))) & 1)
                            << index for index, (i, j) in enumerate(_pairs(p)))
            self.assertEqual(recovered, k)
        # Independently inspect every translation in this physical orbit.
        for translation in range(count):
            for k, bank in enumerate(orbit):
                self.assertEqual(_central(p, bank, translation), orbit[k ^ translation])

    def test_physical_orbit_phase_mixes_all_central_characters(self):
        p, top = 4, 15
        count = 1 << len(_pairs(p))
        phases = []
        for k in range(count):
            lowered = _central(p, 1 << top, k)
            # Complemented labels give the algebra product g_K, including its
            # formal unit. This is an arbitrary-input basis witness, not work
            # asserted to be initialized by the compiler.
            bank = sum(((lowered >> s) & 1) << (top ^ s) for s in range(1 << p))
            phase = (-1) ** _pfaffian(_rows(p, k), tuple(range(p)))
            physical = apply_word([Gate('S', (top,)), Gate('S', (top,))], {bank: CO})
            self.assertEqual(physical, {bank: CQ2(phase)})
            phases.append(phase)
        # Exact unnormalized Walsh transform of the literal physical action.
        coefficients = list(phases)
        stride = 1
        while stride < count:
            for start in range(0, count, 2 * stride):
                for j in range(start, start + stride):
                    a, b = coefficients[j], coefficients[j + stride]
                    coefficients[j], coefficients[j + stride] = a + b, a - b
            stride *= 2
        self.assertEqual({abs(value) for value in coefficients}, {8})
        rank_four = _wedge(p, 1, 2) ^ _wedge(p, 4, 8)
        self.assertEqual(_rank(_rows(p, rank_four)), 4)
        self.assertEqual(Fraction(abs(coefficients[rank_four]), count), Fraction(1, 8))

    def test_coherent_query_fourier_support_and_direct_error(self):
        p, size = 4, 4
        count = 1 << len(_pairs(p))
        h = [[(-1) ** ((i & j).bit_count() & 1) for j in range(size)]
             for i in range(size)]  # H tensor H, with denominator two.
        first = [_wedge(p, u, v) for u, v in ((1, 2), (4, 8), (3, 12), (5, 10))]
        second = [_wedge(p, u, v) for u, v in ((1, 4), (2, 8), (9, 6), (7, 8))]
        one, two = [], []
        for character in range(count):
            d1 = _oracle_diagonal(character, first)
            d2 = _oracle_diagonal(character, second)
            one.append(_imatmul(_imatmul(h, d1), h))  # denominator four
            two.append(_imatmul(_imatmul(_imatmul(_imatmul(h, d2), h), d1), h))
        for matrices, allowed, ceiling in (
                (one, set(first), 2),
                (two, {a ^ b for a in first for b in second}, 4)):
            observed = set()
            for frequency in range(count):
                coefficient = [[sum(_sign(s, frequency) * matrices[s][i][j]
                                    for s in range(count))
                                for j in range(size)] for i in range(size)]
                if any(entry for row in coefficient for entry in row):
                    observed.add(frequency)
                    self.assertLessEqual(_rank(_rows(p, frequency)), ceiling)
            self.assertEqual(observed, allowed)
        target = _wedge(p, 1, 2) ^ _wedge(p, 4, 8)
        self.assertEqual(_rank(_rows(p, target)), 4)
        squared_errors = [sum((one[s][i][0] - (4 * _sign(s, target) if i == 0 else 0)) ** 2
                              for i in range(size)) for s in range(count)]
        self.assertEqual(sum(squared_errors), 2 * count * 16)
        # A fixed central multiplier shifts Fourier frequencies; omitting this
        # shift would incorrectly classify an already supplied target as hard.
        shifted_support = {frequency ^ target for frequency in first}
        measured_support = set()
        for frequency in range(count):
            if any(sum(_sign(s, frequency ^ target) * one[s][i][j] for s in range(count))
                   for i in range(size) for j in range(size)):
                measured_support.add(frequency)
        self.assertEqual(measured_support, shifted_support)

    def test_fixed_offset_cannot_remove_the_worst_case_full_rank(self):
        # Pfaffians use explicit perfect-matchings recursion; rank uses row
        # elimination. Thus their agreement is an independent algebra guard.
        p = 4
        for offset in range(1 << len(_pairs(p))):
            pfaffians = []
            for x in range(1 << p):
                rows = _rows(p, offset ^ _clique(p, x))
                pf = _pfaffian(rows, tuple(range(p)))
                pfaffians.append(pf)
                self.assertEqual(pf, int(_rank(rows) == p))
            self.assertEqual(sum(pfaffians) & 1, 1)  # full multilinear coefficient
        # Check the adjacent odd dimension too, with every constant offset.
        p = 5
        for offset in range(1 << len(_pairs(p))):
            self.assertEqual(max(_rank(_rows(p, offset ^ _clique(p, x)))
                                 for x in range(1 << p)), p - 1)

    def test_accepted_gram_cannot_hide_a_missing_rank_six_character(self):
        # A literal four-path H^2--SELECT--H^2 block, evaluated over all
        # characters rather than inferred from its possible path frequencies.
        p = 6
        masks = [0, _wedge(p, 1, 2), _wedge(p, 4, 8), _wedge(p, 16, 32)]
        target = masks[1] ^ masks[2] ^ masks[3]
        self.assertEqual(_rank(_rows(p, target)), 6)
        count = 1 << len(_pairs(p))
        fourier_error = 0
        maximum_error = 0
        for character in range(count):
            amplitude_numerator = sum(_sign(character, mask) for mask in masks)
            gram_numerator = amplitude_numerator ** 2  # denominator sixteen
            target_numerator = 4 * (1 + _sign(character, target))
            error = gram_numerator - target_numerator
            fourier_error += _sign(character, target) * error
            maximum_error = max(maximum_error, abs(error))
        self.assertEqual(fourier_error, -4 * count)
        self.assertGreaterEqual(maximum_error, 4)

    def test_fixed_masks_attain_the_rank_query_count(self):
        for p in range(2, 8):
            masks = []
            for j in range(p // 2):
                tail = ((1 << p) - 1) ^ ((1 << (2 * j + 2)) - 1)
                masks.append(((1 << (2 * j)) | tail,
                              (1 << (2 * j + 1)) | tail))
            for x in range(1 << p):
                result = 0
                for a, b in masks:
                    result ^= _wedge(p, x & a, x & b)
                self.assertEqual(result, _clique(p, x))
            self.assertEqual(_rank(_rows(p, _clique(p, (1 << p) - 1))), 2 * len(masks))

    def test_one_clean_unequal_native_boundaries_give_sharp_half_gram(self):
        # C=SH. Chronological C T C^dagger has its literal e^(i*pi/8)
        # scalar; the inverse rotation supplies the opposite scalar.
        rotation = [('S', 0)] * 3 + [('H', 0), ('T', 0), ('H', 0), ('S', 0)]
        prep_word = rotation + [('Z', 0), ('H', 0)]  # (HZ) C T C^dagger
        post_word = [('S', 0)] * 3 + [('H', 0), ('TDG', 0), ('H', 0), ('S', 0)]
        prep = native_columns(1, prep_word)
        post = native_columns(1, post_word)
        prep_inverse = native_columns(1, _reverse(prep_word))
        post_inverse = native_columns(1, _reverse(post_word))
        balanced = native_columns(1, [('H', 0)])
        self.assertEqual(sum(name in ('T', 'TDG') for name, *_ in prep_word), 1)
        self.assertEqual(sum(name in ('T', 'TDG') for name, *_ in post_word), 1)
        self.assertEqual(multiply(adjoint(prep), prep), eye(2))
        self.assertEqual(multiply(adjoint(post), post), eye(2))
        self.assertEqual(prep_inverse, adjoint(prep))
        self.assertEqual(post_inverse, adjoint(post))
        self.assertNotEqual(prep, adjoint(post))
        p = 4
        left, right = _wedge(p, 1, 2), _wedge(p, 4, 8)
        target = left ^ right
        coefficient = CQ2(Q2(0, Fraction(1, 4)))  # 1/(2 sqrt(2))
        for character in range(1 << len(_pairs(p))):
            a, b = _sign(character, left), _sign(character, right)
            select = [[CQ2(a), CZ], [CZ, CQ2(b)]]
            source = multiply(multiply(post, select), prep)
            self.assertEqual(multiply(adjoint(source), source), eye(2))
            accepted = source[0][0]
            self.assertEqual(accepted, coefficient * (a + b))
            self.assertEqual(accepted.conj() * accepted,
                             CQ2(Fraction(1 + _sign(character, target), 4)))
            balanced_source = multiply(multiply(balanced, select), balanced)
            balanced_accepted = balanced_source[0][0]
            self.assertEqual(balanced_accepted.conj() * balanced_accepted,
                             CQ2(Fraction(1 + _sign(character, target), 2)))
            # The exact inverse retains all rejected amplitudes and phases.
            inverse = multiply(multiply(prep_inverse, select), post_inverse)
            self.assertEqual(inverse, adjoint(source))
            self.assertEqual(multiply(inverse, source), eye(2))


if __name__ == '__main__':
    unittest.main()
