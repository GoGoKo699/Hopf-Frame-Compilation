"""Exact finite-field shear source and simultaneous lookup intertwiner.

The all-size character proof is in FIELD_SHEAR_INTERTWINER.md.  A separate
four-wire native p=2 word complements the dense eigenbasis fixture; it
does not give an all-size synthesis algorithm.  Native checks retain all
helper values and literal phases, including reversed physical words.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CO, CZ, Q2, adjoint, eye, multiply, native_columns, zeros,
)
from tests.test_bent_query_baseline import _compact_toffoli
from tests.test_coherent_feedback_sources import _shared_ccz
from tests.test_dirty_cubic_phase_compiler import _gray_phase_batch
from tests.test_dirty_indicator_depth import _exact_phase_action
from tests.test_joint_lookup_pauli_source import _bands, _inverse
from tests.test_t_depth import _word
from tests.test_weak_pauli_orbit_compression import _binary_rank


# Degree 1, 2, 4 and 8 irreducibles, used only for bounded exact fixtures.
_MODULI = {1: 0b11, 2: 0b111, 3: 0b10011, 4: 0x11b}


def _field_product(left, right, modulus):
    degree = modulus.bit_length() - 1
    answer = 0
    while right:
        if right & 1:
            answer ^= left
        right >>= 1
        left <<= 1
        if left >> degree:
            left ^= modulus
    return answer


def _pencil_columns(p, address):
    k = 1 << (p - 1)
    return [_field_product(address, 1 << bit, _MODULI[p])
            for bit in range(k)]


def _linear_image(columns, value):
    result = 0
    for bit, column in enumerate(columns):
        if value >> bit & 1:
            result ^= column
    return result


def _shear(p, address, basis):
    k = 1 << (p - 1)
    u = basis & ((1 << k) - 1)
    return basis ^ (_linear_image(_pencil_columns(p, address), u) << k)


def _translation(p, address, basis):
    return sum(((basis >> wire) & 1) << (wire ^ address)
               for wire in range(1 << p))


def _permutation_matrix(permutation):
    result = zeros(len(permutation))
    for column, row in enumerate(permutation):
        result[row][column] = CO
    return result


def _query_word(p):
    """All targets share one cubic batch, with private returned masks."""
    k, m = 1 << (p - 1), 1 << p
    triples = []
    for coordinate in range(p):
        for source, column in enumerate(_pencil_columns(p, 1 << coordinate)):
            for target in range(k):
                if column >> target & 1:
                    triples.append((coordinate, p + source, p + k + target))
    helpers = list(range(p + m, p + m + len(triples)))
    wall = [('H', p + k + bit) for bit in range(k)]
    core = _gray_phase_batch(triples, helpers, [1] * len(triples), 3)
    return wall + core + list(reversed(wall)), p + m + len(triples), triples


def _joint_eigenbasis(p, action):
    """Canonical orbit Fourier bases, grouped by their literal character."""
    m, dimension = 1 << p, 1 << (1 << p)
    sectors = [[] for _ in range(m)]
    unseen = set(range(dimension))
    while unseen:
        representative = min(unseen)
        points = {}
        stabilizer = []
        for address in range(m):
            point = action(p, address, representative)
            points.setdefault(point, address)
            if point == representative:
                stabilizer.append(address)
        unseen.difference_update(points)
        # Only p <= 2 is used: each orbit has size 1, 2, or 4.
        scale = {1: Q2(1), 2: Q2(0, Fraction(1, 2)),
                 4: Q2(Fraction(1, 2))}[len(points)]
        for character in range(m):
            if any((character & address).bit_count() % 2
                   for address in stabilizer):
                continue
            column = [CZ] * dimension
            for point, address in points.items():
                sign = (-1) ** ((character & address).bit_count() % 2)
                column[point] = CO * scale * sign
            sectors[character].append(column)
    columns = [column for sector in sectors for column in sector]
    return [list(row) for row in zip(*columns)], tuple(map(len, sectors))


def _p2_dense_intertwiner():
    source, source_counts = _joint_eigenbasis(2, _shear)
    target, target_counts = _joint_eigenbasis(2, _translation)
    assert source_counts == target_counts == (7, 3, 3, 3)
    return multiply(target, adjoint(source))


def _p1_native_source():
    # Address, u, v, arbitrary helper. E=CX(v;u) on either side.
    wall = [('CX', 2, 1)]
    return wall + _word(_compact_toffoli(0, 1, 2, 3)) + wall


def _native_x(wire):
    return [('H', wire), ('S', wire), ('S', wire), ('H', wire)]


def _native_swap(left, right):
    return [('CX', left, right), ('CX', right, left), ('CX', left, right)]


def _field_toffoli(a, b, target):
    return [('H', target)] + _shared_ccz(a, [(b, target)]) + [('H', target)]


def _controlled_companion(control, r, s, inverse=False):
    # C: (r,s)->(s,r+s), chronologically lower then upper shear.
    pair = [(r, s), (s, r)]
    if inverse:
        pair.reverse()
    return [gate for source, target in pair
            for gate in _field_toffoli(control, source, target)]


def _field_basis(d0=0, d1=1, r=2, s=3):
    # F_d=[d,Cd]=I,C,C^2 for d=1,2,3; F_0=I.
    # Its exponent is 2*d0+(d0 xor d1) modulo three.
    return (_controlled_companion(d0, r, s, inverse=True)
            + [('CX', d0, d1)] + _controlled_companion(d1, r, s)
            + [('CX', d0, d1)])


def _conditional_h_nonzero(d0=0, d1=1, r=2):
    # R=S H T H S^-1 satisfies R Z R^-1=H, with literal phase.
    rotation = [('S', r)] * 3 + [('H', r), ('T', r), ('H', r), ('S', r)]
    controlled_twice = (_inverse(rotation) + _shared_ccz(d0, [(d1, r)])
                        + rotation)
    # H * H^((1+d0)(1+d1)) = H^(d0 OR d1).
    negative = _native_x(d0) + _native_x(d1)
    return [('H', r)] + negative + controlled_twice + negative


def _p2_native_decoder():
    """Fixed four-qubit decoder, without initialized or borrowed wires."""
    # Input u0,u1,v0,v1; output physical wire labels 0,1,2,3.
    word = [('H', 2), ('H', 3)]
    word += _field_basis() + _native_swap(2, 3) + _field_toffoli(0, 1, 3)
    word += _conditional_h_nonzero()
    # a=r*u+s*C*u; b=a+u. Convert (u,a) to (a,b).
    word += _field_basis()
    word += [('CX', 2, 0), ('CX', 3, 1)]
    word += _native_swap(0, 2) + _native_swap(1, 3)
    # Bell transform on the two physical pairs.
    word += [('H', 0), ('CX', 0, 1), ('H', 2), ('CX', 2, 3)]
    return word


class FieldShearIntertwinerTests(unittest.TestCase):
    def test_nonzero_pencils_invertible_and_sparse_occurrences(self):
        for p in _MODULI:
            k, m = 1 << (p - 1), 1 << p
            # This exhausts the field too, so irreducibility of every
            # displayed modulus is checked without relying on its name.
            for element in range(1, 1 << k):
                columns = [_field_product(element, 1 << bit, _MODULI[p])
                           for bit in range(k)]
                self.assertEqual(_binary_rank(columns), k)
            occurrences = 0
            for coordinate in range(p):
                columns = _pencil_columns(p, 1 << coordinate)
                count = sum(column.bit_count() for column in columns)
                self.assertLessEqual(count, k + (k - 1) * coordinate)
                occurrences += count
            self.assertLessEqual(occurrences,
                                 k * p + (k - 1) * p * (p - 1) // 2)
            for address in range(1, m):
                self.assertEqual(_binary_rank(_pencil_columns(p, address)), k)

    def test_complete_group_characters_and_different_common_fixed_sets(self):
        for p in (1, 2, 3):
            k, m = 1 << (p - 1), 1 << p
            dimension = 1 << m
            characters = []
            for address in range(m):
                fixed = []
                for action in (_shear, _translation):
                    fixed.append(sum(action(p, address, basis) == basis
                                     for basis in range(dimension)))
                self.assertEqual(fixed[0], fixed[1])
                self.assertEqual(fixed[0], dimension if address == 0 else 1 << k)
                characters.append(fixed[0])
            multiplicities = [sum((-1) ** ((s & x).bit_count() % 2) * value
                                  for x, value in enumerate(characters)) // m
                              for s in range(m)]
            expected = [(dimension + (m - 1) * (1 << k)) // m]
            expected += [(dimension - (1 << k)) // m] * (m - 1)
            self.assertEqual(multiplicities, expected)
            common = [sum(all(action(p, x, basis) == basis for x in range(m))
                          for basis in range(dimension))
                      for action in (_shear, _translation)]
            self.assertEqual(common, [1 << k, 2])

    def test_all_native_phase_columns_and_actual_inverse_p2(self):
        # Evaluate the complete cubic core on every address, carrier and
        # helper column. H on all v wires conjugates this exact phase to Q.
        p, k = 2, 2
        word, width, triples = _query_word(p)
        for native in (word, _inverse(word)):
            self.assertEqual(native[:k], [('H', p + k + bit) for bit in range(k)])
            self.assertEqual(set(native[-k:]), set(native[:k]))
            core = native[k:-k]
            for basis in range(1 << width):
                cubic = sum(all(basis >> wire & 1 for wire in term)
                            for term in triples)
                self.assertEqual(_exact_phase_action(basis, core),
                                 (basis, 4 * (cubic % 2)))

    def test_native_shear_all_columns_p1_and_uniform_paid_ledger(self):
        word, width, _ = _query_word(1)
        expected = _permutation_matrix([
            column ^ ((_shear(1, column & 1, (column >> 1) & 3)
                       ^ ((column >> 1) & 3)) << 1)
            for column in range(1 << width)])
        self.assertEqual(native_columns(width, word), expected)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(expected))
        for p in _MODULI:
            word, width, triples = _query_word(p)
            k, m, occurrences = 1 << (p - 1), 1 << p, len(triples)
            self.assertEqual(width, p + m + occurrences)
            for native in (word, _inverse(word)):
                counts = Counter(gate[0] for gate in native)
                self.assertEqual(counts['T'] + counts['TDG'], 8 * occurrences)
                self.assertEqual(counts['CX'], 8 * occurrences)
                self.assertEqual(counts['H'], 2 * k)
                self.assertEqual(len(native), 16 * occurrences + 2 * k)
                self.assertEqual(len(_bands(native)), 8)
                for band in _bands(native):
                    self.assertEqual(len({gate[1] for gate in band}), len(band))

    def test_p2_exact_joint_eigenbasis_intertwiner_all_operators(self):
        source, counts = _joint_eigenbasis(2, _shear)
        target, target_counts = _joint_eigenbasis(2, _translation)
        self.assertEqual(counts, target_counts)
        self.assertEqual(counts, (7, 3, 3, 3))
        self.assertEqual(multiply(adjoint(source), source), eye(16))
        self.assertEqual(multiply(adjoint(target), target), eye(16))
        intertwiner = _p2_dense_intertwiner()
        self.assertEqual(multiply(adjoint(intertwiner), intertwiner), eye(16))
        x0 = _permutation_matrix([basis ^ 1 for basis in range(16)])
        for address in range(4):
            query = _permutation_matrix([_shear(2, address, basis)
                                         for basis in range(16)])
            translation = _permutation_matrix([_translation(2, address, basis)
                                               for basis in range(16)])
            actual = multiply(multiply(intertwiner, query), adjoint(intertwiner))
            self.assertEqual(actual, translation)
            self.assertEqual(multiply(multiply(actual, x0), adjoint(actual)),
                             _permutation_matrix([basis ^ (1 << address)
                                                  for basis in range(16)]))

    def test_clifford_similarity_obstruction_already_at_p2(self):
        for p in (2, 3):
            m = 1 << p
            # X-block deviations N_x=A_x-I of Q_x multiply to zero.
            for x, y in ((1, 2), (1, 3), (2, 3)):
                for bit in range(m):
                    vector = _shear(p, y, 1 << bit) ^ (1 << bit)
                    self.assertEqual(_shear(p, x, vector) ^ vector, 0)
                # W deviations have a nonzero four-point mixed product.
                vector = _translation(p, y, 1) ^ 1
                self.assertEqual(_translation(p, x, vector) ^ vector,
                                 1 ^ (1 << x) ^ (1 << y) ^ (1 << (x ^ y)))
                self.assertNotEqual(_translation(p, x, vector) ^ vector, 0)

    def test_p1_native_intertwiner_source_and_complete_indicator(self):
        source = _p1_native_source()
        expected = _permutation_matrix([
            column ^ ((_translation(1, column & 1, (column >> 1) & 3)
                       ^ ((column >> 1) & 3)) << 1)
            for column in range(16)])
        self.assertEqual(native_columns(4, source), expected)
        self.assertEqual(native_columns(4, _inverse(source)), adjoint(expected))
        probe = [('H', 1), ('S', 1), ('S', 1), ('H', 1)]
        indicator = _inverse(source) + probe + source
        target = _permutation_matrix([column ^ (1 << (1 + (column & 1)))
                                      for column in range(16)])
        self.assertEqual(native_columns(4, indicator), target)
        self.assertEqual(native_columns(4, _inverse(indicator)), adjoint(target))
        self.assertEqual(len(source), 28)
        self.assertEqual(len(_bands(source)), 2)
        counts = Counter(gate[0] for gate in source)
        self.assertEqual(counts, Counter({'CX': 18, 'H': 2, 'T': 4, 'TDG': 4}))

    def test_p2_native_decoder_all_operators_and_actual_inverse(self):
        word = _p2_native_decoder()
        actual = native_columns(4, word)
        reverse = native_columns(4, _inverse(word))
        self.assertEqual(reverse, adjoint(actual))
        self.assertEqual(multiply(reverse, actual), eye(16))
        x0 = _permutation_matrix([basis ^ 1 for basis in range(16)])
        for address in range(4):
            query = _permutation_matrix([_shear(2, address, basis)
                                         for basis in range(16)])
            translation = _permutation_matrix([_translation(2, address, basis)
                                               for basis in range(16)])
            source = multiply(multiply(actual, query), reverse)
            self.assertEqual(source, translation)
            self.assertEqual(multiply(multiply(source, x0), adjoint(source)),
                             _permutation_matrix([basis ^ (1 << address)
                                                  for basis in range(16)]))
        # The fixed pulled-back seed R=E^-1 X0 E is genuinely non-Clifford:
        # conjugating X0 by R produces two nonzero entries in column zero.
        seed = multiply(multiply(reverse, x0), actual)
        witness = multiply(multiply(seed, x0), seed)
        self.assertEqual([row[0] for row in witness],
                         [CO * Q2(0, Fraction(1, 2)) if row in (1, 5) else CZ
                          for row in range(16)])
        for native in (word, _inverse(word)):
            counts = Counter(gate[0] for gate in native)
            self.assertEqual(counts['T'] + counts['TDG'], 72)
            self.assertEqual(counts['CX'], 117)
            self.assertEqual(counts['H'], 35)
            self.assertEqual(counts['S'], 16)
            self.assertEqual(len(native), 240)
            self.assertEqual(len(_bands(native)), 42)
            self.assertEqual({wire for gate in native for wire in gate[1:]},
                             set(range(4)))
            for band in _bands(native):
                self.assertEqual(len({gate[1] for gate in band}), len(band))

    def test_p2_complete_native_source_and_indicator_ledgers(self):
        query, width, _ = _query_word(2)
        decoder = [(gate[0], *(wire + 2 for wire in gate[1:]))
                   for gate in _p2_native_decoder()]
        # Chronology E^-1,Q,E has operator E Q E^-1.  Full correctness
        # follows from the preceding 16x16 block identities and the exact
        # all-2048-column native cubic-core identity, with identity masks.
        source = _inverse(decoder) + query + decoder
        indicator = _inverse(source) + _native_x(2) + source
        self.assertEqual(width, 11)
        for word, t_count, total, bands in ((source, 184, 564, 92),
                                           (indicator, 368, 1132, 184)):
            for native in (word, _inverse(word)):
                counts = Counter(gate[0] for gate in native)
                self.assertEqual(counts['T'] + counts['TDG'], t_count)
                self.assertEqual(len(native), total)
                self.assertLessEqual(len(_bands(native)), bands)
                self.assertTrue(all(0 <= wire < width
                                    for gate in native for wire in gate[1:]))
                for band in _bands(native):
                    self.assertEqual(len({gate[1] for gate in band}), len(band))


if __name__ == '__main__':
    unittest.main()
