"""Exact local ingredients of the unrestricted three-band source obstruction.

The all-size rank and key-rigidity proof is analytic. These bounded complete
matrices check dressed point defects, literal band normalization, and the
independence/Clifford hypotheses; they do not search for a source or simulate
the exponentially wide p=7 target.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, Q2, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_joint_lookup_pauli_source import _inverse


def _chain(*operators):
    result = eye(len(operators[0]))
    for operator in operators:
        result = multiply(result, operator)
    return result


def _conjugate(unitary, operator):
    return _chain(unitary, operator, adjoint(unitary))


def _difference(left, right):
    return add(left, scale(right, -1))


def _rank(matrix):
    """Eliminate in Q(sqrt(2),i), without numerical rank thresholds."""
    rows = [list(row) for row in matrix]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows))
                      if rows[i][column] != CZ), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        divisor = rows[rank][column]
        rows[rank] = [entry / divisor for entry in rows[rank]]
        for i in range(rank + 1, len(rows)):
            factor = rows[i][column]
            if factor != CZ:
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def _pauli(width, x, z):
    matrix = zeros(1 << width)
    scalar = (CO, CI, -CO, -CI)[(x & z).bit_count() % 4]
    for column in range(1 << width):
        matrix[column ^ x][column] = scalar * (-1) ** ((column & z).bit_count() % 2)
    return matrix


def _support(matrix, width):
    dimension = 1 << width
    support = {}
    for x in range(dimension):
        for z in range(dimension):
            pauli = _pauli(width, x, z)
            coefficient = sum((pauli[row ^ x][row].conj()
                               * matrix[row ^ x][row]
                               for row in range(dimension)), CZ) / dimension
            if coefficient != CZ:
                support[x | (z << width)] = coefficient
    return support


def _pair(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((right & mask) & (left >> width)).bit_count()) % 2


def _diagonal(values):
    result = zeros(len(values))
    for index, value in enumerate(values):
        result[index][index] = CO * value
    return result


def _lookup(p):
    width = p + (1 << p)
    signs = []
    for basis in range(1 << width):
        address = basis & ((1 << p) - 1)
        coefficient = basis >> p
        exponent = sum((coefficient >> subset) & 1
                       for subset in range(1 << p) if subset & address == subset)
        signs.append((-1) ** (exponent % 2))
    return _diagonal(signs)


class ThreeBandLookupSourceObstructionTests(unittest.TestCase):
    def assertHermitianClifford(self, matrix, width):
        self.assertEqual(matrix, adjoint(matrix))
        self.assertEqual(multiply(matrix, matrix), eye(1 << width))
        for bit in range(2 * width):
            mask = 1 << bit
            image = _conjugate(matrix, _pauli(width, mask & ((1 << width) - 1),
                                             mask >> width))
            support = _support(image, width)
            self.assertEqual(len(support), 1)
            self.assertIn(next(iter(support.values())), (CO, -CO))

    def test_dressed_point_rank_identity_with_address_bank_mixing(self):
        p, width = 2, 6
        dimension = 1 << width
        # C Z_0 C† = X_0 Z_5, coupling address and the top coefficient.
        frame_word = [('H', 0), ('H', 5), ('CX', 0, 5), ('H', 5)]
        word = _inverse(frame_word) + [('T', 0)] + frame_word
        g = native_columns(width, word)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(g))
        target = _lookup(p)
        m = _conjugate(g, target)
        probe = _pauli(width, 0, 1)
        ranks = []
        for address in range(1 << p):
            # ANF of [x=address]: product of x_i or (1+x_i).
            mask = sum(1 << subset for subset in range(1 << p)
                       if subset & address == address)
            bank_shift = _pauli(width, mask << p, 0)
            a = _conjugate(g, bank_shift)
            s = _conjugate(g, _diagonal([
                -1 if basis & 3 == address else 1 for basis in range(dimension)]))
            self.assertEqual(_conjugate(a, m), multiply(m, s))
            self.assertEqual(multiply(a, s), multiply(s, a))
            self.assertEqual(multiply(m, s), multiply(s, m))
            self.assertEqual(_rank(_difference(eye(dimension), s)), dimension // 4)
            shifted_probe = _conjugate(a, probe)
            left = _conjugate(a, _conjugate(m, probe))
            right = _conjugate(m, shifted_probe)
            expected = _conjugate(m, _difference(
                _conjugate(s, shifted_probe), shifted_probe))
            self.assertEqual(_difference(left, right), expected)
            rank = _rank(expected)
            self.assertLessEqual(rank, 2 * dimension // 4)
            ranks.append(rank)
        # This checks a genuine defect, not only the commuting zero case.
        self.assertTrue(any(ranks))

    def test_one_phase_trace_rank_floor_and_strict_address_threshold(self):
        width, dimension = 3, 8  # Two spectator qubits remain arbitrary.
        h = native_columns(width, [('H', 0)])
        t = native_columns(width, [('T', 0)])
        right = _conjugate(t, h)
        overlap = sum(multiply(h, right)[i][i] for i in range(dimension)) / dimension
        delta = Q2(Fraction(1, 2), Fraction(-1, 4))
        self.assertEqual(overlap, CO * (1 - delta))
        defect = _difference(h, right)
        hs_squared = sum(entry.abs2() for row in defect for entry in row)
        rank = _rank(defect)
        self.assertEqual(hs_squared, 2 * dimension * delta)
        self.assertLessEqual(hs_squared, 4 * rank)
        self.assertGreaterEqual(Q2(rank), dimension * delta / 2)
        z = _pauli(width, 0, 1)
        self.assertEqual(_rank(_difference(z, _conjugate(t, z))), 0)
        self.assertLess(Q2(Fraction(2, 32)), delta / 2)
        self.assertGreater(Q2(Fraction(2, 16)), delta / 2)
        self.assertEqual(Q2(4) / delta, Q2(16, 8))

    def test_independent_band_has_anticommuting_support_or_literal_fixation(self):
        width = 3
        # C maps Z_1,Z_2 to X_0 Z_1,X_0 Z_2, independently.
        frame = [('H', 0), ('CX', 0, 1), ('CX', 0, 2), ('H', 0)]
        word = _inverse(frame) + [('T', 1), ('TDG', 2)] + frame
        band = native_columns(width, word)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(band))
        axes = [1 | (2 << width), 1 | (4 << width)]
        self.assertEqual(_pair(*axes, width), 0)
        core = 1 << width
        support = _support(_conjugate(band, _pauli(width, 0, 1)), width)
        self.assertEqual(set(support),
                         {core, core ^ axes[0], core ^ axes[1], core ^ axes[0] ^ axes[1]})
        self.assertEqual(_pair(core, core ^ axes[0], width), 1)
        for z in (2, 4, 6):
            pauli = _pauli(width, 0, z)
            self.assertEqual(_conjugate(band, pauli), pauli)

    def test_absorbed_clifford_normal_form_retains_actual_first_prefix(self):
        width = 3
        words = [
            [('S', 0), ('CX', 0, 1)],
            [('H', 0), ('CX', 1, 2)],
            [('H', 2), ('S', 1), ('CX', 2, 0)],
            [('H', 1), ('CX', 0, 2), ('S', 2)],
        ]
        bands = [[('T', 0), ('TDG', 2)], [('TDG', 1)], [('T', 0), ('T', 2)]]
        c0, c1, c2, c3 = [native_columns(width, word) for word in words]
        d1, d2, d3 = [native_columns(width, word) for word in bands]
        complete_word = words[0] + bands[0] + words[1] + bands[1] + words[2] + bands[2] + words[3]
        complete = native_columns(width, complete_word)
        self.assertEqual(native_columns(width, _inverse(complete_word)), adjoint(complete))
        e3 = _conjugate(c3, d3)
        e2 = _conjugate(multiply(c3, c2), d2)
        prefix = _chain(c3, c2, c1, d1, c0)
        self.assertEqual(complete, _chain(e3, e2, prefix))
        probe = _pauli(width, 1, 5)
        k = _conjugate(prefix, probe)
        self.assertHermitianClifford(k, width)
        self.assertEqual(_conjugate(complete, probe),
                         _chain(e3, e2, k, adjoint(e2), adjoint(e3)))

    def test_key_fixing_band_can_have_nonclifford_work_fibers(self):
        width = 3
        # Key is wire0. Axes Z0,Z0X1,X2 include the key and commute with it.
        frame = [('CX', 0, 1), ('H', 1), ('H', 2)]
        word = _inverse(frame) + [('T', 0), ('T', 1), ('TDG', 2)] + frame
        band = native_columns(width, word)
        key = _pauli(width, 0, 1)
        self.assertEqual(_conjugate(band, key), key)
        work_frame = native_columns(width, [('H', 1), ('H', 2)])
        diagonal = _conjugate(adjoint(work_frame), band)
        self.assertTrue(all(entry == CZ for i, row in enumerate(diagonal)
                            for j, entry in enumerate(row) if i != j))
        zero_fiber = [[band[2 * i][2 * j] for j in range(4)] for i in range(4)]
        image = _conjugate(zero_fiber, _pauli(2, 0, 1))
        self.assertGreater(len(_support(image, 2)), 1)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(band))


if __name__ == '__main__':
    unittest.main()
