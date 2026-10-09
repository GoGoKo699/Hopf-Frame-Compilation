"""Exact finite checks of the tree phase gauge and regular Cayley chart.

All columns and literal phases are checked in Q(sqrt(2), i), including
singular original angles. These algebraic identities do not implement the
remaining Cayley transform as a native circuit.
"""
from fractions import Fraction as F
from itertools import product
import unittest

from tests.exact_matrices import (
    CQ2, Q2, CZ, CO, CI, add, adjoint, eye, multiply, scale, zeros,
)


RATIONAL_PHASES = (CQ2(F(3, 5), F(4, 5)), CQ2(F(5, 13), -F(12, 13)),
                   CQ2(-F(7, 25), F(24, 25)), CQ2(-F(20, 29), -F(21, 29)))
SINGULAR_PHASES = (CO, -CO, CI, -CI, RATIONAL_PHASES[0])


def _edges(n):
    return [(p << (n - d), (2 * p + 1) << (n - d - 1))
            for d in range(n) for p in range(1 << d)]


def _partial_word(n, phases):
    matrix = eye(1 << n)
    for (a, b), z in zip(_edges(n), phases):
        first, second = matrix[a][:], matrix[b][:]
        diagonal, off = (1 + z) / 2, (1 - z) / 2
        matrix[a] = [diagonal * x + off * y for x, y in zip(first, second)]
        matrix[b] = [off * x + diagonal * y for x, y in zip(first, second)]
    return matrix


def _coarse_permutation(n, bits):
    permutation = list(range(1 << n))
    for (a, b), bit in zip(_edges(n), bits):
        if bit:
            permutation = [b if x == a else a if x == b else x for x in permutation]
    return permutation


def _chart_coordinates(n, bits, phases):
    sigma, residual, depth_maps = [0], [], []
    for d in range(n):
        prefix = _coarse_permutation(d, bits[:(1 << d) - 1])
        mapping = [prefix[sigma[p]] for p in range(1 << d)]
        start = (1 << d) - 1
        residual.extend([phases[start + p] * (-1) ** bits[start + p] for p in mapping])
        next_sigma = [0] * (2 << d)
        for p, old_parent in enumerate(mapping):
            next_sigma[2 * p] = 2 * sigma[p]
            next_sigma[2 * p + 1] = 2 * old_parent + 1
        sigma = next_sigma
        depth_maps.append(mapping)
    return sigma, residual, depth_maps


def _strict_chart_bit(real_part):
    # A fixture enclosure oracle, with exact boundary values included.
    for precision in range(1, 16):
        radius = Q2(F(1, 1 << precision))
        if real_part + radius < F(1, 2):
            return 0, precision
        if real_part - radius > -F(1, 2):
            return 1, precision
    raise AssertionError('the overlapping strict chart tests must terminate')


class TreePhaseChartTests(unittest.TestCase):
    def test_phase_gauge_all_columns_fresh_markers_and_decrement(self):
        cases = (RATIONAL_PHASES, SINGULAR_PHASES, (CO,), (CI,))
        for n in range(1, 5):
            size = 1 << n
            for phases in cases:
                with self.subTest(n=n, phases=phases):
                    target = eye(size)
                    initial, current = [CO] * size, [CO] * size
                    seen = {0}
                    units = [phases[j % len(phases)] for j in range(size - 1)]
                    for (a, b), unit in zip(_edges(n), units):
                        self.assertEqual(unit.abs2(), 1)
                        self.assertIn(a, seen)
                        self.assertNotIn(b, seen)
                        seen.add(b)
                        initial[b] = current[b] = CI * current[a]
                        current[a] *= unit.conj()
                        current[b] *= unit.conj()
                        first, second = target[a][:], target[b][:]
                        target[a] = [CQ2(unit.real) * x - CQ2(unit.imag) * y
                                     for x, y in zip(first, second)]
                        target[b] = [CQ2(unit.imag) * x + CQ2(unit.real) * y
                                     for x, y in zip(first, second)]
                    partial = _partial_word(n, [unit * unit for unit in units])
                    gauged = [[current[i] * partial[i][j] * initial[j].conj()
                               for j in range(size)] for i in range(size)]
                    self.assertEqual(gauged, target)
                    self.assertEqual(seen, set(range(size)))
                    self.assertTrue(all(sum(row, CZ) == CO for row in partial))
                    self.assertEqual(multiply(partial, adjoint(partial)), eye(size))
                    if phases == (CI,):
                        self.assertEqual(partial, [[CQ2(i == (j - 1) % size)
                                                  for j in range(size)] for i in range(size)])

    def test_polynomial_resolvent_order_including_singular_phases(self):
        for n in (2, 3, 4):
            size, count = 1 << n, (1 << n) - 1
            incidence = zeros(size, count)
            for j, (a, b) in enumerate(_edges(n)):
                incidence[a][j], incidence[b][j] = CO, -CO
            gram = scale(multiply(adjoint(incidence), incidence), F(1, 2))
            lower = [[gram[i][j] if i > j else CZ for j in range(count)]
                     for i in range(count)]
            for values in (RATIONAL_PHASES, SINGULAR_PHASES):
                phases = [values[j % len(values)] for j in range(count)]
                delta = zeros(count)
                for j, z in enumerate(phases):
                    delta[j][j] = z - 1
                dl = multiply(delta, lower)
                coefficient = [row[:] for row in delta]
                # Forward substitution is the finite nilpotent polynomial;
                # it remains defined at z=1 without a division/equality test.
                for i in range(count):
                    for j in range(count):
                        for k in range(i):
                            if dl[i][k] != CZ and coefficient[k][j] != CZ:
                                coefficient[i][j] += dl[i][k] * coefficient[k][j]
                self.assertEqual(coefficient, add(delta, multiply(dl, coefficient)))
                represented = add(eye(size), scale(multiply(multiply(incidence, coefficient),
                                                            adjoint(incidence)), F(1, 2)))
                self.assertEqual(represented, _partial_word(n, phases))
                optical_identity = add(add(coefficient, adjoint(coefficient)),
                                       multiply(multiply(coefficient, gram), adjoint(coefficient)))
                self.assertEqual(optical_identity, zeros(count))
                if n == 3 and values == RATIONAL_PHASES:
                    # Transposing the chronological coefficient produces a
                    # different word; a reversed product cannot pass silently.
                    wrong = add(eye(size), scale(multiply(multiply(incidence, list(map(list, zip(*coefficient)))),
                                                         adjoint(incidence)), F(1, 2)))
                    self.assertNotEqual(wrong, represented)

    def test_every_small_coarse_pattern_and_individual_edge_conjugacy(self):
        values = RATIONAL_PHASES + SINGULAR_PHASES
        for n in range(1, 5):
            size, count = 1 << n, (1 << n) - 1
            phases = [values[j % len(values)] for j in range(count)]
            original = _partial_word(n, phases)
            patterns = product(range(2), repeat=count) if n <= 3 else (
                [0] * count, [1] * count,
                [j % 2 for j in range(count)],
                [int((j * j + 3 * j + 7) % 5 < 2) for j in range(count)],
            )
            for bits in patterns:
                sigma, residual, depth_maps = _chart_coordinates(n, bits, phases)
                permutation = _coarse_permutation(n, bits)
                represented = zeros(size)
                partial = _partial_word(n, residual)
                for i in range(size):
                    for j in range(size):
                        represented[permutation[sigma[i]]][sigma[j]] = partial[i][j]
                self.assertEqual(represented, original)
                prefix = list(range(size))
                transformed = []
                for (a, b), bit in zip(_edges(n), bits):
                    inverse = [0] * size
                    for i, j in enumerate(prefix):
                        inverse[j] = i
                    transformed.append((inverse[a], inverse[b]))
                    if bit:
                        prefix = [b if x == a else a if x == b else x for x in prefix]
                for d, mapping in enumerate(depth_maps):
                    start = (1 << d) - 1
                    for p, mapped in enumerate(mapping):
                        a, b = _edges(n)[start + p]
                        self.assertEqual((sigma[a], sigma[b]), transformed[start + mapped])

    def test_strict_chart_selection_and_bounded_cotangent(self):
        for real in (Q2(-1), Q2(-F(1, 2)), Q2(0), Q2(F(1, 2)), Q2(1)):
            bit, steps = _strict_chart_bit(real)
            self.assertLess(steps, 16)
            self.assertLess((-1) ** bit * real, F(1, 2))
        for phase in RATIONAL_PHASES + SINGULAR_PHASES:
            bit, _ = _strict_chart_bit(phase.real)
            residual = (-1) ** bit * phase
            self.assertLess(residual.real, F(1, 2))
            self.assertGreater((residual - 1).abs2(), 1)
            cotangent = CI * (residual + 1) / (residual - 1)
            self.assertEqual(cotangent.imag, 0)
            self.assertLess(cotangent.real * cotangent.real, 3)
            self.assertEqual((cotangent + CI) / (cotangent - CI), residual)

    def test_tree_cut_congruence_and_haar_diagonal_difference(self):
        for n in range(1, 5):
            size = 1 << n
            basis, support, intervals = zeros(size), [], []
            root = CQ2(Q2(0, F(1, 2)))
            uniform = CO
            for _ in range(n):
                uniform *= root
            for x in range(size):
                basis[x][0] = uniform
            for a, b in _edges(n):
                length = b - a
                support.append(set(range(b, b + length)))
                intervals.append(set(range(a, b + length)))
                coefficient = CO
                for _ in range((2 * length).bit_length() - 1):
                    coefficient *= root
                for x in range(a, b):
                    basis[x][b] = -coefficient
                for x in range(b, b + length):
                    basis[x][b] = coefficient
            self.assertEqual(multiply(adjoint(basis), basis), eye(size))
            projection = [[CQ2(int(i == j) - F(1, size)) for j in range(size)] for i in range(size)]
            incidence, cut_inverse = zeros(size, size - 1), zeros(size - 1, size)
            for j, ((a, b), subtree) in enumerate(zip(_edges(n), support)):
                incidence[a][j], incidence[b][j] = root, -root
                for x in range(size):
                    cut_inverse[j][x] = -2 * root * (F(int(x in subtree)) - F(len(subtree), size))
            self.assertEqual(multiply(incidence, cut_inverse), projection)
            self.assertEqual(multiply(cut_inverse, incidence), eye(size - 1))
            left, centered = zeros(size), zeros(size)
            diagonal_a, diagonal_b = [F(0)] * size, [F(0)] * size
            for j, subtree in enumerate(support):
                coefficient = F((7 * j) % 9 - 4, 5)
                for x in subtree:
                    diagonal_a[x] += 2 * coefficient * len(subtree)
                    for y in subtree:
                        left[x][y] += CQ2(2 * coefficient)
                for (_, b), interval in zip(_edges(n), intervals):
                    if interval <= subtree:
                        diagonal_b[b] += 2 * coefficient * len(subtree)
                centered_cut = [F(int(x in subtree)) - F(len(subtree), size)
                                for x in range(size)]
                for x in range(size):
                    for y in range(size):
                        centered[x][y] += CQ2(2 * coefficient * centered_cut[x] * centered_cut[y])
            first, second = zeros(size), zeros(size)
            for j in range(size):
                first[j][j], second[j][j] = CQ2(diagonal_a[j]), CQ2(diagonal_b[j])
            self.assertEqual(left, add(first, scale(multiply(multiply(basis, second), adjoint(basis)), -1)))
            self.assertEqual(centered, multiply(multiply(projection, left), projection))
            self.assertEqual(multiply(centered, [[CO] for _ in range(size)]), zeros(size, 1))


if __name__ == '__main__':
    unittest.main()
