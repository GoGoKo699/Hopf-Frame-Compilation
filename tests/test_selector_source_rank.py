"""Exact fixtures for scoped selector-source rank and uniform-vector screens.

The rank argument concerns complete origin/empty-bank fibers and common
interlayers that preserve those fibers.  These tests include arbitrary
selector inputs, a dirty private wire, and a clean signal in a small coherent
example.  They do not establish a general native-depth obstruction.  In
particular, the rank test is inconclusive for three p=3 stage calls.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.test_co_moving_lookup_stage import _Matching


def _row_swaps(word, swaps):
    """Execute physical row swaps, independently of Matching.relabel."""
    for first, second in swaps:
        different = ((word >> first) ^ (word >> second)) & 1
        word ^= different * ((1 << first) | (1 << second))
    return word


def _origin_permutation(matching):
    return tuple(_row_swaps(s, matching.swaps) if s & 1 else s
                 for s in range(1 << matching.m))


def _involution_defect_rank(permutation):
    """Count antisymmetric two-cycles, rather than forming a dense matrix."""
    if any(permutation[image] != source
           for source, image in enumerate(permutation)):
        raise AssertionError("the permutation is not an involution")
    return sum(source < image for source, image in enumerate(permutation))


def _hadamards(vector, wires):
    """Exact even-sized Hadamard wall by a fast butterfly transform."""
    wires = tuple(wires)
    if len(wires) % 2:
        raise ValueError("the rational fixture uses even Hadamard walls")
    result = list(vector)
    for wire in wires:
        bit = 1 << wire
        for first in range(len(result)):
            if first & bit:
                continue
            second = first | bit
            a, b = result[first], result[second]
            result[first], result[second] = a + b, a - b
    denominator = 1 << (len(wires) // 2)
    return [value / denominator for value in result]


def _permute(vector, image):
    result = [Fraction(0)] * len(vector)
    for source, amplitude in enumerate(vector):
        result[image(source)] = amplitude
    return result


def _cell_ccz(vector, m, omitted):
    """CCZ on three corners of every four-row cell in label bits 0,1."""
    masks = [sum(1 << (base + corner) for corner in range(4)
                 if corner != omitted) for base in range(0, m, 4)]
    return [value * (-1) ** sum((s & mask) == mask for mask in masks)
            for s, value in enumerate(vector)]


def _fourier_ccz(vector, m, omitted):
    return _hadamards(_cell_ccz(_hadamards(vector, range(m)), m, omitted),
                      range(m))


def _literal_origin_vector(vector, matching, bank=1, inverse=False):
    result = [Fraction(0)] * len(vector)
    for selectors, amplitude in enumerate(vector):
        x, y, output, borrowed, phase = matching.stage_literal(
            (0, bank, selectors, 1, 0), inverse=inverse)
        if (x, y, borrowed, phase) != (0, bank, 1, 0):
            raise AssertionError("literal stage left its phase-free fiber")
        result[output] += amplitude
    return result


def _rref_and_kernel(columns):
    """Exact rank and one nonzero null vector of a short rectangular map."""
    matrix = [list(row) for row in zip(*columns)]
    height, width = len(matrix), len(columns)
    pivots = []
    for column in range(width):
        pivot = next((row for row in range(len(pivots), height)
                      if matrix[row][column]), None)
        if pivot is None:
            continue
        row_index = len(pivots)
        matrix[row_index], matrix[pivot] = matrix[pivot], matrix[row_index]
        scale = matrix[row_index][column]
        matrix[row_index] = [value / scale for value in matrix[row_index]]
        for row in range(height):
            if row == row_index or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [value - factor * entry
                           for value, entry in zip(matrix[row], matrix[row_index])]
        pivots.append(column)
    free = next((column for column in range(width) if column not in pivots), None)
    if free is None:
        return len(pivots), None
    vector = [Fraction(0)] * width
    vector[free] = Fraction(1)
    for row, pivot in enumerate(pivots):
        vector[pivot] = -matrix[row][free]
    return len(pivots), vector


class SelectorSourceRankTests(unittest.TestCase):
    def test_literal_origin_empty_fibers_and_actual_inverse(self):
        for p, pairs in ((2, ((0, 1),)), (3, ((2, 0),)),
                         (4, ((0, 1), (2, 3)))):
            matching = _Matching(p, pairs)
            dimension = 1 << matching.m
            if p < 4:
                selectors = range(dimension)
            else:
                # Every p=4 selector column is counted below; literal phase
                # and borrowed-wire checks here use varied sparse/dense inputs.
                selectors = sorted({0, dimension - 1, 0x5555, 0xA55A}
                                   | {1 << row for row in range(matching.m)}
                                   | {(dimension - 1) ^ (1 << row)
                                      for row in range(matching.m)})
            for selector, bank, borrowed, phase, inverse in product(
                    selectors, (0, 1), (0, 1), (0, 3), (False, True)):
                state = (0, bank, selector, borrowed, phase)
                actual = matching.stage_literal(state, inverse=inverse)
                output = (_row_swaps(selector, matching.swaps)
                          if bank and selector & 1 else selector)
                self.assertEqual(actual, (0, bank, output, borrowed, phase))
                self.assertEqual(matching.stage_literal(actual, inverse=not inverse),
                                 state)

    def test_exact_permutation_cycle_ranks_include_full_matchings(self):
        for p, pairs, expected_rank in (
                (2, ((0, 1),), 2),
                (3, ((0, 1),), 48),
                (3, ((1, 2),), 48),
                (4, ((0, 1),), 15360),
                (4, ((0, 1), (2, 3)), 16128),
                (4, ((0, 2), (1, 3)), 16128)):
            matching = _Matching(p, pairs)
            q, m, dimension = len(pairs), matching.m, 1 << matching.m
            fixed_rows = sum(row == image for row, image in enumerate(matching.images))
            self.assertEqual(fixed_rows, m // (1 << q))
            k = (m - fixed_rows) // 2
            self.assertEqual(len(matching.swaps), k)
            actual_rank = _involution_defect_rank(_origin_permutation(matching))
            self.assertEqual(actual_rank, expected_rank)
            self.assertEqual(actual_rank, dimension // 4 - dimension // (1 << (k + 2)))

    def test_uniform_vector_survives_fourier_ccz_dressings(self):
        for p, families in (
                (2, (((0, 1),), ((1, 0),))),
                (3, (((0, 1),), ((1, 2),), ((2, 0),)))):
            m, dimension = 1 << p, 1 << (1 << p)
            plus = [Fraction(1, 1 << (m // 2))] * dimension
            zero = [Fraction(1)] + [Fraction(0)] * (dimension - 1)
            self.assertEqual(_hadamards(plus, range(m)), zero)
            for omitted in range(3):
                self.assertEqual(_fourier_ccz(plus, m, omitted), plus)
                # Mixing the two Fourier frames really changes these witnesses;
                # no source success is inferred from merely changing a vector.
                self.assertNotEqual(_cell_ccz(plus, m, omitted), plus)
                self.assertNotEqual(_fourier_ccz(zero, m, omitted), zero)
                for pairs, bank, inverse in product(families, (0, 1), (False, True)):
                    matching = _Matching(p, pairs)
                    dressed = _fourier_ccz(plus, m, omitted)
                    dressed = _literal_origin_vector(dressed, matching, bank, inverse)
                    dressed = _fourier_ccz(dressed, m, omitted)
                    self.assertEqual(dressed, plus)
            # Arbitrary complex path coefficients are a relaxation of native
            # signal coefficients.  These exact cases include saturation and
            # destructive interference, with a common output on the two fibers.
            coefficient_sets = (
                [(Fraction(1, 16), Fraction(0))] * 8,
                [(Fraction((-1) ** j, 8), Fraction(0)) for j in range(8)],
                [(Fraction(0), Fraction(1, 16))] * 8,
                [(Fraction(1, 8), Fraction((-1) ** j, 16)) for j in range(8)],
            )
            for coefficients in coefficient_sets:
                real = sum(pair[0] for pair in coefficients)
                imag = sum(pair[1] for pair in coefficients)
                output_norm = sum((real * value) ** 2 + (imag * value) ** 2
                                  for value in plus)
                self.assertEqual(output_norm, real * real + imag * imag)
                self.assertGreaterEqual(max(output_norm, abs(output_norm - Fraction(1, 2))),
                                        Fraction(1, 4))
            self.assertEqual(sum(pair[0] for pair in coefficient_sets[0]) ** 2,
                             Fraction(1, 4))

    def test_rank_budget_boundary_is_not_a_three_bit_exclusion(self):
        # The clean-code estimate is D - 2^k sum(r_j); a private dimension
        # multiplies both terms and cannot improve the sufficient criterion.
        for p, h, expected_margin in ((2, 2, 8), (2, 3, 4),
                                      (3, 2, 64), (3, 3, -32)):
            matching = _Matching(p, ((0, 1),))
            dimension = 1 << matching.m
            rank = _involution_defect_rank(_origin_permutation(matching))
            margin = dimension - 2 * h * rank  # one clean signal
            self.assertEqual(margin, expected_margin)
            for private_dimension in (1, 2, 5):
                self.assertEqual(dimension * private_dimension
                                 - 2 * private_dimension * h * rank,
                                 private_dimension * margin)
        self.assertLess(expected_margin, 0)  # inconclusive, never "excluded"

    def test_coherent_signal_private_interlayers_have_common_code_vector(self):
        # Wires: four selectors, initialized signal 4, arbitrary private wire 5.
        # The common interlayers entangle all three kinds of register.  This
        # fixture is outside a selector-only conjugation/branch polynomial.
        dimension = 64
        stages = [_origin_permutation(_Matching(2, pairs))
                  for pairs in (((0, 1),), ((1, 0),))]

        def cx(vector, control, target):
            return _permute(vector, lambda bits: bits ^ (((bits >> control) & 1) << target))

        def stage(vector, permutation):
            return _permute(vector, lambda bits: (bits & ~15) | permutation[bits & 15])

        def circuit(vector, marked):
            vector = _hadamards(vector, (4, 0))
            vector = cx(vector, 1, 5)
            if marked:
                vector = stage(vector, stages[0])
            vector = _hadamards(vector, (1, 5))
            vector = cx(vector, 2, 4)
            mask = (1 << 4) | (1 << 0) | (1 << 5)
            vector = [value * (-1 if bits & mask == mask else 1)
                      for bits, value in enumerate(vector)]
            if marked:
                vector = stage(vector, stages[1])
            vector = _hadamards(vector, (4, 5))
            return cx(vector, 5, 3)

        code_indices = [selectors | (private << 5)
                        for private, selectors in product(range(2), range(16))]
        marked_columns, empty_columns = [], []
        for index in code_indices:
            basis = [Fraction(0)] * dimension
            basis[index] = Fraction(1)
            marked_columns.append(circuit(basis, True))
            empty_columns.append(circuit(basis, False))
        differences = [[a - b for a, b in zip(marked, empty)]
                       for marked, empty in zip(marked_columns, empty_columns)]
        rank, code_vector = _rref_and_kernel(differences)
        self.assertGreater(rank, 0)
        self.assertLessEqual(rank, 16)  # 2^k W (r_C+r_D)=2*2*(2+2)
        self.assertIsNotNone(code_vector)
        initial = [Fraction(0)] * dimension
        for index, value in zip(code_indices, code_vector):
            initial[index] = value
        marked, empty = circuit(initial, True), circuit(initial, False)
        self.assertEqual(marked, empty)
        norm = sum(value * value for value in initial)
        self.assertGreater(norm, 0)
        self.assertEqual(sum(value * value for value in marked), norm)
        accepted = sum(value * value for bits, value in enumerate(marked)
                       if not (bits >> 4) & 1) / norm
        self.assertGreaterEqual(max(accepted, abs(accepted - Fraction(1, 2))),
                                Fraction(1, 4))


if __name__ == '__main__':
    unittest.main()
