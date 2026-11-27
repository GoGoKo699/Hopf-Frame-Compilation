"""Exact small checks of Pauli-transfer normalization and the precision witness."""

from collections import Counter
from fractions import Fraction
from itertools import product
import unittest

from tests.exact_matrices import Q2, CZ, CO, add, multiply, native_columns
from tests.test_two_layer_galois_return import _conjugate, _sigma, _tau


def _normalized_trace_product(left, right):
    return sum((a * b for row, column in zip(left, zip(*right))
                for a, b in zip(row, column)), CZ) / len(left)


def _valuation(value):
    """The 2-adic valuation of the norm equals the sqrt(2)-valuation."""
    norm = abs(value.a**2 - 2 * value.b**2)
    assert norm
    numerator, denominator = norm.numerator, norm.denominator
    return ((numerator & -numerator).bit_length()
            - (denominator & -denominator).bit_length())


class PrecisionCountDepthTradeoffTests(unittest.TestCase):
    def test_scaled_single_band_diagonalizes_over_the_integral_ring(self):
        paulis = [_tau(label, 1) for label in (0, 2, 1, 3)]
        root = Q2(0, 1)
        for name, orientation in (("T", 1), ("TDG", -1)):
            unitary = native_columns(1, [(name, 0)])
            images = [_conjugate(unitary, pauli) for pauli in paulis]
            band = [[root * _normalized_trace_product(a, b).real
                     for b in images] for a in paulis]
            band[3] = [a - orientation * b for a, b in zip(band[3], band[2])]
            for row in band:
                row[3] += orientation * row[2]
            diagonal = [root, root, Q2(1), Q2(2)]
            self.assertEqual(band, [[value if i == j else Q2()
                                     for j in range(4)]
                                    for i, value in enumerate(diagonal)])
            exponents = [_valuation(value) for value in diagonal]
            self.assertEqual(Counter(exponents), {0: 1, 1: 2, 2: 1})
            self.assertEqual(Counter(a + b for a, b in product(exponents, repeat=2)),
                             {0: 1, 1: 4, 2: 6, 3: 4, 4: 1})

    def test_initialized_transfer_block_includes_a_leaked_dirty_column(self):
        # Target 0 rotates, and dirty spectator 1 is copied into clean 2.
        word = [("H", 0), ("T", 0), ("H", 0), ("CX", 1, 2)]
        unitary = native_columns(3, word)
        system = [_tau(label, 2) for label in (4, 6, 14, 12)]
        physical_labels = (8, 10, 26, 24)
        target_paulis = [_tau(label, 3) for label in physical_labels]
        trace_block = []
        copied_block = []
        for label in physical_labels:
            code_probe = add(_tau(label, 3), _tau(label ^ 32, 3))
            full_image = _conjugate(unitary, code_probe)
            trace_block.append([_normalized_trace_product(pauli, full_image)
                                for pauli in target_paulis])
        kraus = [[row[:4] for row in unitary[4 * outcome:4 * (outcome + 1)]]
                 for outcome in (0, 1)]
        for pauli in system:
            image = add(*(_conjugate(block, pauli) for block in kraus))
            copied_block.append([_normalized_trace_product(out, image) for out in system])
        self.assertEqual(trace_block, copied_block)
        transfer = Q2(0, Fraction(1, 2))
        expected = [[CO * transfer if i == j and i in (0, 3) else CZ
                     for j in range(4)] for i in range(4)]
        self.assertEqual(trace_block, expected)
        self.assertEqual(_sigma(trace_block), [[-item for item in row] for row in expected])

        # The two nonzero diagonal positions contribute the nontrivial minor.
        det = (3 * transfer - 1) * (3 * transfer - 1)
        self.assertEqual(_valuation(det), -2)
        self.assertEqual(det * Q2(det.a, -det.b), Fraction(49, 4))
