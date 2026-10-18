"""Bounded evidence for the abstract partial-normalizer counterexample.

The growing-width degree statement is analytic. Here every work column
is checked for small banks, and the charged Fredkin source is expanded
into the retained literal seven-T, five-layer CCZ construction.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_read_only_three_layer_rigidity import H, S, _degree
from tests.test_returned_diagonal_echo_degree import _ccz, _column, _inverse


def _bank_translation(column, width, shift):
    address_mask = (1 << width) - 1
    result = column & address_mask
    for source in range(1 << width):
        if column >> (width + source) & 1:
            result |= 1 << (width + (source ^ shift))
    return result


def _fredkin(column, control, left, right):
    if column >> control & 1 and ((column >> left) ^ (column >> right)) & 1:
        column ^= (1 << left) | (1 << right)
    return column


def _source(column, width):
    for coordinate in range(width):
        for label in range(1 << width):
            if not (label >> coordinate) & 1:
                column = _fredkin(column, coordinate, width + label,
                                  width + (label ^ (1 << coordinate)))
    return column


def _cx(control, target):
    matrix = np.zeros((8, 8), dtype=complex)
    for column in range(8):
        image = column ^ ((1 << target) if (column >> control) & 1 else 0)
        matrix[image, column] = 1
    return matrix


def _h(target):
    matrix = np.zeros((8, 8), dtype=complex)
    for column in range(8):
        for bit in range(2):
            row = (column & ~(1 << target)) | (bit << target)
            matrix[row, column] = H[bit, (column >> target) & 1]
    return matrix


class PartialNormalizerScopeTests(unittest.TestCase):
    def test_full_work_word_and_full_tableau_degree(self):
        for width in range(1, 4):
            work = width + (1 << width)
            address_mask = (1 << width) - 1
            images = []
            for column in range(1 << work):
                actual = _source(column, width)
                expected = _bank_translation(column, width, column & address_mask)
                self.assertEqual(actual, expected)
                self.assertEqual(_source(actual, width), column)
                images.append(actual)
            self.assertEqual(len(set(images)), 1 << work)
            # The extra arbitrary spectator is flipped by B=U tensor X.
            def balanced(column):
                spectator = column & (1 << work)
                return _source(column & ((1 << work) - 1), width) | (
                    spectator ^ (1 << work))

            for column in range(1 << (work + 1)):
                for x in range(1 << width):
                    image = balanced(balanced(column) ^ x)
                    spectator = column & (1 << work)
                    expected = _bank_translation(
                        (column & ((1 << work) - 1)) ^ x, width, x) | spectator
                    self.assertEqual(image, expected)
            # Input X_z^x has fixed identity symplectic tableau and linear
            # varying signs. The image's data-wire permutation occurs in
            # both symplectic blocks; these are FULL-tableau coordinates.
            for coordinate in range(width):
                self.assertEqual(_degree([(x >> coordinate) & 1
                                          for x in range(1 << width)]), 1)
            for source in range(1 << width):
                for target in range(1 << width):
                    entry = [int(target == source ^ x)
                             for x in range(1 << width)]
                    self.assertEqual(_degree(entry), width)

    def test_native_fredkin_actual_inverse_and_source_ledger(self):
        schedule = _ccz(0, 1, 2)
        layers = [gates for kind, gates in schedule if kind == 'T']
        self.assertEqual((sum(map(len, layers)), len(layers)), (7, 5))
        cliffords = sum(len(gates) for kind, gates in schedule if kind == 'C')
        self.assertEqual(cliffords, 10)
        for layer in layers:
            self.assertEqual(len({gate[1] for gate in layer}), len(layer))
        diagonal = []
        for column in range(8):
            self.assertEqual(_column(schedule, column),
                             (column, 4 * int(column == 7)))
            self.assertEqual(_column(schedule + _inverse(schedule), column),
                             (column, 0))
            diagonal.append(-1 if column == 7 else 1)
        ccx = _h(2) @ np.diag(diagonal) @ _h(2)
        actual = _cx(2, 1) @ ccx @ _cx(2, 1)
        expected = np.zeros((8, 8), dtype=complex)
        for column in range(8):
            expected[_fredkin(column, 0, 1, 2), column] = 1
        np.testing.assert_allclose(actual, expected, atol=2e-15)
        np.testing.assert_allclose(actual.conj().T @ actual, np.eye(8), atol=2e-15)
        # Two external CNOTs and two Toffoli H gates are charged too.
        self.assertEqual(cliffords + 2 + 2, 14)
        np.testing.assert_allclose(H @ S @ S @ H,
                                   np.array([[0, 1], [1, 0]]), atol=2e-15)
        for width in range(1, 5):
            swaps = width * (1 << (width - 1))
            self.assertEqual(sum(1 for coordinate in range(width)
                                 for label in range(1 << width)
                                 if not (label >> coordinate) & 1), swaps)
            ledger = (7 * swaps, 14 * swaps + 4, 5 * swaps,
                      width + (1 << width) + 1, 0)
            self.assertEqual(ledger, (7 * width * 2 ** (width - 1),
                                      14 * width * 2 ** (width - 1) + 4,
                                      5 * width * 2 ** (width - 1),
                                      width + 2 ** width + 1, 0))


if __name__ == '__main__':
    unittest.main()
