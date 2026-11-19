"""Finite complete-operator and binary-pencil checks for prescribed scalar tables."""

from itertools import product
import unittest

import numpy as np

from compiler_robust_hopf.compiler_boundaries import normalized_hadamard
from compiler_robust_hopf.paid_logical_encodings import binary_rank


def subset_shear(address, bits):
    """I plus the edge from a subset with its highest set bit removed."""
    size = 1 << bits
    matrix = np.eye(size, dtype=int)
    for subset in range(1, size):
        highest = 1 << (subset.bit_length() - 1)
        matrix[subset, subset ^ highest] = int(bool(address & highest))
    return matrix


def marked_edge_graph(table, left_bits, right_bits, left, right):
    """Return the raw-table adjacency; column is the edge source."""
    left_size, right_size = table.shape
    right_root = 2 * left_size - 1
    size = right_root + 2 * right_size - 1
    edges = np.zeros((size, size), dtype=int)
    for offset, bits, address, towards_root in (
        (0, left_bits, left, True),
        (right_root, right_bits, right, False),
    ):
        for depth in range(bits):
            child_bit = (address >> (bits - depth - 1)) & 1
            for parent in range((1 << depth) - 1, (1 << (depth + 1)) - 1):
                child = 2 * parent + 1 + child_bit
                if towards_root:
                    edges[offset + parent, offset + child] = 1
                else:
                    edges[offset + child, offset + parent] = 1
    edges[left_size - 1:right_root, right_root + right_size - 1:] = table
    edges[right_root, 0] = 1
    return edges, right_root


def nullity(matrix):
    columns = [
        sum(int(bit) << site for site, bit in enumerate(column))
        for column in matrix.T
    ]
    return len(columns) - binary_rank(columns)


class PrescribedScalarTablesTests(unittest.TestCase):
    def test_shifted_nearest_even_tables_recover_all_two_bit_signs(self):
        # The nearest even-parity table to an odd table flips entry zero.
        # Translating the address first moves that single defect to the shift.
        def nearest_even_after_shift(table, shift):
            return table ^ ((table.bit_count() % 2) << shift)

        for target in range(16):
            lifted = sum(((target >> x) & 1) << (2 * x + 1) for x in range(4))
            derivatives = []
            for mask, shift in ((0b10100100, 0), (0b10100101, 3), (0b10001000, 5)):
                first = nearest_even_after_shift(mask, shift)
                second = nearest_even_after_shift(mask ^ lifted, shift)
                self.assertEqual((first.bit_count() % 2, second.bit_count() % 2), (0, 0))
                combined = first ^ second
                derivatives.append([
                    ((combined >> (2 * x)) ^ (combined >> (2 * x + 1))) & 1
                    for x in range(4)
                ])
            for x in range(4):
                self.assertEqual(sum(row[x] for row in derivatives) >= 2, (target >> x) & 1)

    def test_packed_walsh_characters_implement_each_full_xor(self):
        hadamard = normalized_hadamard(3)
        for selector, address in product(range(3), range(8)):
            input_sites = (selector, (selector + 1) % 3)
            target_sites = ((selector + 1) % 3, selector)
            control, data = ((address >> bit) & 1 for bit in input_sites)
            active_row = (control & data, control ^ data)
            word = sum(value << site for site, value in zip(target_sites, active_row))
            phases = []
            for character in range(8):
                dual = [(character >> bit) & 1 for bit in target_sites]
                phases.append((-1) ** sum(y * v for y, v in zip(dual, active_row)))
            actual = hadamard @ np.diag(phases) @ hadamard
            expected = np.eye(8)[:, [y ^ word for y in range(8)]]
            np.testing.assert_allclose(actual, expected, atol=1e-12)

    def test_coherent_derivative_majority_cancels_record_phases(self):
        # Each two-entry word is H D_A D_B H on one zero vote bit.
        hadamard = normalized_hadamard(1)
        records = []
        for left, right in product(range(4), repeat=2):
            signs = [(-1) ** (((left ^ right) >> t) & 1) for t in range(2)]
            recording = hadamard @ np.diag(signs) @ hadamard
            derivative = int(signs[0] != signs[1])
            records.append((recording, derivative))
        majority = np.diag([(-1) ** (v.bit_count() >= 2) for v in range(8)])
        for first, second, third in product(records, repeat=3):
            recording = np.kron(np.kron(first[0], second[0]), third[0])
            returned = recording.conj().T @ majority @ recording[:, 0]
            expected = np.zeros(8)
            expected[0] = (-1) ** (first[1] + second[1] + third[1] >= 2)
            np.testing.assert_allclose(returned, expected, atol=1e-12)

    def test_corank_and_clock_entry_for_every_small_table(self):
        for left_bits, right_bits in ((0, 0), (0, 1), (1, 1), (1, 2)):
            rows, cols = 1 << left_bits, 1 << right_bits
            for word in range(1 << (rows * cols)):
                table = np.array([(word >> i) & 1 for i in range(rows * cols)])
                table = table.reshape(rows, cols)
                # Each ANF coefficient is the parity on its subset downset.
                coefficients = np.zeros_like(table)
                for u, v in product(range(rows), range(cols)):
                    coefficients[u, v] = sum(
                        table[i, j] for i in range(rows) for j in range(cols)
                        if (i & u) == i and (j & v) == j
                    ) % 2
                for u, v in product(range(rows), range(cols)):
                    pencil = np.block([
                        [subset_shear(u, left_bits).T, coefficients],
                        [np.zeros((cols, rows), dtype=int), subset_shear(v, right_bits)],
                    ])
                    pencil[rows, 0] ^= 1
                    graph, root = marked_edge_graph(table, left_bits, right_bits, u, v)
                    power = np.linalg.matrix_power(graph, left_bits + right_bits + 2)
                    self.assertEqual(nullity(pencil), table[u, v])
                    self.assertEqual(nullity(np.eye(len(graph), dtype=int) ^ graph), table[u, v])
                    self.assertEqual(power[root, root], table[u, v])


if __name__ == "__main__":
    unittest.main()
