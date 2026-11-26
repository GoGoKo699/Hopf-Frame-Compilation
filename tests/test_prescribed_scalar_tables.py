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


def shear_table_completion(table):
    """The fixed binary shear whose upper-right block is the supplied table."""
    size = len(table)
    identity = np.eye(size, dtype=int)
    return np.block([[identity, table], [np.zeros_like(identity), identity]])


def bit_flip_matrix(width, bit):
    size = 1 << width
    return np.eye(size, dtype=int)[:, [basis ^ (1 << bit) for basis in range(size)]]


def four_call_scalar_difference(source, row_block_flip, column_block_flip):
    """Chronological V, X_a, V-dagger, X_b, V, X_a, V-dagger, X_b."""
    reverse = source.conj().T
    return (column_block_flip @ reverse @ row_block_flip @ source
            @ column_block_flip @ reverse @ row_block_flip @ source)


class PrescribedScalarTablesTests(unittest.TestCase):
    def test_shear_completion_and_four_calls_on_every_two_bit_table(self):
        # The source address is (alpha, u, beta, v), followed by one
        # arbitrary dirty spectator. Full matrices include every block-bit
        # and spectator input, rather than only the selected upper-right block.
        row_flip = bit_flip_matrix(5, 4)
        column_flip = bit_flip_matrix(5, 2)
        for word in range(16):
            table = np.array([(word >> bit) & 1 for bit in range(4)]).reshape(2, 2)
            with self.subTest(table=word):
                completed = shear_table_completion(table)
                np.testing.assert_array_equal(completed @ completed % 2, np.eye(4))
                nilpotent = completed ^ np.eye(4, dtype=int)
                np.testing.assert_array_equal(nilpotent @ nilpotent % 2,
                                              np.zeros((4, 4)))
                source = np.kron(np.diag(1 - 2 * completed.reshape(-1)), np.eye(2, dtype=int))
                actual = four_call_scalar_difference(source, row_flip, column_flip)
                expected = np.kron(np.diag(1 - 2 * np.tile(table, (2, 2)).reshape(-1)),
                                   np.eye(2, dtype=int))
                np.testing.assert_array_equal(actual, expected)

    def test_shear_completion_four_bit_and_odd_arity_padding(self):
        # These four-bit tables cover the zero, identity and dense nonsymmetric
        # cases. The final table duplicates an unused row-address bit of a
        # three-bit function; the complete output must return that pad too.
        tables = [
            np.zeros((4, 4), dtype=int),
            np.eye(4, dtype=int),
            np.array([(0xACE1 >> bit) & 1 for bit in range(16)]).reshape(4, 4),
        ]
        three_bit_table = np.array([(0xB4 >> bit) & 1 for bit in range(8)]).reshape(2, 4)
        tables.append(np.vstack((three_bit_table, three_bit_table)))
        row_flip = bit_flip_matrix(6, 5)
        column_flip = bit_flip_matrix(6, 2)
        for index, table in enumerate(tables):
            with self.subTest(table=index):
                completed = shear_table_completion(table)
                np.testing.assert_array_equal(completed @ completed % 2, np.eye(8))
                source = np.diag(1 - 2 * completed.reshape(-1))
                actual = four_call_scalar_difference(source, row_flip, column_flip)
                expected = np.diag(1 - 2 * np.tile(table, (2, 2)).reshape(-1))
                np.testing.assert_array_equal(actual, expected)
                if index == len(tables) - 1:
                    pad_flip = bit_flip_matrix(6, 4)
                    np.testing.assert_array_equal(pad_flip @ actual @ pad_flip, actual)

    def test_four_call_source_adjoint_with_address_motion_and_clean_leakage(self):
        # This is a source-interface fixture, not a claimed native emitter.
        # Physical order: alpha, u, beta, v, clean, arbitrary dirty work.
        # J includes every logical and dirty basis column. Its operator-norm
        # bounds also hold after tensoring any external reference.
        table = np.array([[0, 1], [1, 1]])
        completed = shear_table_completion(table)
        ideal_source = np.kron(np.diag(1 - 2 * completed.reshape(-1)), np.eye(4))
        target = np.kron(np.diag(1 - 2 * np.tile(table, (2, 2)).reshape(-1)), np.eye(4))
        identity = np.eye(64)
        initialization = identity[:, [basis for basis in range(64) if not (basis & 2)]]
        clean_projector = initialization @ initialization.T
        row_flip = bit_flip_matrix(6, 5)
        column_flip = bit_flip_matrix(6, 3)
        row_z = np.diag([(-1) ** ((basis >> 5) & 1) for basis in range(64)])
        column_z = np.diag([(-1) ** ((basis >> 3) & 1) for basis in range(64)])
        logical_z = np.diag([(-1) ** ((basis >> 4) & 1) for basis in range(64)])
        perturbation = (row_z @ column_z @ bit_flip_matrix(6, 4)
                        @ bit_flip_matrix(6, 1) @ bit_flip_matrix(6, 0))
        angle, phase = 0.05, 0.37
        source = (np.exp(1j * phase)
                  * (np.cos(angle) * identity + 1j * np.sin(angle) * perturbation)
                  @ ideal_source)
        np.testing.assert_allclose(source.conj().T @ source, identity, atol=1e-12)
        error = np.linalg.norm((np.exp(-1j * phase) * source - ideal_source)
                               @ initialization, 2)
        inverse_error = np.linalg.norm((np.exp(1j * phase) * source.conj().T - ideal_source)
                                       @ initialization, 2)
        self.assertLessEqual(inverse_error, error + 1e-12)
        self.assertGreater(np.linalg.norm((identity - clean_projector) @ source
                                          @ initialization, 2), 0.01)
        self.assertGreater(np.linalg.norm(source @ logical_z - logical_z @ source, 2), 0.01)

        actual = four_call_scalar_difference(source, row_flip, column_flip)
        returned_error = np.linalg.norm((actual - target) @ initialization, 2)
        self.assertGreater(returned_error, error)
        self.assertLessEqual(returned_error, 4 * error + 1e-12)
        # Four forward calls retain the source's fixed phase and are not the
        # actual-adjoint wrapper. This guards against that mistaken replacement.
        wrong_forward_word = (column_flip @ source @ row_flip @ source
                              @ column_flip @ source @ row_flip @ source)
        wrong_error = np.linalg.norm((wrong_forward_word - target) @ initialization, 2)
        self.assertGreater(wrong_error, max(4 * error, 0.5))

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
