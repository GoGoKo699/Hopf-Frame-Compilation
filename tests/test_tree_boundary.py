"""Finite complete boundary-feedback checks; no native-emitter claim.

The original chronological edge word, inherited-charge recurrence, literal
rotation word, and all-port routing provide independent reference paths.
"""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.tree_boundary import (
    boundary_core,
    boundary_feedback,
    canonical_reversal_permutation,
    feedback_entrance,
    full_tree_coin_bank,
    prefix_encoder,
    boundary_incidence,
    coefficient_matrix,
    endpoint_block_ledger,
    propagation_block_encoding,
    propagation_dilation,
    propagation_matrices,
    shift_dilation_permutation,
)


ATOL = 2e-11


def _edges(height):
    return [(prefix << (height - depth),
             (2 * prefix + 1) << (height - depth - 1))
            for depth in range(height) for prefix in range(1 << depth)]


def _tables(height):
    count = (1 << height) - 1
    unequal = np.resize([0., .25, -1.6, .9, 0., -.7, 1.2], count)
    return (np.zeros(count), np.ones(count), unequal,
            np.resize([3., -2., 0., .125, -4.5], count))


def _chronological_phase_word(height, table):
    # Left multiplication updates the original physical edge rows. These
    # edges are not the differentiated boundary endpoints used by F.
    word = np.eye(1 << height, dtype=complex)
    for (left, right), value in zip(_edges(height), table):
        phase = (value + 1j) / (value - 1j)
        diagonal, off = (1 + phase) / 2, (1 - phase) / 2
        first, second = word[left].copy(), word[right].copy()
        word[left] = diagonal * first + off * second
        word[right] = off * first + diagonal * second
    return word


def _charge_recurrence(table):
    count = len(table)
    h = table / (1 + 1j * table)
    charges = np.zeros((count, count), complex)
    output = np.zeros_like(charges)
    for parent in range(count):
        output[parent] = -1j * h[parent] * charges[parent]
        output[parent, parent] += h[parent]
        left, right = 2 * parent + 1, 2 * parent + 2
        if right < count:
            charges[left] = -output[parent]
            charges[right] = charges[parent] + output[parent]
    return output


def _native_shift_destinations(height):
    size = 1 << height
    destinations = []
    for source in range(2 * size):
        bits = format(source, f"0{height + 1}b")
        word = int(bits[1:] + bits[:1], 2)
        # CNOT(port, lsb), zero-node-controlled X(port), CNOT(port, lsb).
        if word & size:
            word ^= 1
        if word % size == 0:
            word ^= size
        if word & size:
            word ^= 1
        destinations.append(word)
    return tuple(destinations)


def _permutation(destinations):
    matrix = np.zeros((len(destinations), len(destinations)))
    matrix[destinations, np.arange(len(destinations))] = 1
    return matrix


class TreeBoundaryTests(unittest.TestCase):
    def test_boundary_differentiation_and_chronological_triangular_k(self):
        for height in range(1, 6):
            with self.subTest(height=height):
                size = 1 << height
                cuts, boundary, triangular, decrement = boundary_incidence(height)
                expected_cuts = np.zeros((size, size - 1))
                expected_boundary = np.zeros_like(expected_cuts)
                for column, (left, right) in enumerate(_edges(height)):
                    length = right - left
                    expected_cuts[:, column] = -length / size
                    expected_cuts[right:right + length, column] += 1
                    expected_boundary[right - 1, column] = -1 / np.sqrt(2)
                    expected_boundary[right + length - 1, column] = 1 / np.sqrt(2)
                np.testing.assert_array_equal(cuts, expected_cuts)
                np.testing.assert_array_equal(boundary, expected_boundary)
                np.testing.assert_array_equal(np.argmax(decrement, axis=0),
                                              (np.arange(size) - 1) % size)
                baseline = (np.eye(size) - decrement) / 2
                np.testing.assert_allclose(np.sqrt(2) * baseline @ cuts,
                                           boundary, atol=ATOL, rtol=0)
                np.testing.assert_allclose(triangular, cuts.T @ baseline @ cuts,
                                           atol=ATOL, rtol=0)
                np.testing.assert_allclose(
                    triangular, np.eye(size - 1) / 2
                    + np.tril(boundary.T @ boundary, -1), atol=ATOL, rtol=0)
                np.testing.assert_allclose(cuts.sum(axis=0), 0, atol=ATOL, rtol=0)

    def test_complete_core_equals_original_chronological_edge_word(self):
        for height in range(1, 6):
            size = 1 << height
            uniform = np.ones(size) / np.sqrt(size)
            for index, table in enumerate(_tables(height)):
                with self.subTest(height=height, table=index):
                    actual = boundary_core(table)
                    expected = _chronological_phase_word(height, table)
                    np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(actual.conj().T @ actual, np.eye(size),
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(actual @ uniform, uniform,
                                               atol=ATOL, rtol=0)
            np.testing.assert_allclose(boundary_core(np.zeros(size - 1)),
                                       boundary_incidence(height)[3], atol=ATOL, rtol=0)

    def test_variable_coefficient_recurrence_including_zero_entries(self):
        for height in range(1, 6):
            count = (1 << height) - 1
            for index, table in enumerate(_tables(height)):
                with self.subTest(height=height, table=index):
                    coefficient = coefficient_matrix(table)
                    propagation, complement = propagation_matrices(table)
                    h = table / (1 + 1j * table)
                    resolvent = np.linalg.inv(np.eye(count) - propagation)
                    represented = (np.diag(h.real)
                                   - 1j * np.diag(h) @ resolvent @ np.diag(h.conj())
                                   - 1j * np.diag(h) @ resolvent @ complement @ np.diag(h))
                    np.testing.assert_allclose(coefficient, _charge_recurrence(table),
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(coefficient, represented, atol=ATOL, rtol=0)
                    self.assertLessEqual(np.linalg.norm(coefficient, 2),
                                         1 + np.sqrt(2) * height + ATOL)

    def test_uniform_root_chain_compression_and_norm_witness(self):
        for height in range(1, 6):
            count = (1 << height) - 1
            for tau in (-1., .5, 1.):
                with self.subTest(height=height, tau=tau):
                    table = np.full(count, tau)
                    coefficient = coefficient_matrix(table)
                    propagation, complement = propagation_matrices(table)
                    h = tau / (1 + 1j * tau)
                    resolvent = np.linalg.inv(np.eye(count) - propagation)
                    state = np.eye(count, dtype=complex)[:, 0]
                    columns = []
                    for _ in range(height):
                        columns.append(state)
                        state = propagation @ state
                    chain = np.column_stack(columns)
                    expected = (h.real * np.eye(height)
                                - 1j * abs(h)**2 * np.tril(np.ones((height, height))))
                    np.testing.assert_allclose(chain.conj().T @ chain, np.eye(height),
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(chain.conj().T @ coefficient @ chain,
                                               expected, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(chain.conj().T @ resolvent @ complement,
                                               0, atol=ATOL, rtol=0)
                    lower_squared = h.real**2 + abs(h)**4 * (height + 1) * (2 * height + 1) / 6
                    witness = chain @ (np.ones(height) / np.sqrt(height))
                    self.assertGreaterEqual(np.linalg.norm(coefficient @ witness)**2 + ATOL,
                                            lower_squared)
                    self.assertGreaterEqual(np.linalg.norm(coefficient, 2)**2 + ATOL,
                                            lower_squared)
                    self.assertAlmostEqual(np.linalg.norm(coefficient[:, 0])**2,
                                           abs(h)**2 + 2 * (height - 1) * abs(h)**4,
                                           delta=ATOL)

    def test_propagation_orthogonality_nilpotence_and_all_table_spectrum(self):
        for height in range(1, 6):
            count = (1 << height) - 1
            internal = (1 << (height - 1)) - 1
            active = np.diag([1.] * internal + [0.] * (count - internal))
            lengths = [height]
            for length in range(1, height):
                lengths.extend([length] * (1 << (height - length - 1)))
            expected = sorted(2 * np.sin((2 * k - 1) * np.pi / (4 * length + 2))
                              for length in lengths for k in range(1, length + 1))
            self.assertEqual(len(expected), count)
            for index, table in enumerate(_tables(height)):
                with self.subTest(height=height, table=index):
                    propagation, complement = propagation_matrices(table)
                    np.testing.assert_allclose(propagation.conj().T @ propagation,
                                               active, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(complement.conj().T @ complement,
                                               active, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(propagation.conj().T @ complement,
                                               0, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(np.linalg.matrix_power(propagation, height),
                                               0, atol=ATOL, rtol=0)
                    singular = np.linalg.svd(np.eye(count) - propagation, compute_uv=False)
                    np.testing.assert_allclose(sorted(singular), expected, atol=ATOL, rtol=0)
                    self.assertAlmostEqual(min(singular), 2 * np.sin(np.pi / (4 * height + 2)),
                                           delta=ATOL)

    def test_shift_dilation_routes_every_port_and_dummy_by_native_word(self):
        for height in range(1, 7):
            with self.subTest(height=height):
                size = 1 << height
                destinations = shift_dilation_permutation(height)
                self.assertEqual(destinations, _native_shift_destinations(height))
                self.assertEqual(sorted(destinations), list(range(2 * size)))
                self.assertEqual(destinations[0], size + 1)
                self.assertEqual(destinations[size + size // 2], 0)
                expected_shift = np.zeros((size, size))
                for parent in range(1, size // 2):
                    expected_shift[2 * parent, parent] = 1
                np.testing.assert_array_equal(_permutation(destinations)[:size, :size],
                                              expected_shift)

    def test_full_propagation_dilation_and_literal_rx_rz_factors(self):
        xz = np.array([[0., -1.], [1., 0.]])
        for height in range(1, 5):
            size = 1 << height
            route = _permutation(_native_shift_destinations(height))
            for index, table in enumerate(_tables(height)):
                with self.subTest(height=height, table=index):
                    bank = np.eye(size, dtype=complex)
                    for parent in range(1, size // 2):
                        alpha = np.arctan(table[parent - 1])
                        cosine, sine = np.cos(alpha), np.sin(alpha)
                        rx = np.array([[cosine, -1j * sine], [-1j * sine, cosine]])
                        rz = np.diag([np.exp(-1j * alpha), np.exp(1j * alpha)])
                        bank[2 * parent:2 * parent + 2, 2 * parent:2 * parent + 2] = xz @ rx @ rz
                    actual = propagation_dilation(table)
                    np.testing.assert_allclose(actual, np.kron(np.eye(2), bank) @ route,
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(actual.conj().T @ actual, np.eye(2 * size),
                                               atol=ATOL, rtol=0)
                    padded = np.zeros((size, size), complex)
                    padded[1:, 1:] = propagation_matrices(table)[0]
                    np.testing.assert_allclose(actual[:size, :size], padded, atol=ATOL, rtol=0)
                    # Undoing the route exposes identity on both dummy pairs.
                    unrouted = actual @ route.T
                    for port in range(2):
                        start = port * size
                        np.testing.assert_array_equal(unrouted[start:start + 2, start:start + 2],
                                                      np.eye(2))

    def test_two_flag_lcu_word_on_every_port_and_accepted_block(self):
        hadamard = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
        for height in range(1, 5):
            size = 1 << height
            hc = np.kron(hadamard, np.eye(2 * size))
            zc = np.kron(np.diag([1., -1.]), np.eye(2 * size))
            for index, table in enumerate(_tables(height)):
                with self.subTest(height=height, table=index):
                    controlled = np.eye(4 * size, dtype=complex)
                    controlled[2 * size:, 2 * size:] = propagation_dilation(table)
                    actual = propagation_block_encoding(table)
                    np.testing.assert_allclose(actual, hc @ zc @ controlled @ hc,
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(actual.conj().T @ actual, np.eye(4 * size),
                                               atol=ATOL, rtol=0)
                    padded = np.zeros((size, size), complex)
                    padded[1:, 1:] = propagation_matrices(table)[0]
                    np.testing.assert_allclose(actual[:size, :size], (np.eye(size) - padded) / 2,
                                               atol=ATOL, rtol=0)
                    # The compiled component has rejected-port amplitude;
                    # accepting the block does not return the two flags.
                    self.assertGreater(np.linalg.norm(actual[size:, 0]), .5)

    def test_complete_feedback_and_initialized_encoder_cancellation(self):
        for height in range(1, 6):
            size = 1 << height
            dim = 2 * size
            canonical = _permutation(canonical_reversal_permutation(height))
            for index, table in enumerate(_tables(height)):
                with self.subTest(height=height, table=index):
                    bank = full_tree_coin_bank(table)
                    encoder = prefix_encoder(table)
                    # Construct weighted chains independently from the recurrence.
                    propagation = np.zeros((dim, dim), dtype=complex)
                    for node in range(1, size):
                        propagation[:, node] = bank[:, 2 * node]
                    independent = np.eye(dim, dtype=complex)
                    for node in range(1, dim):
                        age = (node & -node).bit_length() - 1
                        root = node >> age
                        vector = (np.eye(dim)[:, 1].astype(complex) if root == 1
                                  else bank[:, root].copy())
                        for _ in range(age):
                            vector = propagation @ vector
                        independent[:, node] = vector
                    np.testing.assert_allclose(encoder, independent, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(encoder.conj().T @ encoder, np.eye(dim),
                                               atol=ATOL, rtol=0)
                    reversal = encoder @ canonical @ encoder.conj().T
                    np.testing.assert_allclose(reversal @ reversal, np.eye(dim),
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(reversal, reversal.conj().T,
                                               atol=ATOL, rtol=0)
                    incoming, diagonal = feedback_entrance(table)
                    entrance = bank @ incoming @ np.kron(np.eye(2), diagonal)
                    np.testing.assert_allclose(
                        (encoder.conj().T @ entrance)[:, :size],
                        (incoming @ np.kron(np.eye(2), diagonal))[:, :size],
                        atol=ATOL, rtol=0)
                    actual = boundary_feedback(table)
                    expected = _chronological_phase_word(height, table)
                    np.testing.assert_allclose(actual[:, :size],
                                               np.vstack((expected, np.zeros_like(expected))),
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(actual.conj().T @ actual, np.eye(dim),
                                               atol=ATOL, rtol=0)

    def test_canonical_reversal_by_endpoint_bit_reversals(self):
        for height in range(1, 9):
            expected = [0]
            for source in range(1, 2 << height):
                bits = list(format(source, f"0{height + 1}b")[::-1])
                first = bits.index('1')
                last = len(bits) - 1 - bits[::-1].index('1')
                if last > first:
                    bits[first + 1:last] = bits[first + 1:last][::-1]
                expected.append(int(''.join(bits), 2))
            self.assertEqual(canonical_reversal_permutation(height), tuple(expected))
            self.assertEqual(sorted(expected), list(range(2 << height)))

    def test_repeated_block_width_and_full_operator_error(self):
        for height in (3, 4, 8, 20):
            for overhead in range(height - 2):
                for calls in sorted({1 << overhead, max(1, (1 << overhead) - 1)}):
                    actual_overhead = (calls - 1).bit_length()
                    for controls in (0, height - 3 - actual_overhead):
                        with self.subTest(height=height, calls=calls, controls=controls):
                            size = 1 << height
                            ledger = endpoint_block_ledger(height, calls, controls)
                            self.assertEqual(ledger.source_width, size + actual_overhead + 9)
                            self.assertEqual(ledger.core_qubits, ledger.source_width + 1)
                            self.assertGreaterEqual(ledger.core_qubits, 2 * (height - 1))
                            self.assertEqual(ledger.selector_qubits, 0)
                            self.assertEqual(ledger.predicate_helper_qubits, 0)
                            self.assertEqual(ledger.borrowed_signal_qubits, 0)
                            self.assertEqual(ledger.clean_signal_qubits, 2)
                            self.assertEqual(ledger.dirty_qubits, ledger.core_qubits + controls)
                            self.assertLessEqual(ledger.dirty_qubits, size + height + 7)
                            self.assertEqual(ledger.total_width, height + 2 + ledger.dirty_qubits)
                            self.assertEqual(ledger.sectors_per_bank, 1)
                            self.assertEqual(ledger.variable_banks, 2)
                            self.assertEqual(ledger.predicate_literals, 1 + controls)
                            self.assertEqual(ledger.free_table_rows, size // 2)
                            self.assertEqual(ledger.error_ratio,
                                             Fraction(86 * calls, 1 << (actual_overhead + 9)))
                            self.assertLessEqual(ledger.error_ratio, Fraction(43, 256))

    def test_invalid_tables_and_heights_are_rejected(self):
        for function in (coefficient_matrix, propagation_matrices, boundary_core,
                         propagation_dilation, propagation_block_encoding,
                         full_tree_coin_bank, prefix_encoder, feedback_entrance, boundary_feedback):
            for table in ([], [1., 2.], [np.nan], [np.inf], [-np.inf], [1j]):
                with self.subTest(function=function.__name__, table=table):
                    with self.assertRaises(ValueError):
                        function(table)
        for function in (boundary_incidence, shift_dilation_permutation,
                         canonical_reversal_permutation):
            for height in (-1, 0):
                with self.subTest(function=function.__name__, height=height):
                    with self.assertRaises(ValueError):
                        function(height)
        for height in (-1, 0, 1, 2):
            with self.subTest(ledger_height=height):
                with self.assertRaises(ValueError):
                    endpoint_block_ledger(height)
        for height, calls, controls in ((3, 2, 0), (4, 2, 1), (5, 0, 0), (5, 1, -1)):
            with self.assertRaises(ValueError):
                endpoint_block_ledger(height, calls, controls)


if __name__ == '__main__':
    unittest.main()
