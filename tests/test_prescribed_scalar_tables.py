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


def cnot_matrix(width, control, target):
    """Literal CNOT permutation, with bit zero the least significant wire."""
    size = 1 << width
    permutation = [basis ^ (((basis >> control) & 1) << target)
                   for basis in range(size)]
    return np.eye(size, dtype=int)[:, permutation]


def one_clean_decision_dilation(compression):
    """Matrix fixture with J-dagger V-dagger Z_c V J equal to compression.

    The clean signal is the first wire. This canonical unitary completion
    is used to test the recovery identity, not as an unpaid native source.
    """
    values, vectors = np.linalg.eigh(compression)
    accepted = ((vectors * np.sqrt(np.clip((1 + values) / 2, 0, 1)))
                @ vectors.conj().T)
    rejected = ((vectors * np.sqrt(np.clip((1 - values) / 2, 0, 1)))
                @ vectors.conj().T)
    return np.block([[accepted, -rejected], [rejected, accepted]])


def four_call_scalar_difference(source, row_block_flip, column_block_flip):
    """Chronological V, X_a, V-dagger, X_b, V, X_a, V-dagger, X_b."""
    reverse = source.conj().T
    return (column_block_flip @ reverse @ row_block_flip @ source
            @ column_block_flip @ reverse @ row_block_flip @ source)


def involution_rotation(involution, angle):
    """Exact matrix fixture for exp(i angle S), with S Hermitian and S^2=I."""
    return (np.cos(angle) * np.eye(len(involution))
            + 1j * np.sin(angle) * involution)


def cubic_pair_reflection_word(signal, reflection):
    """A matrix QSP fixture, not a native emitter or a recursive compiler."""
    word = involution_rotation(reflection, -np.pi / 3)
    for angle in (np.pi / 3, np.pi / 3, -np.pi / 3):
        word = word @ signal @ involution_rotation(reflection, angle)
    return word


class PrescribedScalarTablesTests(unittest.TestCase):
    def test_one_t_clean_decision_has_a_different_rejected_sector_sign(self):
        # Each source is an actual word of CNOT, H, X and exactly one T.
        # Physical order: address, clean c. Extra address CNOTs and an X
        # extend the two-wire construction to nonconstant affine parities.
        hadamard = normalized_hadamard(1)
        gamma = 1 / np.sqrt(2)
        for address_bits, parity_mask, offset in ((1, 1, 0), (2, 3, 0), (3, 5, 1)):
            with self.subTest(address_bits=address_bits, mask=parity_mask, offset=offset):
                width = address_bits + 1
                size = 1 << width
                identity = np.eye(size)
                initialization = identity[:, ::2]
                clean_z = np.diag([(-1) ** (basis & 1) for basis in range(size)])
                source_c = cnot_matrix(width, 1, 0) @ cnot_matrix(width, 0, 1)
                local_h = np.kron(np.eye(1 << (address_bits - 1)),
                                  np.kron(hadamard, np.eye(2)))
                source_k = cnot_matrix(width, 1, 0) @ local_h
                local_t = np.diag([np.exp(1j * np.pi / 4) if basis & 2 else 1
                                   for basis in range(size)])
                address_change = identity
                for address_bit in range(1, address_bits):
                    if parity_mask & (1 << address_bit):
                        address_change = (cnot_matrix(width, address_bit + 1, 1)
                                          @ address_change)
                if offset:
                    address_change = bit_flip_matrix(width, 1) @ address_change
                source = source_c @ source_k @ local_t @ source_k.conj().T @ address_change
                signs = np.array([(-1) ** ((address & parity_mask).bit_count() + offset)
                                  for address in range(1 << address_bits)])
                target = np.diag(signs)
                full_target = np.kron(target, np.eye(2))

                np.testing.assert_allclose(source.conj().T @ source, identity, atol=1e-12)
                reflection = source.conj().T @ clean_z @ source
                address_z = np.diag([(-1) ** ((basis >> 1) & 1)
                                     for basis in range(size)])
                address_y_signal_x = np.kron(
                    np.eye(1 << (address_bits - 1)),
                    np.kron(np.array([[0, -1j], [1j, 0]]), np.array([[0., 1.], [1., 0.]])))
                expected_reflection = (address_change.conj().T
                                       @ (gamma * (address_z + address_y_signal_x))
                                       @ address_change)
                np.testing.assert_allclose(reflection, expected_reflection, atol=1e-12)
                compression = initialization.conj().T @ reflection @ initialization
                accepted = initialization.conj().T @ source @ initialization
                np.testing.assert_allclose(compression, gamma * target, atol=1e-12)
                np.testing.assert_allclose(compression,
                                           2 * accepted.conj().T @ accepted - np.eye(len(target)),
                                           atol=1e-12)
                pair = reflection @ clean_z
                hermitian_pair = (pair + pair.conj().T) / 2
                np.testing.assert_allclose(hermitian_pair,
                                           gamma * full_target @ clean_z, atol=1e-12)
                self.assertAlmostEqual(np.linalg.norm(
                    hermitian_pair / gamma - full_target, 2), 2, places=12)
                clean_projector = initialization @ initialization.conj().T
                self.assertGreater(np.linalg.norm(
                    (identity - clean_projector) @ source @ initialization, 2), 0.8)
                self.assertGreater(np.linalg.norm(
                    source @ full_target - full_target @ source, 2), 0.2)

    def test_clean_decision_amplifies_nonflat_work_with_address_motion(self):
        # Physical order: clean c, address x0,x1, arbitrary dirty d.
        # M mixes x1 and depends on d, but has the exact sign Z_x0.
        # A right change acts only on the rejected input sector, while a
        # left change moves addresses and dirty outputs on the clean code.
        # Both are allowed; the raw source neither returns its clean bit
        # nor satisfies the prescribed sign on the whole physical space.
        pauli_x = np.array([[0., 1.], [1., 0.]])
        pauli_z = np.diag([1., -1.])
        hadamard = normalized_hadamard(1)
        code_identity = np.eye(8)
        identity = np.eye(16)
        initialization = identity[:, :8]
        clean_z = np.kron(pauli_z, code_identity)
        target = np.kron(pauli_z, np.eye(4))
        first_axis = np.kron(np.eye(2), np.kron(pauli_x, pauli_z))
        second_axis = np.kron(np.eye(2), np.kron(pauli_z, pauli_x))
        magnitude = 0.75 * code_identity + 0.15 * first_axis + 0.05 * second_axis
        expected_compression = target @ magnitude
        rejected_change = np.kron(hadamard, np.eye(4))
        off_code_change = np.block([[code_identity, np.zeros((8, 8))],
                                    [np.zeros((8, 8)), rejected_change]])
        output_motion = np.kron(np.eye(2),
                                np.kron(hadamard, np.kron(np.eye(2), hadamard)))
        source = (output_motion @ one_clean_decision_dilation(expected_compression)
                  @ off_code_change)
        reflection = source.conj().T @ clean_z @ source
        compression = initialization.conj().T @ reflection @ initialization
        pair = reflection @ clean_z
        hermitian_pair = (pair + pair.conj().T) / 2
        full_target = np.kron(np.eye(2), target)
        address_z1 = np.kron(np.eye(2), np.kron(pauli_z, np.eye(2)))
        np.testing.assert_allclose(source.conj().T @ source, identity, atol=1e-12)
        np.testing.assert_allclose(compression, expected_compression, atol=1e-12)
        np.testing.assert_allclose(hermitian_pair @ initialization,
                                   initialization @ compression, atol=1e-12)
        np.testing.assert_allclose(compression @ target, target @ compression, atol=1e-12)
        spectrum = np.linalg.eigvalsh(target @ compression)
        self.assertAlmostEqual(spectrum[0], 0.55, places=12)
        self.assertAlmostEqual(spectrum[-1], 0.95, places=12)
        self.assertGreater(np.linalg.norm(
            compression @ address_z1 - address_z1 @ compression, 2), 0.2)
        self.assertGreater(np.linalg.norm(
            source @ full_target - full_target @ source, 2), 1)
        self.assertGreater(np.linalg.norm(
            hermitian_pair @ full_target - full_target @ hermitian_pair, 2), 1)

        amplified = cubic_pair_reflection_word(pair, clean_z)
        accepted = initialization.conj().T @ amplified @ initialization
        polynomial = (3 * compression - np.linalg.matrix_power(compression, 3)) / 2
        error = amplified @ initialization - initialization @ target
        np.testing.assert_allclose((accepted + accepted.conj().T) / 2,
                                   polynomial, atol=1e-12)
        np.testing.assert_allclose(accepted @ compression, compression @ accepted, atol=1e-12)
        np.testing.assert_allclose(error.conj().T @ error,
                                   2 * (code_identity - target @ polynomial), atol=1e-12)
        self.assertLess(np.linalg.norm(error, 2),
                        np.linalg.norm(pair @ initialization - initialization @ target, 2))
        # This is a full initialized-isometry comparison; it includes the
        # surviving rejected-signal amplitude instead of discarding it.
        projector = initialization @ initialization.conj().T
        self.assertGreater(np.linalg.norm((identity - projector)
                                          @ amplified @ initialization, 2), 0.1)
        self.assertGreater(np.linalg.norm(amplified - full_target, 2), 1)

    def test_basis_decision_bias_does_not_control_coherent_dirty_inputs(self):
        # All computational dirty inputs report the right value with bias
        # 1/4. The dirty |-> state instead has signed bias -1/2. A spectral
        # sign filter strengthens this wrong answer rather than correcting
        # it, so basis probabilities cannot establish the source premise.
        pauli_x = np.array([[0., 1.], [1., 0.]])
        pauli_z = np.diag([1., -1.])
        target = np.kron(pauli_z, np.eye(2))
        compression = np.kron(pauli_z, 0.25 * np.eye(2) + 0.75 * pauli_x)
        source = one_clean_decision_dilation(compression)
        identity = np.eye(8)
        initialization = identity[:, :4]
        clean_z = np.kron(pauli_z, np.eye(4))
        reflection = source.conj().T @ clean_z @ source
        actual_compression = initialization.conj().T @ reflection @ initialization
        np.testing.assert_allclose(actual_compression, compression, atol=1e-12)
        np.testing.assert_allclose(np.diag(target) * np.diag(compression),
                                   np.full(4, 0.25), atol=1e-12)
        self.assertAlmostEqual(np.linalg.eigvalsh(target @ compression)[0], -0.5, places=12)

        pair = reflection @ clean_z
        amplified = cubic_pair_reflection_word(pair, clean_z)
        wrong_input = np.array([1., -1., 0., 0.]) / np.sqrt(2)
        ideal = initialization @ target @ wrong_input
        initial_error = np.linalg.norm(pair @ initialization @ wrong_input - ideal)
        final_error = np.linalg.norm(amplified @ initialization @ wrong_input - ideal)
        self.assertAlmostEqual(initial_error ** 2, 3, places=12)
        self.assertAlmostEqual(final_error ** 2, 27 / 8, places=12)
        self.assertGreater(final_error, initial_error)

    def test_signed_gap_allows_nonflat_spectrum_and_completed_address_motion(self):
        # Physical order: address x0,x1, arbitrary signal t, arbitrary dirty d.
        # V moves x0 across the target sign sectors, but that left motion
        # commutes with X_t and cancels from V-dagger X_t V. Its remaining
        # perturbation mixes x1 and depends on d, so H is nonflat even within
        # a sign sector. This checks a genuinely weaker source premise than
        # completed computational-address preservation. Analytic rotations
        # are a matrix fixture; the theorem separately pays native synthesis.
        dimension = 16
        identity = np.eye(dimension)
        signal = bit_flip_matrix(4, 1)
        address_x1 = bit_flip_matrix(4, 2)
        dirty_z = np.diag([(-1) ** (basis & 1) for basis in range(dimension)])
        projector = (identity + address_x1 @ dirty_z) / 2
        pauli_y = np.array([[0, -1j], [1j, 0]])
        signal_y = np.kron(np.eye(4), np.kron(pauli_y, np.eye(2)))
        angle = np.pi / 6
        perturbation = (identity + (np.cos(angle) - 1) * projector
                        + 1j * np.sin(angle) * projector @ signal_y)
        controlled_sign = np.diag([
            (-1) ** (((basis >> 3) & 1) * ((basis >> 1) & 1))
            for basis in range(dimension)
        ])
        address_motion = np.kron(normalized_hadamard(1), np.eye(8))
        source = address_motion @ controlled_sign @ perturbation
        target = np.diag([(-1) ** ((basis >> 3) & 1) for basis in range(dimension)])
        address_z1 = np.diag([(-1) ** ((basis >> 2) & 1) for basis in range(dimension)])
        np.testing.assert_allclose(source.conj().T @ source, identity, atol=1e-12)
        self.assertGreater(np.linalg.norm(source @ target - target @ source, 2), 1)

        pair = source.conj().T @ signal @ source @ signal
        hermitian_pair = (pair + pair.conj().T) / 2
        expected = target @ (identity + (np.cos(2 * angle) - 1) * projector)
        np.testing.assert_allclose(hermitian_pair, expected, atol=1e-12)
        self.assertGreater(np.linalg.norm(
            hermitian_pair @ address_z1 - address_z1 @ hermitian_pair, 2), 0.4)
        np.testing.assert_allclose(hermitian_pair @ target,
                                   target @ hermitian_pair, atol=1e-12)
        signed_spectrum = np.linalg.eigvalsh(target @ hermitian_pair)
        self.assertAlmostEqual(signed_spectrum[0], 0.5, places=12)
        self.assertAlmostEqual(signed_spectrum[-1], 1, places=12)

        amplified = cubic_pair_reflection_word(pair, signal)
        polynomial = (3 * hermitian_pair - np.linalg.matrix_power(hermitian_pair, 3)) / 2
        np.testing.assert_allclose((amplified + amplified.conj().T) / 2,
                                   polynomial, atol=1e-12)
        np.testing.assert_allclose(amplified @ target, target @ amplified, atol=1e-12)
        np.testing.assert_allclose(
            (amplified - target).conj().T @ (amplified - target),
            2 * (identity - target @ polynomial), atol=1e-12)
        self.assertLess(np.linalg.norm(amplified - target, 2),
                        np.linalg.norm(pair - target, 2))

    def test_cubic_pair_reflection_on_a_full_dirty_source(self):
        # Physical order: two address bits x, auxiliary bit t, arbitrary
        # dirty bit d. Every source column is included. The source preserves
        # x, but its entangling perturbation moves both t and d. These are
        # complete operator norms, so an external reference is also allowed.
        identity = np.eye(2)
        pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
        pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        pauli_z = np.diag([1., -1.])
        table = (0, 1, 1, 1)
        source = np.zeros((16, 16), dtype=complex)
        ideal_source = np.zeros_like(source)
        phase = 0.37
        for address, value in enumerate(table):
            ideal_block = np.kron(pauli_z if value else identity, identity)
            perturbation = (
                involution_rotation(np.kron(pauli_z, identity), 0.03 + 0.008 * address)
                @ involution_rotation(np.kron(pauli_z, pauli_z), 0.02 - 0.003 * address)
                @ involution_rotation(np.kron(pauli_x, pauli_y), 0.045 + 0.005 * address)
            )
            block = slice(4 * address, 4 * (address + 1))
            source[block, block] = np.exp(1j * phase) * perturbation @ ideal_block
            ideal_source[block, block] = ideal_block

        full_identity = np.eye(16)
        reflection = np.kron(np.eye(4), np.kron(pauli_x, identity))
        target = np.kron(np.diag(1 - 2 * np.array(table)), np.eye(4))
        auxiliary_z = np.kron(np.eye(4), np.kron(pauli_z, identity))
        dirty_z = np.kron(np.eye(8), pauli_z)
        np.testing.assert_allclose(source.conj().T @ source, full_identity, atol=1e-12)
        for address_z in (np.kron(pauli_z, np.eye(8)),
                          np.kron(identity, np.kron(pauli_z, np.eye(4)))):
            np.testing.assert_allclose(source @ address_z, address_z @ source, atol=1e-12)
        self.assertGreater(np.linalg.norm(source @ auxiliary_z - auxiliary_z @ source, 2), 0.05)
        self.assertGreater(np.linalg.norm(source @ dirty_z - dirty_z @ source, 2), 0.05)

        source_error = np.linalg.norm(np.exp(-1j * phase) * source - ideal_source, 2)
        signal = source.conj().T @ reflection @ source @ reflection
        phase_free = np.exp(-1j * phase) * source
        np.testing.assert_allclose(
            signal, phase_free.conj().T @ reflection @ phase_free @ reflection, atol=1e-12)
        np.testing.assert_allclose(reflection @ signal @ reflection, signal.conj().T, atol=1e-12)
        initial_error = np.linalg.norm(signal - target, 2)
        self.assertGreater(initial_error, 0.05)
        self.assertLessEqual(initial_error, 2 * source_error + 1e-12)

        word = cubic_pair_reflection_word(signal, reflection)
        hermitian_signal = (signal + signal.conj().T) / 2
        polynomial = (3 * hermitian_signal - np.linalg.matrix_power(hermitian_signal, 3)) / 2
        np.testing.assert_allclose((word + word.conj().T) / 2, polynomial, atol=1e-12)
        np.testing.assert_allclose(word.conj().T @ word, full_identity, atol=1e-12)
        # The all-size construction factorizes one degree-k QSP word. It
        # cannot recursively reuse this cubic output with the same S: the
        # pair-reflection relation has not been retained by that output.
        self.assertGreater(np.linalg.norm(reflection @ word @ reflection - word.conj().T, 2),
                           0.001)
        # The signal is not just a two-dimensional scalar-angle test copied
        # over spectators: its Hermitian part depends on the dirty state.
        first_block = hermitian_signal[:4, :4]
        self.assertGreater(np.linalg.norm(first_block - np.trace(first_block) * np.eye(4) / 4, 2),
                           0.001)

        final_error = np.linalg.norm(word - target, 2)
        self.assertLess(final_error, initial_error / 4)
        self.assertLessEqual(final_error, np.sqrt(3) * initial_error ** 2 / 2 + 1e-12)
        self.assertAlmostEqual(final_error ** 2,
                               2 * np.linalg.norm(full_identity - target @ polynomial, 2), places=12)
        # Forward substitution would retain the source's unknown common phase.
        wrong_signal = source @ reflection @ source @ reflection
        self.assertGreater(np.linalg.norm(wrong_signal - target, 2), 0.5)

    def test_cubic_pair_reflection_keeps_one_dimensional_endpoints(self):
        # The simultaneous one-dimensional reflection blocks can have either
        # reflection sign and either W endpoint. No extra phase is permitted.
        for eigenvalue, reflection_sign in product((-1, 1), repeat=2):
            with self.subTest(eigenvalue=eigenvalue, reflection=reflection_sign):
                signal = np.array([[eigenvalue]], dtype=complex)
                reflection = np.array([[reflection_sign]], dtype=complex)
                word = cubic_pair_reflection_word(signal, reflection)
                np.testing.assert_allclose(word, signal, atol=1e-12)

    def test_cubic_pair_reflection_cannot_remove_an_address_basis_tilt(self):
        # A small coherent address rotation gives an exact involution with
        # the wrong spectral subspaces. Endpoint-preserving QSP keeps those
        # subspaces, even though V is close to the intended scalar source.
        # Physical order: x, t, arbitrary dirty d.
        angle = 0.07
        address_rotation = np.array([[np.cos(angle), -np.sin(angle)],
                                     [np.sin(angle), np.cos(angle)]])
        tilt = np.kron(address_rotation, np.eye(4))
        ideal_source = np.diag([(-1) ** (((basis >> 2) & 1) * ((basis >> 1) & 1))
                                for basis in range(8)])
        source = np.exp(0.37j) * tilt @ ideal_source @ tilt.conj().T
        reflection = bit_flip_matrix(3, 1)
        target = np.diag([(-1) ** ((basis >> 2) & 1) for basis in range(8)])
        self.assertGreater(np.linalg.norm(source @ target - target @ source, 2), 0.1)
        self.assertLess(np.linalg.norm(np.exp(-0.37j) * source - ideal_source, 2), 0.15)

        signal = source.conj().T @ reflection @ source @ reflection
        np.testing.assert_allclose(signal, tilt @ target @ tilt.conj().T, atol=1e-12)
        np.testing.assert_allclose(signal @ signal, np.eye(8), atol=1e-12)
        word = cubic_pair_reflection_word(signal, reflection)
        np.testing.assert_allclose(word, signal, atol=1e-12)
        initial_error = np.linalg.norm(signal - target, 2)
        final_error = np.linalg.norm(word - target, 2)
        self.assertAlmostEqual(initial_error, 2 * np.sin(angle), places=12)
        self.assertAlmostEqual(final_error, initial_error, places=12)
        self.assertGreater(final_error, np.sqrt(3) * initial_error ** 2 / 2)

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
