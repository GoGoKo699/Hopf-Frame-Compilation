"""Convention-sensitive fixtures for the fixed-boundary decision screen.

These finite operator checks do not prove the asymptotic information bound.
The Choi vectors are mathematical witnesses, not source ancilla resources.
"""

import unittest

import numpy as np


def _choi_state(program):
    dimension = len(program)
    maximally_entangled = np.eye(dimension).reshape(-1) / np.sqrt(dimension)
    return np.kron(np.eye(dimension), program) @ maximally_entangled


class FixedBoundaryDecisionTests(unittest.TestCase):
    def test_complex_boundary_requires_the_transposed_choi_factor(self):
        pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
        pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        signal = np.diag([1., -1.])
        hadamard = (pauli_x + signal) / np.sqrt(2)
        phase = np.diag([1., 1j])
        t_gate = np.diag([1., np.exp(1j * np.pi / 4)])
        boundary = phase @ hadamard
        program = hadamard @ t_gate @ hadamard
        source = boundary.conj().T @ program @ boundary
        initialization = np.eye(2)[:, :1]
        compression = (initialization.conj().T @ source.conj().T
                       @ signal @ source @ initialization)
        transformed_signal = boundary @ signal @ boundary.conj().T
        np.testing.assert_allclose(transformed_signal, pauli_y, atol=1e-12)

        bias = np.trace(compression)
        correlation = np.trace(program.conj().T @ transformed_signal
                               @ program @ transformed_signal) / 2
        choi = _choi_state(program)
        measurement = np.kron(transformed_signal.T, transformed_signal)
        decoded_bias = np.vdot(choi, measurement @ choi)
        np.testing.assert_allclose([bias, correlation, decoded_bias],
                                   np.full(3, 1 / np.sqrt(2)), atol=1e-12)
        # Omitting the transpose, or using the adjoint instead, reverses
        # the answer: the transformed signal is Hermitian but imaginary.
        wrong_measurement = np.kron(transformed_signal, transformed_signal)
        self.assertAlmostEqual(np.vdot(choi, wrong_measurement @ choi).real,
                               -1 / np.sqrt(2), places=12)

    def test_code_motion_nonpauli_boundary_and_outer_clifford_absorption(self):
        # Physical order is clean c, arbitrary dirty d. The two addresses
        # give M_0=I and M_1=-I, even though F_0 moves the clean-zero code.
        identity = np.eye(2, dtype=complex)
        pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
        pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        pauli_z = np.diag([1., -1.])
        phase = np.diag([1., 1j])
        cnot = np.array([[1, 0, 0, 0], [0, 1, 0, 0],
                         [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)
        signal = np.kron(pauli_z, identity)
        initialization = np.eye(4)[:, :2]
        code_projector = initialization @ initialization.conj().T
        fixed_boundary = (np.cos(np.pi / 8) * np.eye(4)
                          + 1j * np.sin(np.pi / 8) * np.kron(pauli_x, pauli_x))
        outer_clifford = np.kron(identity, phase) @ cnot
        boundary = outer_clifford @ fixed_boundary
        generator = np.kron(pauli_z, pauli_z)
        mapped_generator = outer_clifford.conj().T @ generator @ outer_clifford
        np.testing.assert_allclose(mapped_generator, np.kron(identity, pauli_z),
                                   atol=1e-12)
        transformed_fixed_signal = fixed_boundary @ signal @ fixed_boundary.conj().T
        # Two orthogonal Pauli components certify that this signal is not
        # Pauli; the code-moving boundary is outside the isospectral screen.
        np.testing.assert_allclose(transformed_fixed_signal,
                                   (signal + np.kron(pauli_y, pauli_x)) / np.sqrt(2),
                                   atol=1e-12)
        self.assertAlmostEqual(np.linalg.norm((np.eye(4) - code_projector)
                                              @ fixed_boundary @ initialization, 2),
                               np.sin(np.pi / 8), places=12)

        absorbed_program = transformed_fixed_signal
        program = outer_clifford @ absorbed_program @ outer_clifford.conj().T
        for address in (0, 1):
            with self.subTest(address=address):
                gauge = np.linalg.matrix_power(generator, address)
                mapped_gauge = np.linalg.matrix_power(mapped_generator, address)
                source = (boundary.conj().T @ gauge @ program
                          @ gauge.conj().T @ boundary)
                absorbed_source = (fixed_boundary.conj().T @ mapped_gauge
                                   @ absorbed_program @ mapped_gauge.conj().T
                                   @ fixed_boundary)
                np.testing.assert_allclose(source, absorbed_source, atol=1e-12)
                compression = (initialization.conj().T @ source.conj().T
                               @ signal @ source @ initialization)
                np.testing.assert_allclose(compression, (-1) ** address * identity,
                                           atol=1e-12)
                transformed_signal = (gauge.conj().T @ boundary @ signal
                                      @ boundary.conj().T @ gauge)
                compressed_bias = np.trace(compression) / 2
                correlation = np.trace(program.conj().T @ transformed_signal
                                       @ program @ transformed_signal) / 4
                choi = _choi_state(program)
                decoded_bias = np.vdot(choi, np.kron(transformed_signal.T,
                                                   transformed_signal) @ choi)
                np.testing.assert_allclose([compressed_bias, correlation, decoded_bias],
                                           np.full(3, (-1) ** address), atol=1e-12)


if __name__ == "__main__":
    unittest.main()
