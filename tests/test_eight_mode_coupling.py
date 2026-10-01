"""Eight-mode limits of the childwise SO(4) magic-basis factorization.

Matrices have dimension at most eight. Parameterized Pauli rotations are
charged synthesis targets, not free native gates. One special angle has
an explicit four-T word including its literal global-phase correction.
The product-gate distance concerns the prefix/child tensor cut only.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_one_clean_compiler import _pauli_t_word
from tests.test_operator_source_compiler import (
    I2, X, Y, Z, _adjoint, _pauli, _rotation, _word_matrix,
)


ATOL = 2e-12
PAULIS = {'X': X, 'Y': Y, 'Z': Z}
ROOT_TERMS = (({2: 'Y'}, 1), ({2: 'Y', 1: 'X', 0: 'X'}, 1),
              ({2: 'Y', 1: 'Y', 0: 'Y'}, -1),
              ({2: 'Y', 1: 'Z', 0: 'Z'}, 1))


def _magic_word():
    # The child MSB is wire 1 here; wire 0 is the low integer bit.
    # Chronological CZ, S on both, H on the MSB, then CNOT MSB -> LSB.
    return [('H', 0), ('CX', 1, 0), ('H', 0), ('S', 1), ('S', 0),
            ('H', 1), ('CX', 1, 0)]


def _root_rotation(theta):
    result = np.eye(8, dtype=complex)
    result[np.ix_([0, 4], [0, 4])] = _rotation(theta)
    return result


def _pair_generator(low, high):
    result = np.zeros((8, 8), complex)
    result[high, low], result[low, high] = 1j, -1j
    return result


def _terms():
    return [(sign, _pauli(3, {wire: PAULIS[label] for wire, label in factors.items()}))
            for factors, sign in ROOT_TERMS]


class EightModeCouplingTests(unittest.TestCase):
    def test_actual_magic_clifford_and_four_commuting_rotations(self):
        magic = _word_matrix(2, _magic_word())
        bell = np.array([1, 0, 0, 1]) / np.sqrt(2)
        np.testing.assert_allclose(magic[:, 0], bell, atol=ATOL, rtol=0)
        projector = np.outer(bell, bell.conj())
        np.testing.assert_allclose(projector,
                                   (np.eye(4) + np.kron(X, X) - np.kron(Y, Y) + np.kron(Z, Z)) / 4,
                                   atol=ATOL, rtol=0)
        full_magic = np.kron(I2, magic)
        generator = full_magic @ _pair_generator(0, 4) @ full_magic.conj().T
        np.testing.assert_allclose(generator, np.kron(Y, projector), atol=ATOL, rtol=0)
        terms = _terms()
        np.testing.assert_allclose(generator, sum(sign * pauli for sign, pauli in terms) / 4,
                                   atol=ATOL, rtol=0)
        for _, left in terms:
            for _, right in terms:
                np.testing.assert_array_equal(left @ right, right @ left)
        for theta in (0, .15, np.pi / 4, np.pi / 2):
            emitted = np.eye(8, dtype=complex)
            for sign, pauli in terms:
                emitted = (np.cos(theta / 4) * np.eye(8)
                           - 1j * sign * np.sin(theta / 4) * pauli) @ emitted
            expected = full_magic @ _root_rotation(theta) @ full_magic.conj().T
            np.testing.assert_allclose(emitted, expected, atol=ATOL, rtol=0)

    def test_special_angle_native_word_keeps_the_literal_scalar_phase(self):
        factors = []
        for paulis, sign in ROOT_TERMS:
            factors += _pauli_t_word(paulis, sign)
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in factors), 4)
        phase_correction = [('H', 0), ('SDG', 0)] * 3
        np.testing.assert_allclose(_word_matrix(1, phase_correction),
                                   np.exp(-1j * np.pi / 4) * I2, atol=ATOL, rtol=0)
        magic = _magic_word()
        actual = _word_matrix(3, magic + factors + phase_correction + _adjoint(magic))
        np.testing.assert_allclose(actual, _root_rotation(np.pi / 2), atol=ATOL, rtol=0)
        without_correction = _word_matrix(3, magic + factors + _adjoint(magic))
        np.testing.assert_allclose(without_correction,
                                   np.exp(1j * np.pi / 4) * _root_rotation(np.pi / 2),
                                   atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(without_correction - actual, 2), .7)

    def test_product_distance_schmidt_witness_and_attaining_product(self):
        magic = _word_matrix(2, _magic_word())
        full_magic = np.kron(I2, magic)
        phi, chi = magic[:, 0], magic[:, 1]
        product_input = np.kron([1, 0], (phi + chi) / np.sqrt(2))
        for theta in (0, .03, .15, np.pi / 4, np.pi / 2):
            with self.subTest(theta=theta):
                target = full_magic @ _root_rotation(theta) @ full_magic.conj().T
                output = target @ product_input
                singular = np.linalg.svd(output.reshape(2, 4), compute_uv=False)
                np.testing.assert_allclose(singular, [np.cos(theta / 2), np.sin(theta / 2)],
                                           atol=ATOL, rtol=0)
                # Every product gate maps this input to a product vector;
                # the largest Schmidt coefficient bounds its overlap.
                lower_squared = 2 - 2 * singular[0]
                expected = 2 * np.sin(theta / 4)
                self.assertAlmostEqual(lower_squared, expected * expected, delta=ATOL)
                attaining = np.kron(_rotation(theta / 2), np.eye(4))
                self.assertAlmostEqual(np.linalg.norm(target - attaining, 2), expected, delta=ATOL)

    def test_root_and_left_child_commutator_retains_cross_branch_action(self):
        magic = np.kron(I2, _word_matrix(2, _magic_word()))
        root = magic @ _pair_generator(0, 4) @ magic.conj().T
        child = magic @ _pair_generator(0, 2) @ magic.conj().T
        expected_child = -np.kron(I2 + Z, np.kron(Z, I2) + np.kron(I2, Z)) / 4
        np.testing.assert_allclose(child, expected_child, atol=ATOL, rtol=0)
        commutator = root @ child - child @ root
        np.testing.assert_allclose(commutator,
                                   -1j * magic @ _pair_generator(2, 4) @ magic.conj().T,
                                   atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(commutator, 2), 1, delta=ATOL)


if __name__ == '__main__':
    unittest.main()
