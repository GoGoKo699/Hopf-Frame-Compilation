"""Exact finite witnesses for the returned-dirty compression proof.

These tests audit literal phases, projected-Clifford attenuation and the
error-Gram identity. They do not prove the all-width compression theorem.
All matrix and threshold comparisons use rational Q(sqrt(2), i) arithmetic.
"""
from fractions import Fraction
import unittest

import numpy as np

from verification.fault_tolerant.exact_arithmetic import CQ2, Q2


def _matrix(rows):
    return np.array([[CQ2.parse(value) for value in row] for row in rows], dtype=object)


def _identity(size):
    return _matrix(np.eye(size, dtype=int))


def _adjoint(matrix):
    return _matrix([[CQ2.parse(value).conj() for value in row] for row in matrix.T])


def _gram(matrix):
    return _adjoint(matrix) @ matrix


def _commutator(left, right):
    return left @ right @ _adjoint(left) @ _adjoint(right)


class JointNativeWordTests(unittest.TestCase):
    def setUp(self):
        self.i2 = _identity(2)
        self.halfroot = Q2(0, Fraction(1, 2))
        self.omega = CQ2(self.halfroot, self.halfroot)
        self.h = _matrix([[self.halfroot, self.halfroot],
                          [self.halfroot, -self.halfroot]])
        self.s = _matrix([[1, 0], [0, CQ2(0, 1)]])
        self.t = _matrix([[1, 0], [0, self.omega]])
        self.zero = _matrix([[1], [0]])

    def test_signed_identity_axes_keep_the_literal_native_phase(self):
        def native_t(pauli):
            return ((self.i2 + pauli) * Fraction(1, 2)
                    + (self.i2 - pauli) * (self.omega * Fraction(1, 2)))

        np.testing.assert_array_equal(native_t(self.i2), self.i2)
        np.testing.assert_array_equal(native_t(-self.i2), self.i2 * self.omega)
        hs = self.h @ self.s
        np.testing.assert_array_equal(hs @ hs @ hs, self.i2 * self.omega)

        # Include Y and negative axes, where dropping native scalar phases
        # or a conjugation sign would change the literal matrix equality.
        for pauli in (_matrix([[0, 1], [1, 0]]),
                      _matrix([[0, CQ2(0, -1)], [CQ2(0, 1), 0]]),
                      _matrix([[-1, 0], [0, 1]])):
            gate = native_t(pauli)
            np.testing.assert_array_equal(_gram(gate), self.i2)
            np.testing.assert_array_equal(_adjoint(gate @ gate) @ gate,
                                          _adjoint(gate))

    def test_projected_clifford_attenuation_and_full_isometry_differ(self):
        # Active wire first, removed wire last. Projecting H on the removed
        # wire gives B=I/sqrt(2). The accepted-block distance attains the
        # attenuation lower bound delta_0, whereas full leakage costs more.
        insert = np.kron(self.i2, self.zero)
        clifford = np.kron(self.i2, self.h)
        block = _adjoint(insert) @ clifford @ insert
        np.testing.assert_array_equal(block, self.i2 * self.halfroot)
        np.testing.assert_array_equal(_gram(block), self.i2 * Fraction(1, 2))
        delta_squared = (1 - self.halfroot) * (1 - self.halfroot)
        np.testing.assert_array_equal(_gram(block - self.i2), self.i2 * delta_squared)
        full_squared = Q2(2, -1)
        np.testing.assert_array_equal(_gram(clifford @ insert - insert),
                                      self.i2 * full_squared)
        self.assertGreater(full_squared, delta_squared)

        # A zero block and a two-H block cover cancellation and h=2 without
        # numerical singular-value thresholds.
        x = _matrix([[0, 1], [1, 0]])
        vanished = _adjoint(insert) @ np.kron(self.i2, x) @ insert
        np.testing.assert_array_equal(vanished, self.i2 * 0)
        insert_two = np.kron(insert, self.zero)
        two_h = np.kron(clifford, self.h)
        block_two = _adjoint(insert_two) @ two_h @ insert_two
        np.testing.assert_array_equal(_gram(block_two), self.i2 * Fraction(1, 4))

    def test_source_leakage_preserves_the_complete_error_gram(self):
        # Order: initialized clean wire, system wire, arbitrary dirty wire.
        # CNOT(clean -> dirty) has the non-full projected source Pi=|0><0|_C.
        # Q0 deliberately leaves that source; removing the dirty wire must
        # preserve the complete isometry error, not just its accepted block.
        tx = self.h @ self.t @ self.h
        q_clean = _commutator(tx, _commutator(self.t, tx))
        z = _matrix([[1, 0], [0, -1]])
        q0 = np.kron(q_clean, z)
        i4 = _identity(4)
        i8 = _identity(8)
        clifford = i8[:, [0, 1, 2, 3, 5, 4, 7, 6]]
        insert_removed = np.kron(i4, self.zero)
        insert_clean = np.kron(self.zero, self.i2)
        pi = np.kron(self.zero @ _adjoint(self.zero), self.i2)
        block = _adjoint(insert_removed) @ clifford @ insert_removed
        np.testing.assert_array_equal(block, pi)
        np.testing.assert_array_equal(_gram(block), pi)
        np.testing.assert_array_equal(_gram(q0), i4)

        ideal = insert_clean @ z
        np.testing.assert_array_equal(pi @ ideal, ideal)
        leakage = (i4 - pi) @ q0 @ insert_clean
        self.assertTrue(any(value != CQ2() for value in leakage.flat))

        original_restricted = clifford @ insert_removed @ q0 @ insert_clean
        original_ideal = insert_removed @ ideal
        compressed = q0 @ insert_clean  # F=I is a literal Clifford extension.
        np.testing.assert_array_equal(_adjoint(original_restricted) @ original_ideal,
                                      _adjoint(compressed) @ ideal)
        restricted_error = _gram(original_restricted - original_ideal)
        compressed_error = _gram(compressed - ideal)
        np.testing.assert_array_equal(restricted_error, compressed_error)

        # The original circuit satisfies the theorem's strict small-error
        # premise on EVERY dirty input. Its exact error is (sqrt(2)-1)/2.
        insert_clean_full = np.kron(self.zero, i4)
        original_full = clifford @ np.kron(q0, self.i2) @ insert_clean_full
        ideal_full = insert_clean_full @ np.kron(z, self.i2)
        error_squared = Q2(Fraction(3, 4), Fraction(-1, 2))
        np.testing.assert_array_equal(_gram(original_full - ideal_full),
                                      i4 * error_squared)
        np.testing.assert_array_equal(compressed_error, self.i2 * error_squared)
        self.assertGreater(error_squared, 0)
        self.assertLess(error_squared, (1 - self.halfroot) * (1 - self.halfroot))


if __name__ == "__main__":
    unittest.main()
