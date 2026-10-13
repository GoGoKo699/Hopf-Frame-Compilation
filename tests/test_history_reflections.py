"""History reflections, a stable anchored extension, and oracle sensitivity.

These are bounded full-matrix diagnostics and exact scalar identities, not
native synthesis routines. The query comparison permits only parameter-
independent inter-query operations; it is not a native T-count lower bound.
"""
from fractions import Fraction
import unittest

import numpy as np

from tests.test_depth_balanced_rejection import (
    ATOL, PAULIS, cap_squared, success_squared,
)
from tests.test_joint_history_compression import extraction_fixture, random_unitary


class HistoryReflectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = extraction_fixture()

    def test_compressor_columns_control_entire_reflection_with_arbitrary_dirty_action(self):
        f = self.fixture
        dirty_dimension = 2
        identity_dirty = np.eye(dirty_dimension)
        cap = np.kron(f['cap'], identity_dirty)
        history = np.kron(f['history'], identity_dirty)
        normalized = 8*history
        dimension = f['dimension']*dirty_dimension
        eye = np.eye(dimension, dtype=complex)
        r_cap = eye-2*cap @ cap.conj().T
        r_ideal = eye-2*normalized @ normalized.conj().T
        other_modes = list(range(2*f['size']))+list(range(3*f['size'], f['dimension']))
        indices = [mode*dirty_dimension+b for mode in other_modes
                   for b in range(dirty_dimension)]
        arbitrary = eye.copy()
        arbitrary[np.ix_(indices, indices)] = random_unitary(
            len(indices), np.random.default_rng(8127))
        h_base = arbitrary @ np.kron(f['compressor'], identity_dirty)
        success = np.kron(f['success'], identity_dirty)
        self.assertGreater(np.linalg.norm(h_base @ success-success, 2), 1)
        dirty_z = np.kron(np.eye(f['dimension']), PAULIS[2])
        self.assertGreater(np.linalg.norm(arbitrary @ dirty_z-dirty_z @ arbitrary, 2), 1)

        for angle in (0, 1e-4, .03, .3):
            with self.subTest(angle=angle):
                perturb = np.eye(f['dimension'], dtype=complex)
                rotation = [[np.cos(angle), -np.sin(angle)],
                            [np.sin(angle), np.cos(angle)]]
                for mode in range(f['size']):
                    labels = [mode, 2*f['size']+mode]
                    perturb[np.ix_(labels, labels)] = rotation
                h = np.kron(perturb, identity_dirty) @ h_base
                epsilon = np.linalg.norm(h @ history-cap/8, 2)
                x = h.conj().T @ cap
                source_error = np.linalg.norm(x-normalized, 2)
                self.assertAlmostEqual(source_error, 8*epsilon, delta=ATOL)
                actual = h.conj().T @ r_cap @ h
                np.testing.assert_allclose(actual, eye-2*x @ x.conj().T,
                                           atol=ATOL, rtol=0)
                np.testing.assert_allclose(actual @ actual, eye, atol=ATOL, rtol=0)
                error = np.linalg.norm(actual-r_ideal, 2)
                self.assertLessEqual(error, 2*source_error+ATOL)
                self.assertLessEqual(error, 16*epsilon+ATOL)

    def test_anchored_householder_is_stable_and_swaps_all_history_columns(self):
        f = self.fixture
        x, cap, projector = 8*f['history'], f['cap'], f['projector']
        overlap = cap.conj().T @ x
        diagonal = np.real(np.diag(overlap))
        np.testing.assert_allclose(overlap, np.diag(diagonal), atol=ATOL, rtol=0)
        active = np.array(f['support'])
        inactive = np.array([mode for mode in range(f['size']) if mode not in active])
        np.testing.assert_allclose(diagonal[inactive], 1, atol=ATOL, rtol=0)
        self.assertGreaterEqual(np.min(diagonal[active]), 2.0**(-f['f']['depth_count']))
        self.assertLessEqual(np.max(diagonal[active]), .5)
        denominator_squared = 2*(1-diagonal[active])
        self.assertTrue(np.all(denominator_squared >= 1))
        self.assertTrue(np.all(denominator_squared < 2))
        inverse_root = np.zeros(f['size'])
        inverse_root[active] = 1/np.sqrt(denominator_squared)
        q = (x-cap) @ np.diag(inverse_root)
        np.testing.assert_allclose(q.conj().T @ q, projector, atol=ATOL, rtol=0)
        np.testing.assert_allclose(q[:, inactive], 0, atol=ATOL, rtol=0)
        householder = np.eye(f['dimension'])-2*q @ q.conj().T
        np.testing.assert_allclose(householder.conj().T, householder, atol=ATOL, rtol=0)
        np.testing.assert_allclose(householder @ householder, np.eye(f['dimension']),
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(householder @ x, cap, atol=ATOL, rtol=0)
        np.testing.assert_allclose(householder @ cap, x, atol=ATOL, rtol=0)
        # Both success and baseline are fixed on every occupied column.
        np.testing.assert_allclose(householder[:, :2*f['size']],
                                   np.eye(f['dimension'])[:, :2*f['size']],
                                   atol=ATOL, rtol=0)
        # The statement tensors with arbitrary dirty inputs; no dirty
        # initialization is used in the supplied dense extension.
        np.testing.assert_allclose(np.kron(householder, np.eye(2)) @ np.kron(x, np.eye(2)),
                                   np.kron(cap, np.eye(2)), atol=ATOL, rtol=0)

        # A nontrivial stage-only perturbation preserves the cap overlap A
        # and therefore uses the same bounded normalization operator.
        perturb = np.eye(f['dimension'], dtype=complex)
        labels = list(f['f']['stages'][-1]['bank'][:2])
        angle = .019
        perturb[np.ix_(labels, labels)] = [[np.cos(angle), -np.sin(angle)],
                                           [np.sin(angle), np.cos(angle)]]
        y = perturb @ x
        np.testing.assert_allclose(y.conj().T @ y, np.eye(f['size']), atol=ATOL, rtol=0)
        np.testing.assert_allclose(cap.conj().T @ y, overlap, atol=ATOL, rtol=0)
        q_y = (y-cap) @ np.diag(inverse_root)
        householder_y = np.eye(f['dimension'])-2*q_y @ q_y.conj().T
        np.testing.assert_allclose(householder_y @ y, cap, atol=ATOL, rtol=0)
        source_change = np.linalg.norm(x-y, 2)
        reflection_change = np.linalg.norm(householder-householder_y, 2)
        self.assertGreater(source_change, .001)
        self.assertGreater(reflection_change, .001)
        self.assertLessEqual(reflection_change, 2*source_change+ATOL)

    def test_root_family_reflection_is_less_sensitive_than_its_source(self):
        for depth_count in (2, 3, 5, 10):
            for axis_position in range(3):
                with self.subTest(depth_count=depth_count, axis_position=axis_position):
                    a = np.sqrt(axis_position+2)*2.0**(-depth_count)
                    b = np.sqrt(1-a*a)
                    theta = 4.0**(-depth_count)/1024
                    residual = np.cos(theta)*np.eye(2)-1j*np.sin(theta)*PAULIS[1]
                    # E aggregates cap and the first j+1 rejections; L
                    # aggregates later rejections after the real Y coin.
                    x_theta = np.vstack([a*np.eye(2), b*residual])
                    x_zero = np.vstack([a*np.eye(2), b*np.eye(2)])
                    np.testing.assert_allclose(x_theta.conj().T @ x_theta, np.eye(2),
                                               atol=ATOL, rtol=0)
                    r_theta = np.eye(4)-2*x_theta @ x_theta.conj().T
                    r_zero = np.eye(4)-2*x_zero @ x_zero.conj().T
                    target_error = np.linalg.norm(residual-np.eye(2), 2)
                    source_error = np.linalg.norm(x_theta-x_zero, 2)
                    reflection_error = np.linalg.norm(r_theta-r_zero, 2)
                    self.assertGreaterEqual(target_error, theta/2)
                    # Normalize by target error so exponentially small
                    # values cannot pass an absolute tolerance check.
                    self.assertAlmostEqual(source_error/target_error, b, delta=2e-10)
                    self.assertAlmostEqual(reflection_error/target_error, 2*a*b,
                                           delta=2e-10)

        # Check the same split in the actual n=5, D=2 physical allocation,
        # using both root logical columns and every allocated history bank.
        f = self.fixture
        depth_count = f['f']['depth_count']
        theta = 4.0**(-depth_count)/1024
        residual = np.cos(theta)*np.eye(2)-1j*np.sin(theta)*PAULIS[1]
        np.testing.assert_array_equal(residual.imag, np.zeros((2, 2)))
        roots = (0, f['size']//2)
        x_theta = np.zeros((f['dimension'], 2), dtype=complex)
        x_zero = np.zeros_like(x_theta)
        for column, mode in enumerate(roots):
            x_theta[2*f['size']+mode, column] = 2.0**(-depth_count)
            x_zero[2*f['size']+mode, column] = 2.0**(-depth_count)
        for stage_number, stage in enumerate(f['f']['stages']):
            amplitude = 2.0**(stage['depth']-depth_count)
            rows = [stage['bank'][stage['support'].index(mode)] for mode in roots]
            # X,Y,Z order: root Y is position j=1. Its own rejected
            # block, as well as the preceding X block, still carries I.
            x_theta[rows, :] = amplitude*(np.eye(2) if stage_number <= 1 else residual)
            x_zero[rows, :] = amplitude*np.eye(2)
        np.testing.assert_allclose(x_theta.conj().T @ x_theta, np.eye(2),
                                   atol=ATOL, rtol=0)
        difference = -2*(x_theta @ x_theta.conj().T-x_zero @ x_zero.conj().T)
        a_squared = 3*4.0**(-depth_count)
        target_error = np.linalg.norm(residual-np.eye(2), 2)
        self.assertAlmostEqual(np.linalg.norm(difference, 2)/target_error,
                               2*np.sqrt(a_squared*(1-a_squared)), delta=2e-10)

    def test_exact_chart_weight_and_query_scaling(self):
        for depth_count in (*range(2, 65), 128, 512):
            for axis_position in range(3):
                with self.subTest(depth_count=depth_count, axis_position=axis_position):
                    early_squared = (axis_position+2)*Fraction(4)**(-depth_count)
                    later_squared = ((2-axis_position)*Fraction(4)**(-depth_count)
                                     + 3*sum(Fraction(4)**(depth-depth_count)
                                             for depth in range(1, depth_count)))
                    self.assertEqual(early_squared+later_squared, 1)
                    self.assertGreaterEqual(later_squared, Fraction(1, 4))
                    theta = Fraction(4)**(-depth_count)/1024
                    self.assertEqual(theta, cap_squared(0, depth_count)/16)
                    alpha_squared = (success_squared(0, axis_position+1, depth_count)
                                     / success_squared(0, axis_position, depth_count))
                    # theta <= (1-alpha)/4 at every possible root position,
                    # checked without irrational roots.
                    self.assertLess(alpha_squared, (1-4*theta)**2)
                    ratio_squared = 4*early_squared*later_squared
                    self.assertLessEqual(ratio_squared, 16*Fraction(4)**(-depth_count))
                    # Compressor inverse columns are b*kappa apart. At
                    # epsilon <= b*kappa/32, the query bound is
                    # q >= 1/(4a) >= 2^(D-3), with b cancelling.
                    query_lower_bound_squared = 1/(16*early_squared)
                    self.assertGreaterEqual(query_lower_bound_squared,
                                            Fraction(4)**depth_count/64)
                    self.assertEqual(later_squared/ratio_squared,
                                     Fraction(4)**depth_count/(4*(axis_position+2)))
        for n in (*range(10, 101), 256, 1024):
            size, depth_count = 1 << n, n-8
            # epsilon <= eta/64 and b*kappa >= 2^(-2D-12) imply the
            # needed epsilon <= b*kappa/32 when N >= 2D+11. Compare
            # exponents directly; never form the enormous 2^(-N).
            self.assertGreaterEqual(size, 2*depth_count+11)
            self.assertGreaterEqual(size+6, 2*depth_count+17)
            self.assertEqual(Fraction(2)**(depth_count-3), Fraction(size, 2048))
        for axis_position in range(3):
            for query_count in (1, 3, 9):
                early_squared = (axis_position+2)*Fraction(4)**(-8)
                # At D=8, these fixed query budgets cannot bridge half
                # the oriented-source separation b*kappa: b cancels.
                self.assertLess(query_count**2*4*early_squared, Fraction(1, 4))


if __name__ == '__main__':
    unittest.main()
