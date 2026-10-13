"""Depth-grouped scalar histories and their precision ledger.

Fraction and integer checks are exact. The n=5, D=2 canonical cascade is a
small full-input diagnostic: its bank layout is valid, but it is not the
endpoint allocation D=n-8. No test emits the fine source or supplies a
precision-independent implementation of the prefix-unwinding bank.
"""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.packed_rejection_banks import (
    active_modes, lift_supported_dilation,
)
from tests.test_balanced_rejection_history import (
    ATOL, LAMBDA, bank_modes, canonical_rejection_route,
)


HADAMARD = np.array([[1, 1], [1, -1]], dtype=complex)/np.sqrt(2)
PAULIS = (
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.diag([1, -1]).astype(complex),
)


def cap_squared(depth, depth_count):
    return LAMBDA*Fraction(4)**(depth-depth_count)


def success_squared(depth, axis, depth_count):
    return Fraction(17, 64)-(axis+1)*cap_squared(depth, depth_count)


def routed_cap_gate(n, depth, axis, gate):
    """Conjugate a full predicate-controlled gate by a full affine route.

    The stage route inserts sentinel 1 for source f=0. Toggling that
    insertion on e=1 gives the compressor route for source f=1. All e=0
    inputs, including cap, remain fixed. This constructs the entire
    permutation and predicate gate, including their inactive input columns.
    """
    size, dimension = 1 << n, 4 << n
    permutation, _ = canonical_rejection_route(n, depth, axis)
    on_e = np.flatnonzero(np.arange(dimension) & size)
    permutation[on_e] ^= 1 << (depth+1)
    route = np.zeros((dimension, dimension), dtype=complex)
    route[permutation, np.arange(dimension)] = 1
    controlled = np.eye(dimension, dtype=complex)
    for mode in active_modes(n, depth):
        labels = [2*size+mode, 3*size+mode]
        controlled[np.ix_(labels, labels)] = gate
    return route @ controlled @ route.conj().T, permutation


def routed_depth_walsh(n, depth):
    """H tensor H using cap-to-bank H and SWAP, on complete inputs.

    Bank-to-bank pairs are obtained by swapping one bank with cap,
    applying a routed cap pair, and swapping back. The two butterfly
    layers act on (C,F0),(F1,F2), then (C,F1),(F0,F2).
    """
    pair_h = [routed_cap_gate(n, depth, axis, HADAMARD)[0]
              for axis in range(3)]
    pair_swap = [routed_cap_gate(n, depth, axis, PAULIS[0])[0]
                 for axis in range(2)]
    h12 = pair_swap[1] @ pair_h[2] @ pair_swap[1]
    h02 = pair_swap[0] @ pair_h[2] @ pair_swap[0]
    return h02 @ pair_h[1] @ h12 @ pair_h[0]


def canonical_depth_fixture(rounding_bits=19):
    """Unequal complex coins and nontrivial success-only coarse interlayers.

    These generic canonical angles need not lie in the native near-unit
    source chart. Dyadic entry rounding tests the full canonical comparison
    algebra, not realizability of these particular rows by source masks.
    """
    n, depth_count = 5, 2
    size, dimension = 1 << n, 4 << n
    supports = [active_modes(n, depth) for depth in range(depth_count)]
    entry = np.full(size, depth_count, dtype=int)
    for depth, support in enumerate(supports):
        for mode in support:
            entry[mode] = min(entry[mode], depth)

    constant = np.eye(4, dtype=complex)
    a, b = np.sqrt(17)/8, np.sqrt(47)/8
    constant[np.ix_([0, 1], [0, 1])] = [[a, -b], [b, a]]
    preweight = np.kron(constant, np.eye(size))
    split = np.eye(dimension, dtype=complex)
    for mode, first in enumerate(entry):
        cosine = np.sqrt(float(success_squared(int(first), 0, depth_count)))/a
        sine = np.sqrt(float(cap_squared(int(first), depth_count)))/a
        labels = [mode, 2*size+mode]
        split[np.ix_(labels, labels)] = [[cosine, -sine], [sine, cosine]]
    preweight = split @ preweight

    target = np.eye(size, dtype=complex)
    history = np.eye(dimension, dtype=complex)
    rounded_history = np.eye(dimension, dtype=complex)
    unwind = np.eye(dimension, dtype=complex)
    compression = np.eye(dimension, dtype=complex)
    stages, occupied = [], []
    for depth, support in enumerate(supports):
        rank = len(support)
        for axis in range(3):
            # These coarse words act only on success, preserve P_d, and
            # must occur in the prefix inverse used on this stage's bank.
            coarse = np.eye(rank, dtype=complex)
            logical = np.eye(rank, dtype=complex)
            accepted = np.zeros((rank, rank), dtype=complex)
            defects = np.zeros((rank, rank), dtype=complex)
            alpha = np.sqrt(float(success_squared(depth, axis+1, depth_count)
                                  / success_squared(depth, axis, depth_count)))
            leak = np.sqrt(float(cap_squared(depth, depth_count)
                                 / success_squared(depth, axis, depth_count)))
            radii = []
            for prefix in range(1 << depth):
                sl = slice(2*prefix, 2*prefix+2)
                beta = .04*(axis+1)+.08*prefix+.06*depth
                coarse[sl, sl] = (np.cos(beta)*np.eye(2)
                                  - 1j*np.sin(beta)*PAULIS[(axis+1) % 3])
                theta = .07*(axis+1)+.11*prefix+.03*depth
                logical[sl, sl] = (np.cos(theta)*np.eye(2)
                                   - 1j*np.sin(theta)*PAULIS[axis])
                scale = 1 << rounding_bits
                s = round(alpha*np.cos(theta)*scale)/scale
                t = round(alpha*np.sin(theta)*scale)/scale
                rho_squared = s*s+t*t
                radii.append(np.sqrt(rho_squared))
                accepted[sl, sl] = s*np.eye(2)-1j*t*PAULIS[axis]
                defects[sl, sl] = np.sqrt(1-rho_squared)*np.eye(2)

            coarse_full = np.eye(dimension, dtype=complex)
            coarse_full[np.ix_(support, support)] = coarse
            history = coarse_full @ history
            rounded_history = coarse_full @ rounded_history
            target = coarse_full[:size, :size] @ target
            bank = bank_modes(n, depth, axis)
            occupied.extend(bank)
            unwind[np.ix_(bank, bank)] = target.conj().T[np.ix_(support, support)]
            local = np.block([[alpha*logical, -leak*np.eye(rank)],
                              [leak*np.eye(rank), alpha*logical.conj().T]])
            rounded_local = np.block([[accepted, -defects],
                                      [defects, accepted.conj().T]])
            history = lift_supported_dilation(local, support, bank, dimension) @ history
            rounded_history = (lift_supported_dilation(rounded_local, support, bank,
                                                       dimension) @ rounded_history)
            next_target = np.eye(size, dtype=complex)
            next_target[np.ix_(support, support)] = logical
            target = next_target @ target
            stages.append(dict(depth=depth, axis=axis, support=support, bank=bank,
                               local=local, rounded_local=rounded_local, radii=radii,
                               accepted_error=np.linalg.norm(accepted-alpha*logical, 2),
                               leak=leak))
        compression = routed_depth_walsh(n, depth) @ compression

    desired = np.vstack([target/2, np.sqrt(47)*np.eye(size)/8,
                         np.eye(size)/8, np.zeros((size, size))])
    return dict(n=n, depth_count=depth_count, size=size, dimension=dimension,
                preweight=preweight, history=history, rounded_history=rounded_history,
                unwind=unwind, compression=compression, desired=desired, target=target,
                stages=stages, occupied=occupied, rounding_bits=rounding_bits)


class DepthBalancedRejectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = canonical_depth_fixture()

    def test_fraction_weights_all_entry_depths_and_all_axes(self):
        for depth_count in (1, 2, 3, 8, 19, 40):
            for entry in range(depth_count+1):
                with self.subTest(depth_count=depth_count, entry=entry):
                    cap = cap_squared(entry, depth_count)
                    success = success_squared(entry, 0, depth_count)
                    original_success = success
                    attenuation = Fraction(1)
                    self.assertEqual(success+Fraction(47, 64)+cap, 1)
                    for depth in range(entry, depth_count):
                        self.assertEqual(cap, cap_squared(depth, depth_count))
                        for axis in range(3):
                            before = success_squared(depth, axis, depth_count)
                            after = success_squared(depth, axis+1, depth_count)
                            self.assertEqual(success, before)
                            alpha_squared = after/before
                            self.assertEqual(success*(1-alpha_squared), cap)
                            self.assertGreaterEqual(alpha_squared, Fraction(64, 65))
                            # s_d,a >= 2^(d-D)/sqrt(17), squared exactly.
                            self.assertGreaterEqual(1-alpha_squared,
                                                    Fraction(4)**(depth-depth_count)/17)
                            success *= alpha_squared
                            attenuation *= alpha_squared
                        cap *= 4  # Equal four-vector -> amplitude 2c in cap.
                        self.assertEqual(success,
                                         success_squared(depth+1, 0, depth_count))
                    self.assertEqual(success, Fraction(1, 4))
                    self.assertEqual(cap, LAMBDA)
                    self.assertEqual(attenuation*original_success, Fraction(1, 4))
                    self.assertGreaterEqual(attenuation, Fraction(16, 17))

    def test_exact_capacity_in_finite_layout_and_endpoint_reservation(self):
        for n in range(4, 10):
            size = 1 << n
            for depth_count in range(1, n-2):
                occupied = [mode for depth in range(depth_count) for axis in range(3)
                            for mode in bank_modes(n, depth, axis)]
                self.assertEqual(len(occupied), 6*((1 << depth_count)-1))
                self.assertEqual(len(set(occupied)), len(occupied))
                self.assertTrue(all(3*size <= mode < 4*size for mode in occupied))
                self.assertLess(3*size+len(occupied), 4*size)
        for n in (*range(9, 65), 128, 512, 1024, 4096):
            size, depth_count = 1 << n, n-8
            bank_count = 6*((1 << depth_count)-1)
            self.assertEqual(bank_count, 3*size//128-6)
            self.assertLess(bank_count, size)

    def test_complete_pair_routes_and_walsh_from_h_and_swap(self):
        n, size = self.fixture['n'], self.fixture['size']
        dimension = 4*size
        for depth in range(self.fixture['depth_count']):
            support = active_modes(n, depth)
            for axis in range(3):
                actual, permutation = routed_cap_gate(n, depth, axis, HADAMARD)
                np.testing.assert_array_equal(np.sort(permutation), np.arange(dimension))
                fixed = np.flatnonzero((np.arange(dimension) & size) == 0)
                np.testing.assert_array_equal(permutation[fixed], fixed)
                bank = bank_modes(n, depth, axis)
                self.assertEqual(tuple(permutation[3*size+mode] for mode in support), bank)
                expected = np.eye(dimension, dtype=complex)
                for mode, failure in zip(support, bank):
                    labels = [2*size+mode, failure]
                    expected[np.ix_(labels, labels)] = HADAMARD
                np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
            actual = routed_depth_walsh(n, depth)
            expected = np.eye(dimension, dtype=complex)
            banks = [bank_modes(n, depth, axis) for axis in range(3)]
            for row, mode in enumerate(support):
                labels = [2*size+mode]+[bank[row] for bank in banks]
                expected[np.ix_(labels, labels)] = np.kron(HADAMARD, HADAMARD)
            np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
            np.testing.assert_allclose(actual.conj().T @ actual, np.eye(dimension),
                                       atol=ATOL, rtol=0)

    def test_full_cascade_and_each_unwound_rejection_block(self):
        f, size = self.fixture, self.fixture['size']
        unwound = f['unwind'] @ f['history'] @ f['preweight']
        for stage in f['stages']:
            expected = np.zeros((len(stage['support']), size), dtype=complex)
            amplitude = np.sqrt(float(cap_squared(stage['depth'], f['depth_count'])))
            expected[np.arange(len(stage['support'])), stage['support']] = amplitude
            np.testing.assert_allclose(unwound[stage['bank'], :size], expected,
                                       atol=ATOL, rtol=0)
        complete = f['compression'] @ unwound
        np.testing.assert_allclose(complete.conj().T @ complete, np.eye(f['dimension']),
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(complete[:, :size], f['desired'], atol=ATOL, rtol=0)
        np.testing.assert_allclose(f['desired'].conj().T @ f['desired'], np.eye(size),
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(complete[f['occupied'], :size], 0, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(f['target'].imag), .1)

    def test_omitting_prefix_unwind_is_a_counterexample(self):
        f, size = self.fixture, self.fixture['size']
        wrong = f['compression'] @ f['history'] @ f['preweight']
        np.testing.assert_allclose(wrong[:size, :size], f['target']/2, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(wrong[:, :size]-f['desired'], 2), .01)
        self.assertGreater(np.linalg.norm(wrong[f['occupied'], :size], 2), .001)

    def test_rounded_nonuniform_radii_are_close_but_not_exact_scalar_histories(self):
        f, size = self.fixture, self.fixture['size']
        step = 2.0**(-f['rounding_bits'])
        stage_bound = 0.0
        for stage in f['stages']:
            accepted_error = stage['accepted_error']
            self.assertLessEqual(accepted_error, 2.5*step)
            bound = accepted_error*(1+2/stage['leak'])
            self.assertLessEqual(np.linalg.norm(stage['rounded_local']-stage['local'], 2),
                                 bound+ATOL)
            stage_bound += bound
        self.assertTrue(any(np.ptp(stage['radii']) > 1e-8 for stage in f['stages']))
        rounded = f['compression'] @ f['unwind'] @ f['rounded_history'] @ f['preweight']
        exact = f['compression'] @ f['unwind'] @ f['history'] @ f['preweight']
        np.testing.assert_allclose(rounded.conj().T @ rounded, np.eye(f['dimension']),
                                   atol=ATOL, rtol=0)
        self.assertLessEqual(np.linalg.norm(rounded-exact, 2), stage_bound+ATOL)
        depth_count = f['depth_count']
        ledger = (7.5*depth_count+30*np.sqrt(17)*(2**depth_count-1))*step
        self.assertLessEqual(stage_bound, ledger)
        self.assertLess(np.linalg.norm(rounded[:, :size]-f['desired'], 2), ledger)
        self.assertGreater(np.linalg.norm(rounded[f['occupied'], :size], 2), 1e-8)

    def test_exact_precision_bound_and_literal_width_without_underflow(self):
        # 33/8 strictly bounds sqrt(17). Compare error/eta before evaluating
        # any 2^(-q): q contains N=2^n, which is never used as an exponent.
        self.assertGreater(Fraction(33, 8)**2, 17)
        for n in (*range(9, 101), 128, 512, 1024, 4096):
            size, depth_count = 1 << n, n-8
            q = size+depth_count+12
            inverse_min_defect_sum = 3*sum(1 << k for k in range(1, depth_count+1))
            self.assertEqual(inverse_min_defect_sum, 6*((1 << depth_count)-1))
            coefficient_upper = (Fraction(15, 2)*depth_count
                                 + 30*Fraction(33, 8)*((1 << depth_count)-1))
            self.assertLess(coefficient_upper, 128*(1 << depth_count))
            error_over_eta_upper = coefficient_upper/Fraction(1 << (depth_count+12))
            self.assertLess(error_over_eta_upper, Fraction(1, 32))
            self.assertEqual(q+3, size+n+7)


if __name__ == '__main__':
    unittest.main()
