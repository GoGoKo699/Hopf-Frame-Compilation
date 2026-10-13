"""Exact weights and bounded matrix checks for balanced rejection histories.

Fraction identities and integer permutations are exact. Dense matrices are
small floating diagnostics of ideal canonical dilations, not a native fine
compiler. In particular, the prefix-unwinding bank is supplied explicitly;
these checks do not give it a precision-independent synthesis cost.
"""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.conjugate_recovery import recover
from compiler_robust_hopf.packed_rejection_banks import (
    active_modes, lift_supported_dilation,
)


ATOL = 3e-12
LAMBDA = Fraction(1, 64)


def weights(stage_count):
    return [Fraction(1, 4)+LAMBDA*(1-Fraction(2)**(k-stage_count))
            for k in range(stage_count+1)]


def bank_modes(n, depth, axis):
    """11, type[2], zero suffix, sentinel, prefix, target."""
    return tuple(3*(1 << n)+(axis << (n-2))+(1 << (depth+1))+row
                 for row in range(1 << (depth+1)))


def canonical_rejection_route(n, depth, axis):
    """Full controlled-affine truth table, not an emitted gate word.

    The canonical signal is e in (f,e,logical). On e=1, permute (f,logical)
    and apply fixed flips; on e=0, fix every input. The canonical 01 active
    rejection modes become the chosen bank inside 11.
    """
    size = 1 << n
    suffix = n-depth-1
    destinations = [None]*(n+1)  # Relative bit n represents f, excluding e.
    destinations[suffix] = 0
    for bit in range(depth):
        destinations[n-depth+bit] = bit+1
    destinations[n] = depth+1
    destinations[0] = n
    remaining = sorted(set(range(n+1))-set(destinations))
    for bit, destination in zip(range(1, suffix), remaining):
        destinations[bit] = destination
    flips = (1 << n)+(1 << (depth+1))+(axis << (n-2))
    permutation = np.arange(4*size)
    for value in range(2*size):
        moved = sum(((value >> bit) & 1) << destination
                    for bit, destination in enumerate(destinations)) ^ flips
        before = size+(value & (size-1))+((value >> n) << (n+1))
        after = size+(moved & (size-1))+((moved >> n) << (n+1))
        permutation[before] = after
    return permutation, tuple(destinations)


def canonical_fixture():
    """Five logical qubits, six raw stages, unequal complex SU(2) rows."""
    n = 5
    size, dimension = 1 << n, 4 << n
    stages = [(depth, axis) for depth in range(n-3) for axis in range(3)]
    stage_count = len(stages)
    g = weights(stage_count)
    supports = [active_modes(n, depth) for depth, _ in stages]
    entry = np.full(size, stage_count+1, dtype=int)
    for k, support in enumerate(supports, 1):
        for mode in support:
            entry[mode] = min(entry[mode], k)

    # Two complete flag rotations prepare success, baseline, and cap.
    constant = np.eye(4, dtype=complex)
    a, b = np.sqrt(17)/8, np.sqrt(47)/8
    constant[np.ix_([0, 1], [0, 1])] = [[a, -b], [b, a]]
    preweight = np.kron(constant, np.eye(size))
    split = np.eye(dimension, dtype=complex)
    for mode, e in enumerate(entry):
        cap_squared = LAMBDA*Fraction(2)**(int(e)-stage_count-1)
        cosine = np.sqrt(float(g[e-1]))/a
        sine = np.sqrt(float(cap_squared))/a
        labels = [mode, 2*size+mode]
        split[np.ix_(labels, labels)] = [[cosine, -sine], [sine, cosine]]
    preweight = split @ preweight

    paulis = [np.array([[0, 1], [1, 0]], dtype=complex),
              np.array([[0, -1j], [1j, 0]], dtype=complex),
              np.diag([1, -1]).astype(complex)]
    target = np.eye(size, dtype=complex)
    history = np.eye(dimension, dtype=complex)
    unwind = np.eye(dimension, dtype=complex)
    compression = np.eye(dimension, dtype=complex)
    used_banks = []
    for k, ((depth, axis), support) in enumerate(zip(stages, supports), 1):
        rank = len(support)
        logical = np.eye(rank, dtype=complex)
        for prefix in range(1 << depth):
            theta = .07*(axis+1)+.11*prefix+.03*depth
            rotation = np.cos(theta)*np.eye(2)-1j*np.sin(theta)*paulis[axis]
            logical[2*prefix:2*prefix+2, 2*prefix:2*prefix+2] = rotation
        bank = bank_modes(n, depth, axis)
        used_banks.extend(bank)
        # The full prefix is restricted to an invariant later support.
        prefix_inverse = target.conj().T[np.ix_(support, support)]
        unwind[np.ix_(bank, bank)] = prefix_inverse
        alpha = np.sqrt(float(g[k]/g[k-1]))
        leakage = np.sqrt(1-alpha*alpha)
        local = np.block([[alpha*logical, -leakage*np.eye(rank)],
                          [leakage*np.eye(rank), alpha*logical.conj().T]])
        history = lift_supported_dilation(local, support, bank, dimension) @ history
        next_target = np.eye(size, dtype=complex)
        next_target[np.ix_(support, support)] = logical
        target = next_target @ target
        # Conjugated predicate-controlled H, written as its complete matrix.
        hadamards = np.eye(dimension, dtype=complex)
        for mode, failure in zip(support, bank):
            labels = [2*size+mode, failure]
            hadamards[np.ix_(labels, labels)] = np.array([[1, 1], [1, -1]])/np.sqrt(2)
        compression = hadamards @ compression

    desired = np.vstack([target/2, np.sqrt(47)*np.eye(size)/8,
                         np.eye(size)/8, np.zeros((size, size))])
    return dict(size=size, dimension=dimension, target=target, preweight=preweight,
                history=history, unwind=unwind, compression=compression,
                desired=desired, used_banks=used_banks)


class BalancedRejectionHistoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = canonical_fixture()

    def test_exact_fraction_telescoping_for_every_entry_stage(self):
        for count in (1, 2, 3, 6, 12, 30):
            g = weights(count)
            for entry in range(1, count+2):
                with self.subTest(stage_count=count, entry=entry):
                    cap = LAMBDA*Fraction(2)**(entry-count-1)
                    self.assertEqual(g[entry-1]+Fraction(47, 64)+cap, 1)
                    success, attenuation = g[entry-1], Fraction(1)
                    for k in range(entry, count+1):
                        alpha_squared = g[k]/g[k-1]
                        failure = success*(1-alpha_squared)
                        self.assertEqual(failure, LAMBDA*Fraction(2)**(k-count-1))
                        # Equal positive amplitudes are merged by one exact H:
                        # squared cap doubles and the other output is zero.
                        self.assertEqual(cap, failure)
                        cap *= 2
                        success *= alpha_squared
                        attenuation *= alpha_squared
                    self.assertEqual(success, Fraction(1, 4))
                    self.assertEqual(cap, LAMBDA)
                    self.assertEqual(attenuation*g[entry-1], Fraction(1, 4))
                    self.assertGreaterEqual(attenuation, Fraction(16, 17))
            for k in range(1, count+1):
                self.assertGreaterEqual(g[k]/g[k-1], Fraction(32, 33))

    def test_integer_banks_capacity_and_full_controlled_route(self):
        for n in range(4, 9):
            size = 1 << n
            banks = [(depth, axis, bank_modes(n, depth, axis))
                     for depth in range(n-3) for axis in range(3)]
            occupied = [mode for _, _, bank in banks for mode in bank]
            self.assertEqual(len(set(occupied)), len(occupied))
            self.assertEqual(len(occupied), 3*size//4-6)
            self.assertLessEqual(3*size+len(occupied), 4*size)
            self.assertTrue(all(3*size <= mode < 4*size for mode in occupied))
            for depth, axis, bank in banks:
                support = active_modes(n, depth)
                permutation, destinations = canonical_rejection_route(n, depth, axis)
                self.assertEqual(sorted(destinations), list(range(n+1)))
                np.testing.assert_array_equal(np.sort(permutation), np.arange(4*size))
                # Both e=0 sectors, including the whole cap, are fixed.
                fixed = [mode for mode in range(4*size) if not mode & size]
                np.testing.assert_array_equal(permutation[fixed], fixed)
                self.assertEqual(tuple(permutation[size+mode] for mode in support), bank)
                inverse = np.argsort(permutation)
                canonical_active = set(support) | {size+mode for mode in support}
                inactive_banks = ((set(occupied)-set(bank))
                                  | set(range(size, 3*size)))
                self.assertFalse(canonical_active.intersection(inverse[list(inactive_banks)]))

    def test_full_canonical_cascade_unwind_and_balanced_compression(self):
        f = self.fixture
        complete = f['compression'] @ f['unwind'] @ f['history'] @ f['preweight']
        np.testing.assert_allclose(complete.conj().T @ complete, np.eye(f['dimension']),
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(complete[:, :f['size']], f['desired'], atol=ATOL, rtol=0)
        np.testing.assert_allclose(f['desired'].conj().T @ f['desired'], np.eye(f['size']),
                                   atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(f['target'].imag), .1)
        # These column equalities are operator statements, hence also hold
        # after tensoring arbitrary dirty/reference identities.
        np.testing.assert_allclose(complete[f['used_banks'], :f['size']], 0,
                                   atol=ATOL, rtol=0)

    def test_omitting_unwind_breaks_scalar_history_but_not_accepted_block(self):
        f = self.fixture
        wrong = f['compression'] @ f['history'] @ f['preweight']
        size = f['size']
        np.testing.assert_allclose(wrong[:size, :size], f['target']/2, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(wrong[:, :size]-f['desired'], ord=2), .01)
        self.assertGreater(np.linalg.norm(wrong[f['used_banks'], :size], ord=2), .001)
        # Unwinding is needed for the claimed scalar compressor, not for the
        # general OAA lemma, which permits arbitrary actual rejected action.
        expected_recovery = np.vstack([f['target'], np.zeros((3*size, size))])
        np.testing.assert_allclose(recover(wrong, size)[:, :size], expected_recovery,
                                   atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
