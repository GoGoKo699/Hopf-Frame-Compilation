"""Finite complete-unitary checks for the packed-history interface.

Integer routing is exact. Matrix checks use small floating fixtures, including
genuinely dirty rejected blocks; they are not a native fine-circuit emitter.
"""
import unittest

import numpy as np

from compiler_robust_hopf.conjugate_recovery import recover
from compiler_robust_hopf.packed_rejection_banks import (
    active_modes, diagonal_preweight, lift_supported_dilation,
    nested_attenuation, rejection_modes, routing_affine, routing_permutation,
    three_call_error_bound,
)


ATOL = 2e-12


def unitary(rng, dimension):
    a = rng.normal(size=(dimension, dimension))+1j*rng.normal(
        size=(dimension, dimension))
    return np.linalg.qr(a)[0]


def dilation(accepted, left, right):
    """Full Julia dilation with independently dirty rejection coordinates."""
    u, singular, vh = np.linalg.svd(accepted)
    complement = np.sqrt(np.maximum(0, 1-singular**2))
    left_defect = (u*complement) @ u.conj().T
    right_defect = (vh.conj().T*complement) @ vh
    return np.block([[accepted, -left_defect @ right],
                     [left @ right_defect, left @ accepted.conj().T @ right]])


def embed_logical(block, modes, dimension):
    result = np.eye(dimension, dtype=complex)
    result[np.ix_(modes, modes)] = block
    return result


class PackedRejectionBanksTests(unittest.TestCase):
    def test_exact_affine_routing_and_disjoint_banks(self):
        for n in range(1, 8):
            size = 1 << n
            occupied = set()
            for depth in range(n):
                active = active_modes(n, depth)
                bank = rejection_modes(n, depth)
                perm = routing_permutation(n, depth)
                destinations, _ = routing_affine(n, depth)
                self.assertEqual(sorted(destinations), list(range(n+1)))
                self.assertEqual(sorted(perm.tolist()), list(range(4*size)))
                np.testing.assert_array_equal(perm[:2*size], np.arange(2*size))
                self.assertEqual(tuple(perm[2*size+x] for x in active), bank)
                self.assertFalse(occupied.intersection(bank))
                occupied.update(bank)
                inverse = np.argsort(perm)
                np.testing.assert_array_equal(inverse[perm], np.arange(4*size))
            self.assertEqual(len(occupied), 2*size-2)
            self.assertLessEqual(size+size+len(occupied), 4*size)

    def test_three_axis_banks_have_exact_volume_and_native_coordinate_map(self):
        for n in range(3, 9):
            size = 1 << n
            occupied = set()
            for depth in range(n-2):
                for kind in range(3):
                    active = active_modes(n, depth)
                    bank = rejection_modes(n, depth, kind)
                    perm = routing_permutation(n, depth, kind)
                    self.assertEqual(tuple(perm[2*size+x] for x in active), bank)
                    self.assertEqual(len(set(perm.tolist())), 4*size)
                    self.assertFalse(occupied.intersection(bank))
                    occupied.update(bank)
            self.assertEqual(len(occupied), 3*size//2-6)
            self.assertLess(size+size+len(occupied), 4*size)

    def test_routed_full_word_fixes_every_other_bank_on_arbitrary_dirty_inputs(self):
        rng = np.random.default_rng(141010)
        n, dirty = 3, 2
        size, mode_count = 1 << n, 4 << n
        for depth in range(n):
            active = active_modes(n, depth)
            rank = len(active)*dirty
            v = np.kron(unitary(rng, len(active)), np.eye(dirty))
            block = dilation(.93*v, unitary(rng, rank), unitary(rng, rank))
            canonical = lift_supported_dilation(
                block, active, tuple(2*size+x for x in active), mode_count, dirty)
            permutation = routing_permutation(n, depth)
            order = np.array([int(mode)*dirty+j for mode in permutation for j in range(dirty)])
            inverse = np.argsort(order)
            routed = canonical[np.ix_(inverse, inverse)]
            direct = lift_supported_dilation(
                block, active, rejection_modes(n, depth), mode_count, dirty)
            np.testing.assert_allclose(routed, direct, atol=ATOL, rtol=0)
            untouched = set(range(mode_count))-set(active)-set(rejection_modes(n, depth))
            columns = [mode*dirty+j for mode in untouched for j in range(dirty)]
            np.testing.assert_allclose(routed[:, columns], np.eye(mode_count*dirty)[:, columns],
                                       atol=ATOL, rtol=0)

    def test_complete_history_preweight_and_one_literal_recovery(self):
        rng = np.random.default_rng(141011)
        n, dirty = 3, 2
        size, mode_count = 1 << n, 4 << n
        supports = [active_modes(n, d) for d in range(n)]
        alphas = [.94, .95, .96]
        attenuation = nested_attenuation(size, supports, alphas)
        target = np.eye(size, dtype=complex)
        history = np.eye(mode_count*dirty, dtype=complex)
        accepted_product = np.eye(size, dtype=complex)
        stages = []
        for depth, (active, alpha) in enumerate(zip(supports, alphas)):
            rotation = unitary(rng, len(active))
            logical = embed_logical(rotation, active, size)
            target = logical @ target
            diagonal = np.ones(size)
            diagonal[list(active)] = alpha
            accepted_product = logical @ np.diag(diagonal) @ accepted_product
            rank = len(active)*dirty
            block = dilation(alpha*np.kron(rotation, np.eye(dirty)),
                             unitary(rng, rank), unitary(rng, rank))
            stage = lift_supported_dilation(
                block, active, rejection_modes(n, depth), mode_count, dirty)
            stages.append(stage)
            history = stage @ history
        np.testing.assert_allclose(accepted_product, target @ np.diag(attenuation),
                                   atol=ATOL, rtol=0)
        accepted = size*dirty
        np.testing.assert_allclose(history[:accepted, :accepted],
                                   np.kron(accepted_product, np.eye(dirty)), atol=ATOL, rtol=0)
        # The naive scalar attenuation fails on complete frame columns.
        wrong = np.prod(alphas)*target
        self.assertGreater(np.linalg.norm(accepted_product-wrong, ord=2), .09)
        pre = np.kron(np.diag(diagonal_preweight(attenuation)), np.eye(dirty))
        pre_block = dilation(pre, unitary(rng, accepted), unitary(rng, accepted))
        pre_word = lift_supported_dilation(
            pre_block, range(size), range(size, 2*size), mode_count, dirty)
        complete = history @ pre_word
        expected = np.kron(target, np.eye(dirty))
        np.testing.assert_allclose(complete[:accepted, :accepted], expected/2,
                                   atol=ATOL, rtol=0)
        embedding = np.eye(mode_count*dirty)[:, :accepted]
        actual = recover(complete, accepted)
        np.testing.assert_allclose(actual @ embedding, embedding @ expected,
                                   atol=5e-12, rtol=0)
        # This equality is on every dirty input, hence on entangled references.
        self.assertGreater(np.linalg.norm(complete @ embedding-embedding @ expected), 1)
        np.testing.assert_allclose(complete.conj().T @ complete,
                                   np.eye(mode_count*dirty), atol=ATOL, rtol=0)
        reflector = np.diag([1]*accepted+[-1]*(mode_count*dirty-accepted))
        wrong_inverse = -complete @ reflector @ complete @ reflector @ complete
        self.assertGreater(np.linalg.norm(wrong_inverse @ embedding-embedding @ expected), .5)
        # Omitting the leading -I changes the literal answer, despite no leakage.
        np.testing.assert_allclose((-actual) @ embedding, -embedding @ expected,
                                   atol=5e-12, rtol=0)

    def test_bank_collision_reintroduces_rejected_amplitudes(self):
        alpha = .9
        raw = dilation(alpha*np.eye(2), np.eye(2), np.eye(2))
        first = lift_supported_dilation(raw, [0, 1], [2, 3], 6)
        reused = first @ first
        separate = lift_supported_dilation(raw, [0, 1], [4, 5], 6) @ first
        np.testing.assert_allclose(separate[:2, :2], alpha**2*np.eye(2), atol=ATOL)
        np.testing.assert_allclose(reused[:2, :2], (2*alpha**2-1)*np.eye(2), atol=ATOL)
        self.assertGreater(np.linalg.norm(reused[:2, :2]-alpha**2*np.eye(2)), .2)

    def test_robust_recovery_full_error_is_linear_in_accepted_error(self):
        rng = np.random.default_rng(141012)
        dimension = 6  # Three logical modes times a two-dimensional dirty input.
        target = np.kron(unitary(rng, 3), np.eye(2))
        embedding = np.eye(2*dimension)[:, :dimension]
        for epsilon in (.001, .02, .1, .25):
            # A noncommuting polar perturbation tests phase and leakage together.
            rotation = unitary(rng, dimension)
            generator = rotation @ np.diag(np.linspace(-1, 1, dimension)) @ rotation.conj().T
            ev, basis = np.linalg.eigh(generator)
            perturbation = (basis*np.exp(.4j*epsilon*ev)) @ basis.conj().T
            positive = np.diag(.5+epsilon*np.linspace(-.5, .5, dimension))
            accepted = target @ perturbation @ positive
            error = np.linalg.norm(accepted-target/2, ord=2)
            self.assertLessEqual(error, epsilon)
            word = dilation(accepted, unitary(rng, dimension), unitary(rng, dimension))
            amplified = recover(word, dimension)
            actual_error = np.linalg.norm(amplified @ embedding-embedding @ target, ord=2)
            self.assertLessEqual(actual_error, three_call_error_bound(error)+ATOL)
            self.assertLess(three_call_error_bound(error), 8*error)
            # The accepted cubic alone misses the rejected part of this error.
            accepted_error = np.linalg.norm(amplified[:dimension, :dimension]-target, ord=2)
            self.assertGreater(actual_error, accepted_error)

    def test_invalid_inputs_and_normalization(self):
        for args in ((0, 0), (3, 3), (True, 0), (3, 1, 0), (3, 0, 3)):
            with self.assertRaises(ValueError):
                routing_permutation(*args)
        with self.assertRaises(ValueError):
            nested_attenuation(4, [[0, 1], [1, 2]], [.9, .9])
        with self.assertRaises(ValueError):
            diagonal_preweight([.49, .9])
        with self.assertRaises(ValueError):
            lift_supported_dilation(np.eye(4), [0, 1], [1, 2], 4)
        with self.assertRaises(ValueError):
            lift_supported_dilation(.5*np.eye(4), [0, 1], [2, 3], 4)
        with self.assertRaises(ValueError):
            three_call_error_bound(.251)
        for value in ([], 'small', complex(.01), float('nan'), True):
            with self.assertRaises(ValueError):
                three_call_error_bound(value)
        for values in (['.7'], [complex(.7)], [True], [float('nan')]):
            with self.assertRaises(ValueError):
                diagonal_preweight(values)
        self.assertEqual(three_call_error_bound(0), 0)


if __name__ == '__main__':
    unittest.main()
