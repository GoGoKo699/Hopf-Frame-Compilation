"""Bounded full-matrix tests for structured frames, coherent votes and dirty rotations.

The matrices exercise exact algebra and specified perturbations. They do not
establish an asymptotic circuit bound or implement general synthesis.
"""
from itertools import product
import unittest

import numpy as np


def rotation(angle):
    cosine, sine = np.cos(angle), np.sin(angle)
    return np.array([[cosine, -sine], [sine, cosine]])


def complete_frame(width, angle):
    """Literal chronological real rotations on suffix-zero row pairs."""
    result = np.eye(1 << width)
    for depth in range(width):
        bit = 1 << (width - depth - 1)
        for prefix in range(1 << depth):
            start = prefix * 2 * bit
            rows = [start, start + bit]
            result[rows] = rotation(angle(depth, prefix)) @ result[rows]
    return result


class PrefixRotationProfilesTest(unittest.TestCase):
    def test_guarded_prefix_and_noncontiguous_terminal(self):
        rng = np.random.default_rng(941)
        for width in range(3, 7):
            for prefix_bits in range(1, width - 1):
                tail_bits = width - prefix_bits
                top_angles = [rng.uniform(-2, 2, 1 << j) for j in range(prefix_bits)]
                profiles = rng.uniform(-2, 2, (1 << prefix_bits, tail_bits))

                def prescribed(j, prefix):
                    if j < prefix_bits:
                        return top_angles[j][prefix]
                    return profiles[prefix >> (j - prefix_bits), j - prefix_bits]

                expected = complete_frame(width, prescribed)
                for reserve in range(1, tail_bits):
                    upper = tail_bits - reserve
                    # Top frame on the initial prefix when all lower bits vanish.
                    actual = np.eye(1 << width)
                    top_rows = np.arange(1 << prefix_bits) << tail_bits
                    actual[np.ix_(top_rows, top_rows)] = complete_frame(prefix_bits, prescribed)
                    # Independent middle frames on each prefix, guarding the reserve.
                    for prefix, angles in enumerate(profiles):
                        rows = (prefix << tail_bits) + (np.arange(1 << upper) << reserve)
                        middle = complete_frame(upper, lambda j, _: angles[j])
                        actual[rows] = middle @ actual[rows]
                    # The terminal skips the middle register. Its own top layers are zero.
                    def terminal_angle(j, prefix):
                        if j < prefix_bits:
                            return 0
                        return profiles[prefix >> (j - prefix_bits), upper + j - prefix_bits]

                    terminal = complete_frame(prefix_bits + reserve, terminal_angle)
                    for middle_word in range(1 << upper):
                        rows = [(prefix << tail_bits) + (middle_word << reserve) + low
                                for prefix in range(1 << prefix_bits)
                                for low in range(1 << reserve)]
                        actual[rows] = terminal @ actual[rows]
                    np.testing.assert_allclose(actual, expected, atol=2e-13)

    def test_three_basis_vote_gram(self):
        hadamard = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        phase_h = np.diag([1, 1j]) @ hadamard
        identity = np.eye(8)
        bases = [identity, np.kron(np.kron(hadamard, hadamard), hadamard),
                 np.kron(np.kron(phase_h, phase_h), phase_h)]
        for bad in ([0], [3], [2, 5]):
            projectors = [basis[:, bad] @ basis[:, bad].conj().T for basis in bases]
            all_records = np.zeros_like(identity, dtype=complex)
            wrong_records = np.zeros_like(all_records)
            for votes in product((0, 1), repeat=3):
                branch = identity
                for vote, bad_projector in zip(votes, projectors):
                    factor = bad_projector if vote else identity - bad_projector
                    branch = factor @ branch
                gram = branch.conj().T @ branch
                all_records += gram
                if sum(votes) >= 2:
                    wrong_records += gram
            np.testing.assert_allclose(all_records, identity, atol=2e-14)
            mu = len(bad) / np.sqrt(8)
            full_error = 2 * np.sqrt(np.linalg.eigvalsh(wrong_records)[-1])
            self.assertLessEqual(full_error, 4 * (mu + mu**2) + 1e-13)

    def test_one_candidate_rotation_and_private_twirl(self):
        identity = np.eye(2)
        z = np.diag([1, -1])
        x = np.array([[0, 1], [1, 0]])
        for theta in (-0.81, 0.23):
            # Basis order is logical target, candidate flag, candidate target.
            source = np.eye(4)
            source[:2, :2] = rotation(-theta / 4)
            source[2:, 2:] = rotation(1.13)
            flag_reflection = np.diag([1, -1, 1, 1])
            swap = np.eye(8)
            swap[[1, 4]] = swap[[4, 1]]
            goal = np.kron(rotation(theta), np.eye(4))
            expected = goal.copy()
            invalid = [2, 3, 6, 7]
            expected[:, invalid] = np.eye(8)[:, invalid]
            for impurity in (0, 0.017):
                # The perturbation mixes initialized and occupied flag sectors.
                leaky = np.kron(rotation(impurity), identity) @ source
                reflector = np.kron(identity, leaky.T @ flag_reflection @ leaky)
                echo = reflector @ swap @ reflector @ swap
                logical_z = np.kron(z, np.eye(4))
                result = echo @ logical_z @ echo.T @ logical_z
                if impurity == 0:
                    np.testing.assert_allclose(result, expected, atol=2e-14)
                delta = result - goal
                gram = delta.conj().T @ delta
                logical = np.trace(gram.reshape(2, 4, 2, 4), axis1=1, axis2=3) / 4
                paulis = (identity, x, 1j * x @ z, z)
                masks = [np.kron(identity, np.kron(a, b))
                         for a, b in product(paulis, repeat=2)]
                coherent = np.vstack([p.conj().T @ delta @ p for p in masks]) / 4
                np.testing.assert_allclose(coherent.conj().T @ coherent,
                                           np.kron(logical, np.eye(4)), atol=2e-14)

    def test_nested_and_incomparable_defect_planes(self):
        width = 5
        size = 1 << width
        base_angle = np.pi / 4
        changes = {(0, 0): 0.6, (1, 0): -0.45, (2, 2): 0.71, (4, 3): -0.28}
        positive = complete_frame(width, lambda j, p: base_angle)
        literal = complete_frame(width, lambda j, p: base_angle + changes.get((j, p), 0))
        markers = np.eye(size)
        fathers = {(0, 0): markers[0]}
        factors = []
        for depth in range(width):
            for prefix in range(1 << depth):
                father = fathers[depth, prefix]
                mother = markers[(2 * prefix + 1) << (width - depth - 1)]
                fathers[depth + 1, 2 * prefix] = (father - mother) / np.sqrt(2)
                fathers[depth + 1, 2 * prefix + 1] = (father + mother) / np.sqrt(2)
                if (depth, prefix) in changes:
                    plane = np.column_stack([father, mother])
                    relative = np.eye(size) + plane @ (
                        rotation(changes[depth, prefix]) - np.eye(2)) @ plane.T
                    nesting = sum((j, prefix >> (depth - j)) in changes for j in range(depth))
                    factors.append((nesting, relative))
        actual = np.eye(size)
        for _, relative in sorted(factors, key=lambda item: item[0]):
            actual = relative @ actual
        np.testing.assert_allclose(positive @ actual, literal, atol=3e-14)

    def test_local_father_truncation_at_a_patch_boundary(self):
        # Copy a real marker chain onto a fresh dummy root; all descendant
        # mothers coincide. The discarded component contracts at each level.
        width = 6
        basis = np.eye(1 << width)
        true = basis[0]
        for depth in range(2):
            true = (true - basis[1 << (width - depth - 1)]) / np.sqrt(2)
        local = basis[-1]
        initial_difference = true - local
        for offset in range(3):
            depth = 2 + offset
            mother = basis[1 << (width - depth - 1)]
            np.testing.assert_allclose(true - local,
                                       2**(-offset / 2) * initial_difference, atol=1e-14)
            true = (true - mother) / np.sqrt(2)
            local = (local - mother) / np.sqrt(2)
