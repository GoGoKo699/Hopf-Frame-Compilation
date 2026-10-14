"""Finite algebra audits for the four-layer operator-rigidity obstruction.

These check the one-qubit image orbit exactly, literal native words, the
Hermitian-Clifford one-phase lemma, and a nonlocal common commutant. They
do not infer a width-independent theorem from finite enumeration, or claim
that the rare-sector diagonal examples have bounded native T-depth.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations, permutations, product
import unittest

import numpy as np

from tests.exact_matrices import Q2
from tests.test_four_layer_crossing_rank import (
    _binary_basis, _clifford_vector, _pairing,
)
from tests.test_operator_source_compiler import (
    ATOL, X, Y, Z, _adjoint, _pauli, _word_matrix,
)


def _single_qubit_cliffords():
    """The 24 orientation-preserving signed permutations of Bloch axes."""
    result = []
    for permutation in permutations(range(3)):
        parity = (-1) ** sum(permutation[i] > permutation[j]
                            for i in range(3) for j in range(i + 1, 3))
        for signs in product((-1, 1), repeat=3):
            if parity * np.prod(signs) == 1:
                result.append((permutation, signs))
    return result


def _bloch_clifford(vector, clifford):
    permutation, signs = clifford
    return tuple(signs[i] * vector[permutation[i]] for i in range(3))


def _bloch_t(vector):
    x, y, z = vector
    inverse_root_two = Q2(0, Fraction(1, 2))
    return ((x - y) * inverse_root_two,
            (x + y) * inverse_root_two, z)


def _span(generators):
    result = {0}
    for generator in generators:
        result |= {value ^ generator for value in result}
    return result


def _label_matrix(label, width):
    factors = {}
    for wire in range(width):
        x, z = (label >> wire) & 1, (label >> (width + wire)) & 1
        if x or z:
            factors[wire] = Y if x and z else X if x else Z
    return _pauli(width, factors)


class FourLayerOperatorRigidityTests(unittest.TestCase):
    def test_exact_one_qubit_two_layer_image_orbit_and_nearest_pair(self):
        cliffords = _single_qubit_cliffords()
        self.assertEqual(len(cliffords), 24)
        axes = {tuple(Q2(sign if coordinate == axis else 0)
                      for coordinate in range(3))
                for axis in range(3) for sign in (-1, 1)}
        # T-dagger differs from T by a Clifford. Including all Clifford
        # frames and signed input axes therefore covers either gate sign.
        first = {_bloch_clifford(_bloch_t(axis), clifford)
                 for axis in axes for clifford in cliffords}
        images = {_bloch_clifford(_bloch_t(vector), clifford)
                  for vector in first for clifford in cliffords}
        self.assertEqual(len(images), 42)
        self.assertTrue(axes <= images)
        for vector in images:
            self.assertEqual(sum(entry * entry for entry in vector), Q2(1))
        overlaps = [sum(a * b for a, b in zip(left, right))
                    for left, right in combinations(images, 2)]
        nearest_overlap = max(overlaps)
        self.assertEqual(nearest_overlap, Q2(Fraction(1, 4), Fraction(1, 2)))
        squared_distance = 2 * (1 - nearest_overlap)
        self.assertEqual(squared_distance, Q2(Fraction(3, 2), -1))
        distance = 1 - Q2(0, Fraction(1, 2))
        self.assertEqual(squared_distance, distance * distance)
        rho_squared = Q2(Fraction(1, 4), Fraction(-1, 8))
        self.assertGreater(squared_distance, rho_squared)
        # The witness floor tau=rho/2 is strictly above the L=3 error.
        self.assertGreater(rho_squared / 4, Q2(Fraction(1, 64)))

    def test_literal_four_layer_word_split_and_dirty_tensor_invariance(self):
        source = [('T', 0), ('H', 0), ('S', 0), ('H', 0)]
        first_word = source + [('T', 0), ('H', 0)]
        second_word = [('T', 0)] + _adjoint(source)
        for half in (first_word, second_word):
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in half), 2)
        expected_overlap = (1 + 2 * np.sqrt(2)) / 4
        expected_distance = 1 - 1 / np.sqrt(2)
        rho = np.sqrt((1 - 1 / np.sqrt(2)) / 4)
        base_images = None
        for dirty in range(3):
            width = 1 + dirty
            first = _word_matrix(width, first_word)
            second = _word_matrix(width, second_word)
            complete = _word_matrix(width, first_word + second_word)
            pauli = _pauli(width, {0: X})
            left = second.conj().T @ pauli @ second
            right = first @ pauli @ first.conj().T
            np.testing.assert_allclose(complete, second @ first, atol=ATOL, rtol=0)
            defect = pauli - complete @ pauli @ complete.conj().T
            np.testing.assert_allclose(second @ (left - right) @ second.conj().T,
                                       defect, atol=ATOL, rtol=0)
            overlap = np.trace(left @ right) / (1 << width)
            self.assertLess(abs(overlap - expected_overlap), ATOL)
            self.assertAlmostEqual(np.linalg.norm(left - right, ord=2),
                                   expected_distance, places=12)
            self.assertAlmostEqual(np.linalg.norm(defect, ord=2),
                                   expected_distance, places=12)
            self.assertGreater(expected_distance, rho)
            if base_images is None:
                base_images = left, right
            else:
                for image, base in zip((left, right), base_images):
                    np.testing.assert_allclose(image, np.kron(np.eye(1 << dirty), base),
                                               atol=ATOL, rtol=0)

    def test_one_phase_gap_with_entangled_hermitian_clifford_boundaries(self):
        delta = (1 - 1 / np.sqrt(2)) / 2
        maximum_distinct = -1.0
        for width in (2, 3):
            identity = np.eye(1 << width)
            scramble_word = [('H', 0), ('CX', 0, 1), ('S', 1), ('H', 1)]
            if width == 3:
                scramble_word += [('CX', 1, 2), ('H', 2), ('CX', 2, 0)]
            scramble = _word_matrix(width, scramble_word)
            cores = [_word_matrix(width, word) for word in (
                [], [('H', 0)], [('CX', 0, 1)],
                [('H', 1), ('CX', 0, 1), ('H', 1)],
                [('CX', 0, 1), ('CX', 1, 0), ('CX', 0, 1)])]
            cores += [scramble @ core @ scramble.conj().T for core in cores[1:3]]
            cores.append(-cores[1])
            phases = [identity,
                      _word_matrix(width, [('T', 0)]),
                      _word_matrix(width, [('T' if wire % 2 else 'TDG', wire)
                                           for wire in range(width)])]
            phases.append(scramble @ phases[-1] @ scramble.conj().T)
            for core in cores:
                np.testing.assert_allclose(core, core.conj().T, atol=ATOL, rtol=0)
                np.testing.assert_allclose(core @ core, identity, atol=ATOL, rtol=0)
            for phase, left, core in product(phases, cores, cores):
                right = phase @ core @ phase.conj().T
                overlap = np.trace(left @ right) / (1 << width)
                self.assertLess(abs(overlap.imag), ATOL)
                if abs(overlap.real - 1) < ATOL:
                    np.testing.assert_allclose(left, right, atol=ATOL, rtol=0)
                    continue
                maximum_distinct = max(maximum_distinct, overlap.real)
                self.assertLessEqual(overlap.real, 1 - delta + ATOL)
                normalized_squared = np.linalg.norm(left - right, 'fro') ** 2 / (1 << width)
                self.assertAlmostEqual(normalized_squared, 2 * (1 - overlap.real), places=12)
                self.assertGreaterEqual(normalized_squared + ATOL, 2 * delta)
                self.assertGreaterEqual(np.linalg.norm(left - right, ord=2) ** 2 + ATOL,
                                        normalized_squared)
        # H versus T H T-dagger attains the boundary lemma, even with
        # spectator wires. The enumeration includes entangled boundaries too.
        self.assertAlmostEqual(maximum_distinct, 1 - delta, places=12)

    def test_nonmaximal_active_spans_have_exact_common_commutant(self):
        width = 4
        x = [1 << wire for wire in range(width)]
        z = [1 << (width + wire) for wire in range(width)]
        scramble_word = [('H', 0), ('CX', 0, 2), ('S', 2), ('CX', 3, 1),
                         ('H', 1), ('CX', 2, 3), ('S', 0)]
        left = [_clifford_vector(label, width, scramble_word) for label in z[:3]]
        right = [_clifford_vector(label, width, scramble_word)
                 for label in (x[0], z[1], z[2])]
        for space in (left, right):
            self.assertEqual(len(_binary_basis(space)), 3)
            self.assertTrue(all(_pairing(a, b, width) == 0
                                for a, b in product(space, repeat=2)))
        common = _span(left) & _span(right)
        # Use the full centralizers of the actual active spans, including
        # the inactive fourth coordinate; no Lagrangian completion is made.
        centralizers = [{label for label in range(1 << (2 * width))
                         if all(_pairing(label, generator, width) == 0
                                for generator in space)}
                        for space in (left, right)]
        self.assertEqual([len(space) for space in centralizers], [32, 32])
        commuting_generators = _binary_basis(centralizers[0] | centralizers[1])
        commutant = {label for label in range(1 << (2 * width))
                     if all(_pairing(label, generator, width) == 0
                            for generator in commuting_generators)}
        self.assertEqual(commutant, common)
        self.assertEqual(len(common), 4)
        self.assertEqual(len(commuting_generators), 6)
        self.assertEqual(len(_binary_basis(left + right)), 4)

        scramble = _word_matrix(width, scramble_word)
        identity = np.eye(1 << width)
        common_axes = [scramble @ _pauli(width, {wire: Z}) @ scramble.conj().T
                       for wire in (1, 2)]
        projector = (identity - common_axes[0]) @ (identity + common_axes[1]) / 4
        zeta = np.exp(1j * np.pi / 4)
        central = identity + (zeta - 1) * projector
        np.testing.assert_allclose(central.conj().T @ central, identity, atol=ATOL, rtol=0)
        for generator in commuting_generators:
            pauli = _label_matrix(generator, width)
            np.testing.assert_allclose(central @ pauli, pauli @ central, atol=ATOL, rtol=0)
        support = {label for label in range(1 << (2 * width))
                   if abs(np.trace(_label_matrix(label, width) @ central)) > ATOL}
        self.assertEqual(support, common)
        diagonalized = scramble.conj().T @ central @ scramble
        np.testing.assert_allclose(diagonalized, np.diag(np.diag(diagonalized)),
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(np.diag(diagonalized) ** 8, 1, atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(central - identity, ord=2),
                               np.sqrt(2 - np.sqrt(2)), places=12)

    def test_rare_diagonal_reflection_has_small_trace_defect_but_operator_norm_two(self):
        # These are exact multi-controlled-Z matrices, not asserted native
        # four-layer words or a pair of two-layer Pauli images.
        previous_normalized_distance = np.inf
        for width in (2, 4, 6):
            dimension = 1 << width
            reflection = np.eye(dimension)
            reflection[-1, -1] = -1
            difference = reflection - np.eye(dimension)
            normalized_distance = np.linalg.norm(difference, 'fro') / np.sqrt(dimension)
            self.assertAlmostEqual(np.trace(reflection) / dimension, 1 - 2 / dimension,
                                   places=12)
            self.assertAlmostEqual(normalized_distance, 2 / np.sqrt(dimension), places=12)
            self.assertAlmostEqual(np.linalg.norm(difference, ord=2), 2, places=12)
            self.assertLess(normalized_distance, previous_normalized_distance)
            previous_normalized_distance = normalized_distance


if __name__ == '__main__':
    unittest.main()
