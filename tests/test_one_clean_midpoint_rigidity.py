"""Finite algebra audits for one-clean midpoint rigidity.

The fixtures check generalized phase layers, exact clean return with
intermediate mixing, a non-Pauli Clifford midpoint, and compressed
commutant lifting. They do not infer any infinite-width gap numerically,
or remove the midpoint hypothesis from the research theorem.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import gcd, lcm
import unittest

import numpy as np

from tests.test_four_layer_crossing_rank import _binary_basis, _pairing
from tests.test_operator_source_compiler import (
    ATOL, H, I2, X, Y, Z, _adjoint, _pauli, _word_matrix,
)


def _pauli_basis(width):
    return [_pauli(width, {wire: factor for wire, factor in enumerate(factors)})
            for factors in product((I2, X, Y, Z), repeat=width)]


def _is_pauli(matrix, basis):
    dimension = matrix.shape[0]
    return any(abs(abs(np.trace(pauli.conj().T @ matrix) / dimension) - 1) < ATOL
               for pauli in basis)


def _assert_clifford(case, matrix, width, basis):
    for wire, axis in product(range(width), (X, Z)):
        pauli = _pauli(width, {wire: axis})
        image = matrix @ pauli @ matrix.conj().T
        case.assertTrue(_is_pauli(image, basis))


def _phase_matrix(values):
    return np.diag(np.exp(1j * np.pi * np.asarray(values) / 4))


def _clean_embedding(width, clean):
    columns = [index for index in range(1 << width) if not (index >> clean) & 1]
    return np.eye(1 << width)[:, columns]


class OneCleanMidpointRigidityTests(unittest.TestCase):
    def test_dependent_axis_phase_layer_maps_paulis_to_cliffords(self):
        # The seven nonzero Z parities span only three dimensions. Their
        # signed phase sum is 4*x0*x1*x2; one extra x0 gives T_0 CCZ.
        weights = {1: 2, 2: 1, 4: 1, 3: -1, 5: -1, 6: -1, 7: 1}
        self.assertEqual(len(_binary_basis(weights)), 3)
        values = [sum(weight * ((mask & x).bit_count() % 2)
                      for mask, weight in weights.items()) % 8 for x in range(8)]
        expected = [((x & 1) + 4 * int(x == 7)) % 8 for x in range(8)]
        self.assertEqual(values, expected)
        phase = _phase_matrix(values)
        basis = _pauli_basis(3)
        first_image = phase @ _pauli(3, {0: X}) @ phase.conj().T
        self.assertFalse(_is_pauli(first_image, basis))
        for wire, axis in product(range(3), (X, Y, Z)):
            pauli = _pauli(3, {wire: axis})
            image = phase @ pauli @ phase.conj().T
            np.testing.assert_allclose(image, image.conj().T, atol=ATOL, rtol=0)
            _assert_clifford(self, image, 3, basis)

    def test_diagonal_clifford_derivative_including_imaginary_pauli_phase(self):
        # Test integer phase differences before using floating matrices.
        # x0 is translated by the Pauli; the remaining cubic is invariant.
        width = 4
        embedding = np.zeros((16, 8), dtype=complex)
        for logical in range(8):
            embedding[2 * logical:2 * logical + 2, logical] = 1 / np.sqrt(2)
        pauli = _pauli(width, {0: Y, 2: Z})
        for constant, rest_mask in product(range(4), (0, 3, 7)):
            q, remainder = [], []
            for x in range(16):
                bits = [(x >> wire) & 1 for wire in range(width)]
                q.append((2 * constant * bits[0]
                          + 4 * bits[0] * ((rest_mask & (x >> 1)).bit_count())) % 8)
                remainder.append((bits[1] + 2 * bits[1] * bits[2]
                                  + 4 * bits[1] * bits[2] * bits[3]) % 8)
            values = [(a + b) % 8 for a, b in zip(q, remainder)]
            for x in range(16):
                expected = (2 * constant + 4 * ((constant % 2) * (x & 1)
                             + (rest_mask & (x >> 1)).bit_count())) % 8
                self.assertEqual((values[x ^ 1] - values[x]) % 8, expected)
                self.assertEqual((q[x ^ 1] - q[x]) % 8, expected)
            word = [('S', 0)] * constant
            for wire in range(1, width):
                if (rest_mask >> (wire - 1)) & 1:
                    word += [('H', wire), ('CX', 0, wire), ('H', wire)]
            clifford = _word_matrix(width, word)
            np.testing.assert_allclose(clifford, _phase_matrix(q), atol=ATOL, rtol=0)
            diagonal = _phase_matrix(values)
            np.testing.assert_allclose(diagonal @ pauli @ diagonal.conj().T,
                                       clifford @ pauli @ clifford.conj().T,
                                       atol=ATOL, rtol=0)
            invariant = diagonal @ clifford.conj().T
            clean_x = _pauli(width, {0: X})
            np.testing.assert_allclose(invariant @ clean_x, clean_x @ invariant,
                                       atol=ATOL, rtol=0)
            descended = embedding.conj().T @ invariant @ embedding
            np.testing.assert_allclose(descended, _phase_matrix(remainder[::2]),
                                       atol=ATOL, rtol=0)

    def test_two_native_layers_return_clean_after_genuine_mixing(self):
        # The binary involution A maps x to x+parity(x)*1111. K=A*S^4
        # centralizes M=T^4 X^4 T^-4. After decoding the X^4 code the
        # occupied action is CCZ times three CZ gates, hence non-Clifford.
        linear_word = [('CX', 3, 2), ('CX', 2, 3), ('CX', 2, 1),
                       ('CX', 2, 0), ('CX', 3, 2), ('CX', 1, 0),
                       ('CX', 0, 3), ('CX', 0, 2), ('CX', 0, 1), ('CX', 1, 0)]
        encode = [('H', 3)] + [('CX', 3, wire) for wire in range(3)]
        layer = [('T', wire) for wire in range(4)]
        clifford_word = [('S', wire) for wire in range(4)] + linear_word
        word = encode + layer + clifford_word + _adjoint(layer) + _adjoint(encode)
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 8)
        expected = np.diag([1 if x.bit_count() != 2 else -1 for x in range(8)])
        basis = _pauli_basis(3)
        self.assertFalse(_is_pauli(expected @ _pauli(3, {0: X}) @ expected, basis))
        for dirty in (0, 1):
            width = 4 + dirty
            unitary = _word_matrix(width, word)
            embedding = _clean_embedding(width, 3)
            projector = embedding @ embedding.conj().T
            midpoint = _word_matrix(width, encode + layer)
            leakage = (np.eye(1 << width) - projector) @ midpoint @ embedding
            self.assertAlmostEqual(np.linalg.norm(leakage, ord=2), 1 / np.sqrt(2), places=12)
            clean_z = _pauli(width, {3: Z})
            np.testing.assert_allclose(unitary @ clean_z, clean_z @ unitary,
                                       atol=ATOL, rtol=0)
            occupied = embedding.conj().T @ unitary @ embedding
            np.testing.assert_allclose(occupied, np.kron(np.eye(1 << dirty), expected),
                                       atol=ATOL, rtol=0)
            phase = _word_matrix(width, layer)
            core = _word_matrix(width, clifford_word)
            original = _pauli(width, {wire: X for wire in range(4)})
            encoded = phase @ original @ phase.conj().T
            np.testing.assert_allclose(core @ encoded, encoded @ core, atol=ATOL, rtol=0)

    def test_four_native_layers_with_non_pauli_clifford_midpoint(self):
        # The first layer uses independent commuting axes X_s Z_a,
        # X_s Z_b, X_s Z_a Z_b. The second is a single X_s rotation.
        # Literal Clifford parity routing makes the two physical layers
        # explicit; no extra initialized wire is used for the routing.
        routing = [('CX', 2, 1), ('CX', 2, 0), ('CX', 1, 2), ('CX', 0, 2)]
        first_layer = [('T', 0), ('T', 1), ('TDG', 2)]
        first_word = ([('H', 2)] + routing + first_layer
                      + _adjoint(routing) + [('TDG', 2), ('H', 2)])
        middle_clifford = [('X', 0), ('CX', 1, 2)]
        second_word = middle_clifford + _adjoint(first_word)
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in first_word), 4)
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in second_word), 4)
        labels = [4 | (1 << 3), 4 | (1 << 4), 4 | (1 << 3) | (1 << 4)]
        self.assertEqual(len(_binary_basis(labels)), 3)
        self.assertTrue(all(_pairing(a, b, 3) == 0 for a, b in product(labels, repeat=2)))
        base_occupied = None
        for dirty in (0, 1):
            width = 3 + dirty
            first = _word_matrix(width, first_word)
            second = _word_matrix(width, second_word)
            complete = second @ first
            clean_z = _pauli(width, {2: Z})
            cz = _word_matrix(width, [('H', 1), ('CX', 0, 1), ('H', 1)])
            midpoint = clean_z @ cz
            np.testing.assert_allclose(first @ clean_z @ first.conj().T, midpoint,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(second.conj().T @ clean_z @ second, midpoint,
                                       atol=ATOL, rtol=0)
            basis = _pauli_basis(width)
            self.assertFalse(_is_pauli(midpoint, basis))
            _assert_clifford(self, midpoint, width, basis)
            embedding = _clean_embedding(width, 2)
            projector = embedding @ embedding.conj().T
            encoded_projector = (np.eye(1 << width) + midpoint) / 2
            np.testing.assert_allclose(first @ projector @ first.conj().T,
                                       encoded_projector, atol=ATOL, rtol=0)
            np.testing.assert_allclose(complete @ clean_z, clean_z @ complete,
                                       atol=ATOL, rtol=0)
            self.assertAlmostEqual(np.linalg.norm((np.eye(1 << width) - projector)
                                                 @ first @ embedding, ord=2), 1, places=12)
            for factors in ({0: X}, {1: Y}, {0: Z, 1: X}):
                pauli = _pauli(width, factors)
                left = first @ pauli @ first.conj().T
                right = second.conj().T @ pauli @ second
                compressed = (left - right) @ encoded_projector
                terminal = complete @ (projector @ pauli) @ complete.conj().T - projector @ pauli
                np.testing.assert_allclose(second @ compressed @ second.conj().T,
                                           terminal, atol=ATOL, rtol=0)
                self.assertAlmostEqual(np.linalg.norm(compressed, ord=2),
                                       np.linalg.norm(terminal, ord=2), places=12)
            occupied = embedding.conj().T @ complete @ embedding
            if dirty == 0:
                base_occupied = occupied
                self.assertGreater(np.linalg.norm(occupied - np.eye(4), ord=2), 1)
            else:
                np.testing.assert_allclose(occupied, np.kron(I2, base_occupied),
                                           atol=ATOL, rtol=0)

    def test_compressed_commutant_lifting_does_not_assert_global_commutation(self):
        # N is diagonal on wire 0 and full matrix algebra on wire 1.
        # The Clifford projector is not a Pauli eigenspace in this basis.
        involution = _pauli(2, {1: H})
        projector = (np.eye(4) + involution) / 2
        diagonal = _word_matrix(2, [('T', 0)])
        lifted = diagonal @ projector
        h_plus = np.array([np.cos(np.pi / 8), np.sin(np.pi / 8)])
        inclusion = np.kron(h_plus[:, None], I2)
        logical = inclusion.conj().T @ diagonal @ inclusion
        np.testing.assert_allclose(inclusion @ inclusion.conj().T, projector,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(lifted, inclusion @ logical @ inclusion.conj().T,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(lifted, projector @ diagonal, atol=ATOL, rtol=0)
        for first, second in product((I2, Z), (I2, X, Y, Z)):
            pauli = _pauli(2, {0: first, 1: second})
            compressed = projector @ pauli @ projector
            np.testing.assert_allclose(lifted @ compressed, compressed @ lifted,
                                       atol=ATOL, rtol=0)
        noncommuting = _pauli(2, {1: X})
        self.assertGreater(np.linalg.norm(involution @ noncommuting - noncommuting @ involution), 1)
        self.assertGreater(np.linalg.norm(lifted @ noncommuting - noncommuting @ lifted), 0.5)
        # Off-diagonal projector entries are Clifford ring units divided
        # by two; their row ratios recover occupied eighth-root phases.
        ratios = []
        for row in range(4):
            column = row ^ 2
            self.assertAlmostEqual(abs(projector[row, column]), 1 / (2 * np.sqrt(2)), places=12)
            ratios.append(lifted[row, column] / projector[row, column])
        np.testing.assert_allclose(np.asarray(ratios) ** 8, 1, atol=ATOL, rtol=0)
        np.testing.assert_allclose(ratios, np.diag(diagonal), atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(lifted - projector, ord=2),
                               np.sqrt(2 - np.sqrt(2)), places=12)

    def test_finite_projective_order_orbits_and_target_polynomial(self):
        # Enumerate actual automorphisms of Q(zeta_lcm(m,8)) fixing
        # zeta_8, then restrict them to Q(zeta_m). This audits the orbit
        # fraction without substituting an intersection-field formula.
        attained_indices = set()
        coefficients = (-3, 5, -5, 3)
        for order in range(1, 65):
            ambient = lcm(order, 8)
            rational_orbit = {a for a in range(order) if gcd(a, order) == 1}
            fixed_orbit = {a % order for a in range(ambient)
                           if gcd(a, ambient) == 1 and a % 8 == 1}
            self.assertTrue(fixed_orbit <= rational_orbit)
            self.assertEqual(len(rational_orbit) % len(fixed_orbit), 0)
            index = len(rational_orbit) // len(fixed_orbit)
            self.assertIn(index, (1, 2, 4))
            attained_indices.add(index)
            if order > 1:
                for exponent in rational_orbit:
                    root = np.exp(2j * np.pi * exponent / order)
                    value = sum(c * root ** power for power, c in enumerate(coefficients))
                    self.assertGreater(abs(value), 1e-6)
                    self.assertLessEqual(abs(value), 16 + ATOL)
        self.assertEqual(attained_indices, {1, 2, 4})
        self.assertEqual(sum(abs(c) for c in coefficients), 16)
        self.assertEqual(sum(power * abs(c) for power, c in enumerate(coefficients)), 24)
        theta = np.arccos(1 / 3) / 2
        target = np.cos(theta) * I2 - 1j * np.sin(theta) * Y
        doubled_target = np.kron(target, target.conj().T)
        annihilator = sum(c * np.linalg.matrix_power(doubled_target, power)
                          for power, c in enumerate(coefficients))
        np.testing.assert_allclose(annihilator, 0, atol=ATOL, rtol=0)
        phased = np.exp(1j * np.sqrt(2)) * H
        doubled = np.kron(phased, phased.conj().T)
        np.testing.assert_allclose(doubled, np.kron(H, H.conj().T), atol=ATOL, rtol=0)
        np.testing.assert_allclose(doubled @ doubled, np.eye(4), atol=ATOL, rtol=0)
        floor = Fraction(1, 48 * 16 ** 3)
        self.assertEqual(floor, Fraction(1, 196608))
        self.assertLess(Fraction(1, 2 ** 18), floor)
        self.assertGreater(Fraction(1, 2 ** 17), floor)

    def test_eight_boundary_recurrence_uses_independent_coefficient_scale(self):
        count, coefficient_exponent, active = 8, 4, 0
        self.assertNotEqual(count, 1 << (coefficient_exponent + 1))
        degree = 7 + 3 * active + active * (active - 1) // 2
        multiplier = 2 * (count ** 2 - 1)
        constant = count * degree + coefficient_exponent + 4 * active + 4 * 3 + 1
        self.assertEqual((multiplier, constant), (126, 73))
        threshold = 0
        for _ in range(count):
            following = multiplier * threshold + constant
            hull_codimension = 2 * (count - 1) * threshold
            log_gap = (-(count + 1) * hull_codimension - count * degree
                       - coefficient_exponent - 4 * active - 2 * 3)
            log_cross_bound = 2 * 3 - following
            self.assertEqual(log_cross_bound, log_gap - 1)
            self.assertGreaterEqual(log_gap - 1, -following)
            threshold = following
        self.assertEqual(threshold, 73 * (126 ** 8 - 1) // 125)
        self.assertEqual(threshold, 37100281773115459)


if __name__ == '__main__':
    unittest.main()
