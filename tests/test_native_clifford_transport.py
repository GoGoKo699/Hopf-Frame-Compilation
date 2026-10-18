"""Finite mechanism checks for the diagonal-third-level transport theorem.

The proof supplies arbitrary-width degrees and the robust radius. These
fixtures check signed quotient reduction, actual inverses, modulo-eight
carries and the endpoint hypothesis; they are not returned-query compilers
or a search over all six-layer sources.
"""
from __future__ import annotations

import itertools
import math
import unittest

import numpy as np

from tests.test_five_layer_native_defect import _anf, _degree
from tests.test_one_band_axis_matching import (
    _one_qubit_cliffords, _projective_distance,
)
from tests.test_outer_axis_normalizer import _image, _tableau, _tau
from tests.test_read_only_three_layer_rigidity import H, S, T, X, Z, _cx, _wire


def _diagonal_phase(phase, width):
    return np.diag([np.exp(1j * math.pi * (phase(z) % 8) / 4)
                    for z in range(1 << width)])


def _bit(value, index):
    return (value >> index) & 1


def _dot(left, right):
    return (left & right).bit_count() & 1


def _fixture_phases(z):
    """An active two-wire phase and a signed Clifford on its idle radical."""
    z0, z1, z2 = [_bit(z, i) for i in range(3)]
    return z0 + 2 * z0 * z1, 4 * z1 + 6 * z2 + 4 * z0 * z2


class NativeCliffordTransport(unittest.TestCase):
    def test_signed_radical_matching_needs_the_z_correction(self):
        width = 3
        h = _diagonal_phase(lambda z: _fixture_phases(z)[0], width)
        g = _diagonal_phase(lambda z: _fixture_phases(z)[1], width)
        source = g @ h
        radical = []
        for a in range(8):
            try:
                image = _image(source, a, width)
            except ValueError:
                continue
            radical.append(a)
            self.assertEqual(image, _image(g, a, width))
            np.testing.assert_allclose(h @ _tau(a, width),
                                       _tau(a, width) @ h, atol=2e-14)
        self.assertEqual(radical, [0, 4])
        self.assertEqual(_image(g, 4, width)[1], 1)
        # Dropping a Z correction preserves the label, but reverses its
        # literal sign and prevents the purported quotient from descending.
        unsigned_g = _wire(Z, 2, width) @ g
        self.assertEqual(_image(unsigned_g, 4, width)[0],
                         _image(g, 4, width)[0])
        self.assertNotEqual(_image(unsigned_g, 4, width)[1],
                            _image(g, 4, width)[1])
        wrong_quotient = unsigned_g.conj().T @ source
        self.assertGreater(np.linalg.norm(wrong_quotient @ _tau(4, width)
                                          - _tau(4, width) @ wrong_quotient), 1)
        # This source is a native Clifford plus a two-T phase word: every
        # sign and column agrees, rather than only its Pauli action.
        native_h = (_wire(S, 0, width) @ _wire(T, 1, width)
                    @ _cx(0, 1, width) @ _wire(T.conj().T, 1, width)
                    @ _cx(0, 1, width))
        np.testing.assert_allclose(h, native_h, atol=2e-14)

    def test_native_two_image_family_has_nonmonomial_quotient_reduction(self):
        width = 3
        h = _diagonal_phase(lambda z: _fixture_phases(z)[0], width)
        g = _diagonal_phase(lambda z: _fixture_phases(z)[1], width)
        source = g @ h
        first_band = _wire(H, 2, width) @ _wire(T, 2, width) @ _wire(H, 2, width)
        second_band = _wire(T, 0, width)
        nonmonomial = nontrivial_correction = False
        for u, v, w in itertools.product((0, 1), repeat=3):
            # Two one-band Pauli images and a linear Pauli, with actual
            # adjoints. The first image is nonmonomial when u=1.
            first = first_band @ _tau(u << 5, width) @ first_band.conj().T
            second = second_band @ _tau(v, width) @ second_band.conj().T
            inner = first @ _tau(w << 1, width) @ second
            matrix, _ = _tableau(inner, width)
            self.assertFalse(np.any(matrix[:2, width:]))  # im B is in R.
            nonmonomial |= bool(np.any(np.sum(abs(inner) > 1e-10, axis=0) > 1))
            quotient_images = [_image(inner.conj().T, 1 << (width + j), width)
                               for j in range(2)]
            self.assertTrue(all(label & 7 == 0 for label, _ in quotient_images))

            def phase(z):
                bits = [sign ^ _dot(label >> width, z)
                        for label, sign in quotient_images]
                return bits[0] + 2 * bits[0] * bits[1] - _fixture_phases(z)[0]

            correction = _diagonal_phase(phase, width)
            np.testing.assert_allclose(correction,
                                       inner.conj().T @ h @ inner @ h.conj().T,
                                       atol=3e-14)
            _tableau(correction, width)
            endpoint = source @ inner @ source.conj().T
            _tableau(endpoint, width)
            np.testing.assert_allclose(endpoint, g @ inner @ correction @ g.conj().T,
                                       atol=4e-14)
            nontrivial_correction |= _projective_distance(correction, np.eye(8)) > 0.1
        self.assertTrue(nonmonomial)
        self.assertTrue(nontrivial_correction)
        # A Hadamard in an active quotient direction violates im B subset R;
        # the transported endpoint is then genuinely non-Clifford.
        with self.assertRaises(ValueError):
            _tableau(source @ _wire(H, 0, width) @ source.conj().T, width)

    def test_mod_eight_carry_can_have_degree_four_h(self):
        for h in (1, 2):
            values = []
            for address in range(1 << (4 * h)):
                monomials = [int(all(_bit(address, h * j + k) for k in range(h)))
                             for j in range(4)]
                values.append((sum(monomials) >> 2) & 1)
            polynomial = _anf(values)
            self.assertEqual(polynomial, {(1 << (4 * h)) - 1})
            self.assertEqual(_degree(polynomial), 4 * h)
            self.assertGreater(_degree(polynomial), 3 * h)
        # Weighted quadratic/cubic terms and subtraction are evaluated as
        # integer phases before taking bits; there is no division in F_2.
        tables = [[], [], []]
        for address in range(256):
            b = [_bit(address, 2 * j) & _bit(address, 2 * j + 1)
                 for j in range(4)]
            phase = (3 * b[0] + 7 * b[1] + b[2] - b[3]
                     + 6 * b[0] * b[1] + 2 * b[2] * b[3]
                     - 4 * b[0] * b[2] * b[3]) % 8
            for bit in range(3):
                tables[bit].append((phase >> bit) & 1)
        for bit, values in enumerate(tables):
            self.assertLessEqual(_degree(_anf(values)), 2 * (1 << bit))

    def test_literal_hermitian_sign_needs_binary_sum_and_integer_dot(self):
        width = 3
        needs_carry = False
        for constant, linear, quadratic in [
            (3, (1, 3, 2), (1, 0, 1)),
            (7, (3, 1, 1), (1, 1, 1)),
        ]:
            def phase(z):
                b = [_bit(z, j) for j in range(width)]
                return (constant + 2 * sum(c * bit for c, bit in zip(linear, b))
                        + 4 * sum(c * b[i] * b[j] for c, (i, j)
                                  in zip(quadratic, itertools.combinations(range(3), 2)))) % 8

            correction = _diagonal_phase(phase, width)
            matrix, _ = _tableau(correction, width)
            shear = matrix[width:, :width]
            for a, b in itertools.product(range(8), repeat=2):
                vector = np.array([_bit(a, j) for j in range(width)])
                c = sum(int(bit) << j for j, bit in enumerate(shear @ vector % 2))
                label, sign = _image(correction, a | (b << width), width)
                self.assertEqual(label, a | ((b ^ c) << width))
                raw_phase = (phase(a) - phase(0)) % 8
                literal_phase = (raw_phase + 2 * ((a & b).bit_count()
                                                - (a & (b ^ c)).bit_count())) % 8
                self.assertEqual(literal_phase % 4, 0)
                self.assertEqual(sign, literal_phase // 4)
                needs_carry |= raw_phase != literal_phase
            # Formula (5) recovers diagonal and off-diagonal shear entries.
            for i in range(width):
                ell = ((phase(1 << i) - phase(0)) % 8) // 2
                self.assertEqual(shear[i, i], ell & 1)
            for i, j in itertools.combinations(range(width), 2):
                four_m = (phase((1 << i) | (1 << j)) - phase(1 << i)
                          - phase(1 << j) + phase(0)) % 8
                self.assertEqual(four_m % 4, 0)
                self.assertEqual(shear[i, j], four_m // 4)
        self.assertTrue(needs_carry)

    def test_third_level_rounding_uses_clifford_endpoints_and_literal_actions(self):
        cliffords = _one_qubit_cliffords()
        exact = nonclifford = 0
        for inner in cliffords:
            endpoint = T @ inner @ T.conj().T
            distances = [_projective_distance(endpoint, candidate) for candidate in cliffords]
            nearest = cliffords[int(np.argmin(distances))]
            if min(distances) < 1e-10:
                exact += 1
                for label in range(4):
                    dressed = T @ _tau(label, 1) @ T.conj().T
                    _tableau(dressed, 1)
                    np.testing.assert_allclose(endpoint @ dressed @ endpoint.conj().T,
                                               nearest @ dressed @ nearest.conj().T, atol=4e-14)
            else:
                nonclifford += 1
                self.assertGreaterEqual(min(distances), math.sqrt(2) / 4 - 1e-12)
        self.assertGreater(exact, 0)
        self.assertGreater(nonclifford, 0)
        # An actual cubic source, beyond a product of physical one-wire T's.
        ccz = _diagonal_phase(lambda z: 4 * _bit(z, 0) * _bit(z, 1) * _bit(z, 2), 3)
        for label in range(64):
            _tableau(ccz @ _tau(label, 3) @ ccz.conj().T, 3)

    def test_endpoint_can_be_clifford_without_intermediate_or_c3_closure(self):
        swap = _cx(0, 1) @ _cx(1, 0) @ _cx(0, 1)
        first, second = _wire(T, 0), _wire(T, 1)
        source = first @ second
        intermediate = second @ swap @ second.conj().T
        with self.assertRaises(ValueError):
            _tableau(intermediate, 2)
        np.testing.assert_allclose(source @ swap @ source.conj().T, swap, atol=2e-14)
        # Native two-band sources need not be third-level. In particular
        # this cannot acquire the theorem's hypothesis by C3 multiplication.
        outside = T @ H @ T
        with self.assertRaises(ValueError):
            _tableau(outside @ X @ outside.conj().T, 1)
        np.testing.assert_allclose(outside.conj().T @ outside, np.eye(2), atol=2e-14)


if __name__ == "__main__":
    unittest.main()
