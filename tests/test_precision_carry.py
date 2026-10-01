"""Literal loader gauges, width bridges, and the surviving mask cost.

These checks preserve the native loaders' scalar phases.  Width bridges
act on every dirty input, not a prepared source state.  The carried body
still contains transformed masks; counting its outer loaders alone is not
a compiler resource theorem.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

import numpy as np

from tests.test_one_clean_compiler import _pauli_t_word
from tests.test_operator_source_compiler import (
    X, Y, Z, _adjoint, _pauli, _sign_word, _source_word, _word_matrix,
)


ATOL = 8e-11


def _star_plane_word(j):
    """Literal p_j=e^(-i*pi/8) exp(i*pi*Y0 Z1...Z(j-1) Xj/8)."""
    return _pauli_t_word({0: 'Y', **{k: 'Z' for k in range(1, j)}, j: 'X'}, -1)


def _short_plane_word(j):
    return _pauli_t_word({0: 'Y', j: 'X'}, -1)


def _star_loader_word(m):
    """V_m=p_1...p_(m-1), using a linear Clifford parity carry."""
    word = [('CX', k, 0) for k in range(1, m - 1)]
    for j in range(m - 1, 0, -1):
        word += _short_plane_word(j)
        if j > 1:
            word += [('CX', j - 1, 0)]
    return word


def _star_bridge_word(larger, smaller):
    """V_smaller^dagger V_larger; never the oppositely ordered product."""
    if larger == smaller:
        return []
    word = [('CX', k, 0) for k in range(1, larger - 1)]
    for j in range(larger - 1, smaller - 1, -1):
        word += _short_plane_word(j)
        if j > smaller:
            word += [('CX', j - 1, 0)]
    word += [('CX', k, 0) for k in range(1, smaller)]
    return word


def _star_weights(m):
    return np.array([2.0 ** (1 - m)] * 2
                    + [2.0 ** (j - m) for j in range(2, m)])


def _gate_counts(word):
    t = sum(gate[0] in ('T', 'TDG') for gate in word)
    return t, len(word) - t


def _middle_mask_word(j):
    """Clifford C H_j C^dagger = p_j^dagger Z_j p_j."""
    controlled_pauli = [('SDG', 0), ('CX', j, 0), ('S', 0)]
    for k in range(1, j):
        controlled_pauli += [('H', k), ('CX', j, k), ('H', k)]
    conjugator = controlled_pauli + [('SDG', j)]
    return _adjoint(conjugator) + [('H', j)] + conjugator


def _address_mask_word(address):
    return [('H', 1), ('CX', address, 1), ('H', 1)]


def _scalar_stage_word(m, flag, address, gauged):
    """An addressed scalar plus a later address change; all ports retained."""
    loader = _star_loader_word(m)
    mask = _address_mask_word(address)
    if gauged:
        # V^dagger P V, with literal inverse and unchanged address throughout.
        mask = loader + mask + _adjoint(loader)
        controlled = [('CX', flag, 0)]
        uncontrolled = [('X', 0)]
    else:
        controlled = _adjoint(loader) + [('CX', flag, 0)] + loader
        uncontrolled = _adjoint(loader) + [('X', 0)] + loader
    negative = [('X', flag)] + controlled + [('X', flag)]
    scalar = ([('H', flag)] + controlled + mask + uncontrolled
              + _adjoint(mask) + negative + [('H', flag)])
    return scalar + [('H', address)]


class PrecisionCarryTests(unittest.TestCase):
    def test_existing_chain_bridge_has_the_full_width_transfer_denominator(self):
        for m in range(2, 6):
            with self.subTest(m=m):
                bridge_word = _source_word(m + 1) + _adjoint(_source_word(m))
                bridge = _word_matrix(m + 1, bridge_word)
                prior = _word_matrix(m + 1, _source_word(m))
                generator = _pauli(m + 1, {m - 1: Y, m: X})
                plane = np.exp(-1j * np.pi / 8) * (
                    np.cos(np.pi / 8) * np.eye(1 << (m + 1))
                    + 1j * np.sin(np.pi / 8) * generator)
                np.testing.assert_allclose(bridge, prior.conj().T @ plane @ prior,
                                           atol=ATOL, rtol=0)
                probe = _pauli(m + 1, {0: Z})
                correlation = np.trace(probe @ bridge @ probe @ bridge.conj().T).real / len(bridge)
                target = 1 - (1 - 1 / np.sqrt(2)) * 2.0 ** (1 - m)
                self.assertAlmostEqual(correlation, target, delta=ATOL)
                self.assertEqual(_gate_counts(bridge_word)[0], 2 * m - 1)
                # The nondivisible constant coefficient is exactly one.
                self.assertAlmostEqual(target * np.sqrt(2) ** (2 * m - 1),
                                       1 + (2 ** (m - 1) - 1) * np.sqrt(2), delta=ATOL)
                for k in range(1, m):
                    local_probe = _pauli(m + 1, {k: Z})
                    local_correlation = np.trace(
                        local_probe @ bridge @ local_probe @ bridge.conj().T).real / len(bridge)
                    local_target = 1 - (1 - 1 / np.sqrt(2)) * 2.0 ** (-(m - k))
                    self.assertAlmostEqual(local_correlation, local_target, delta=ATOL)

    def test_reversed_star_native_phases_source_weights_and_linear_clifford_loader(self):
        for m in range(2, 7):
            with self.subTest(m=m):
                word = _star_loader_word(m)
                literal = _word_matrix(m, word)
                independent = _word_matrix(m, sum(
                    (_star_plane_word(j) for j in range(m - 1, 0, -1)), []))
                np.testing.assert_allclose(literal, independent, atol=ATOL, rtol=0)
                self.assertEqual(_gate_counts(word), (m - 1, 10 * m - 12))
                source = literal @ _pauli(m, {0: X}) @ literal.conj().T
                expected = sum(np.sqrt(weight) * _pauli(
                    m, {**{k: Z for k in range(j)}, j: X})
                    for j, weight in enumerate(_star_weights(m)))
                np.testing.assert_allclose(source, expected, atol=ATOL, rtol=0)
                np.testing.assert_allclose(source @ source, np.eye(1 << m), atol=ATOL, rtol=0)

    def test_all_small_native_width_bridges_and_their_orientation(self):
        for larger in range(3, 7):
            large = _word_matrix(larger, _star_loader_word(larger))
            for smaller in range(2, larger):
                with self.subTest(larger=larger, smaller=smaller):
                    small = _word_matrix(larger, _star_loader_word(smaller))
                    word = _star_bridge_word(larger, smaller)
                    bridge = _word_matrix(larger, word)
                    np.testing.assert_allclose(bridge, small.conj().T @ large, atol=ATOL, rtol=0)
                    self.assertEqual(_gate_counts(word),
                                     (larger - smaller, 8 * (larger - smaller) + 2 * larger - 4))
                    # Physical code conversion has a different order and cannot
                    # inherit this cheap bridge merely by renaming the basis.
                    self.assertGreater(np.linalg.norm(bridge - small @ large.conj().T, 2), .1)

    def test_three_width_full_scalar_words_telescope_after_addresses_change(self):
        widths = (5, 3, 2)
        flag, address, total = 5, 6, 7
        physical, bodies, loaders = [], [], []
        for m in widths:
            physical.append(_word_matrix(total, _scalar_stage_word(m, flag, address, False)))
            bodies.append(_word_matrix(total, _scalar_stage_word(m, flag, address, True)))
            loaders.append(_word_matrix(total, _star_loader_word(m)))
            np.testing.assert_allclose(physical[-1],
                                       loaders[-1] @ bodies[-1] @ loaders[-1].conj().T,
                                       atol=ATOL, rtol=0)
        bridges = [_word_matrix(total, _star_bridge_word(widths[j], widths[j + 1]))
                   for j in range(2)]
        carried = (loaders[2] @ bodies[2] @ bridges[1] @ bodies[1]
                   @ bridges[0] @ bodies[0] @ loaders[0].conj().T)
        ordinary = physical[2] @ physical[1] @ physical[0]
        np.testing.assert_allclose(carried, ordinary, atol=ATOL, rtol=0)
        # Every flag and dirty column is present; the intermediate scalar word
        # does not return the core, and the logical address changes between words.
        probe = _pauli(total, {0: Z})
        self.assertGreater(np.linalg.norm(physical[0] @ probe - probe @ physical[0], 2), .1)
        width_only = (_adjoint(_star_loader_word(widths[0]))
                      + _star_bridge_word(widths[0], widths[1])
                      + _star_bridge_word(widths[1], widths[2])
                      + _star_loader_word(widths[2]))
        self.assertEqual(_gate_counts(width_only)[0], 2 * (max(widths) - 1))

    def test_valid_reversed_masks_keep_the_exact_transformed_mask_denominator(self):
        for m in range(2, 7):
            loader = _word_matrix(m, _star_loader_word(m))
            for j in range(1, m):
                with self.subTest(m=m, j=j):
                    source_mask = _pauli(m, {j: Z})
                    target = loader.conj().T @ source_mask @ loader
                    tail = _star_bridge_word(m, j + 1)
                    word = tail + _middle_mask_word(j) + _adjoint(tail)
                    native = _word_matrix(m, word)
                    np.testing.assert_allclose(native, target, atol=ATOL, rtol=0)
                    self.assertEqual(_gate_counts(word)[0], 2 * (m - 1 - j))
                    probe = _pauli(m, {0: Z})
                    correlation = np.trace(probe @ native @ probe @ native.conj().T).real / len(native)
                    self.assertAlmostEqual(correlation, 1 - 2.0 ** (-(m - 1 - j)), delta=ATOL)

    def test_reversed_digits_and_a_genuine_small_group_coefficient(self):
        for m in range(2, 7):
            weights = _star_weights(m)
            for integer in range((1 << (m - 1)) + 1):
                signs = np.array(_sign_word(integer, m - 1)[::-1])
                coefficient = weights @ (1 - 2 * signs)
                self.assertAlmostEqual(coefficient, 1 - 2 * integer / (1 << (m - 1)), delta=ATOL)
                if coefficient >= 0:
                    self.assertEqual(signs[0], 0)
        m = 6
        loader = _word_matrix(m, _star_loader_word(m))
        probe = _pauli(m, {0: X})
        for height in (1, 2, 3):
            dimension = 1 << height
            terms = 1 << int(np.ceil(np.log2(8 * height + 12)))
            desired = Fraction(1, 9 * height)
            sine = float(desired) / (terms * np.sqrt(dimension))
            angle = np.arcsin(sine)
            allowance = 1 / (4 * height * terms * np.sqrt(dimension))
            self.assertLess(angle, allowance)
            # C=I, one changed root, and all other target words identity.
            residual = np.eye(dimension)
            roots = [0, dimension // 2]
            residual[np.ix_(roots, roots)] = [[np.cos(angle), -sine], [sine, np.cos(angle)]]
            self.assertAlmostEqual(terms * np.sqrt(dimension) * residual[roots[1], roots[0]],
                                   float(desired), delta=ATOL)
            grid = 1 << (m - 1)
            integer = round((1 - float(desired)) * grid / 2)
            signs = _sign_word(integer, m - 1)[::-1]
            encoded = Fraction(1) - Fraction(2 * integer, grid)
            self.assertNotEqual(encoded, desired)
            self.assertLessEqual(abs(encoded - desired), Fraction(5, 2) * Fraction(1, 1 << m))
            mask = _pauli(m, {j: Z for j, bit in enumerate(signs) if bit})
            transformed = loader.conj().T @ mask @ loader
            correlation = np.trace(probe @ transformed @ probe @ transformed).real / len(transformed)
            self.assertAlmostEqual(correlation, float(encoded), delta=ATOL)
            self.assertGreater(correlation, 0)
            self.assertLessEqual(correlation, .25 + ATOL)


if __name__ == '__main__':
    unittest.main()
