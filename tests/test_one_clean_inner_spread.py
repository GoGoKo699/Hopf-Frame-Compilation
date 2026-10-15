"""Finite audits of inner-spread rigidity and its failed common-flip route.

The obstruction to a common unitary flip is checked in exact Q(sqrt(2))
Pauli arithmetic. Literal physical layers independently verify native
realizability. Small spread examples do not prove an infinite-width gap.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

import numpy as np

from tests.test_one_clean_midpoint_rigidity import (
    _assert_clifford, _clean_embedding, _is_pauli, _pauli_basis,
)
from tests.test_operator_source_compiler import (
    ATOL, I2, X, Y, Z, _adjoint, _pauli, _word_matrix,
)


# Exact real coefficients a + b*sqrt(2).
def _q2(a=0, b=0):
    return (Fraction(a), Fraction(b))


def _add(left, right):
    return (left[0] + right[0], left[1] + right[1])


def _mul(left, right):
    return (left[0] * right[0] + 2 * left[1] * right[1],
            left[0] * right[1] + left[1] * right[0])


def _scale(value, scalar):
    return (value[0] * scalar, value[1] * scalar)


_ZERO, _ONE, _INV_SQRT2 = _q2(), _q2(1), _q2(0, Fraction(1, 2))


def _pauli_product(left, right):
    """Return the tensor Pauli label and the power of i multiplying it."""
    label, exponent = [], 0
    cyclic = {('X', 'Y'): 'Z', ('Y', 'Z'): 'X', ('Z', 'X'): 'Y'}
    for first, second in zip(left, right):
        if first == 'I':
            label.append(second)
        elif second == 'I':
            label.append(first)
        elif first == second:
            label.append('I')
        elif (first, second) in cyclic:
            label.append(cyclic[first, second])
            exponent += 1
        else:
            label.append(cyclic[second, first])
            exponent += 3
    return ''.join(label), exponent % 4


def _rotate(expansion, axis, sign):
    """Conjugate exactly by exp(-i*sign*pi*axis/8)."""
    result = {}
    for pauli, coefficient in expansion.items():
        rotated, exponent = _pauli_product(axis, pauli)
        if exponent % 2 == 0:
            result[pauli] = _add(result.get(pauli, _ZERO), coefficient)
        else:
            half = _mul(coefficient, _INV_SQRT2)
            result[pauli] = _add(result.get(pauli, _ZERO), half)
            factor = sign * (1 if exponent == 1 else -1)
            result[rotated] = _add(result.get(rotated, _ZERO),
                                   _scale(half, factor))
    return {pauli: value for pauli, value in result.items() if value != _ZERO}


def _conjugate(expansion, chronological_axes):
    for axis, sign in chronological_axes:
        expansion = _rotate(expansion, axis, sign)
    return expansion


def _multiply(left, right):
    """Multiply real Pauli expansions; values are (real, imag) Q2 pairs."""
    result = {}
    for first, coefficient1 in left.items():
        for second, coefficient2 in right.items():
            pauli, exponent = _pauli_product(first, second)
            value = _mul(coefficient1, coefficient2)
            real, imag = result.get(pauli, (_ZERO, _ZERO))
            signed = _scale(value, 1 if exponent < 2 else -1)
            if exponent % 2:
                imag = _add(imag, signed)
            else:
                real = _add(real, signed)
            result[pauli] = (real, imag)
    return {pauli: value for pauli, value in result.items()
            if value != (_ZERO, _ZERO)}


def _matrix(expansion):
    matrices = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
    width = len(next(iter(expansion)))
    result = np.zeros((1 << width, 1 << width), dtype=complex)
    for pauli, coefficient in expansion.items():
        term = np.ones((1, 1), dtype=complex)
        for factor in pauli:
            term = np.kron(term, matrices[factor])
        result += (float(coefficient[0]) + np.sqrt(2) * float(coefficient[1])) * term
    return result


def _counterexample_expansions():
    # Tensor labels list clean a first, data b second. Initial Clifford
    # CA=(SH)_a H_b maps (Za,Zb) to (YI,IX); CB=(SH)_a to (YI,IZ).
    axes_a = [('ZX', 1), ('YY', 1), ('ZY', -1), ('IY', 1)]
    axes_b = [('ZX', 1), ('ZZ', -1), ('IZ', -1)]
    midpoint_a = _conjugate({'YI': _ONE}, axes_a)
    midpoint_b = _conjugate({'YI': _ONE}, axes_b)
    image_a = _conjugate({'IX': _ONE}, axes_a)
    image_b = _conjugate({'IZ': _ONE}, axes_b)
    return midpoint_a, midpoint_b, image_a, image_b


def _native_counterexample():
    # The repository's chronological gate convention has a=1, b=0.
    clean, data = 1, 0
    initial_a = [('H', clean), ('S', clean), ('H', data)]
    initial_b = [('H', clean), ('S', clean)]
    route_a1 = [('H', clean), ('CX', clean, data), ('S', clean),
                ('SDG', data), ('H', data)]
    route_a2 = [('CX', data, clean), ('H', data), ('S', data)]
    route_b1 = [('CX', clean, data), ('H', data)]
    route_b2 = [('CX', data, clean)]
    layers = [
        (route_a1, [('T', clean), ('T', data)]),
        (route_a2, [('TDG', clean), ('T', data)]),
        (route_b1, [('T', data)]),
        (route_b2, [('TDG', clean), ('TDG', data)]),
    ]
    words = [_adjoint(route) + layer + route for route, layer in layers]
    return initial_a + words[0] + words[1], initial_b + words[2] + words[3], layers


class OneCleanInnerSpreadTests(unittest.TestCase):
    def test_exact_counterexample_and_sector_moment_obstruction(self):
        midpoint, other_midpoint, first, second = _counterexample_expansions()
        self.assertEqual(midpoint, other_midpoint)
        half, quarter_root = _q2(Fraction(1, 2)), _q2(0, Fraction(1, 4))
        self.assertEqual(midpoint, {'YI': half, 'XY': half, 'XX': _scale(half, -1),
                                   'XZ': half})
        self.assertEqual(first, {'IX': quarter_root, 'IZ': _scale(quarter_root, -1),
                                 'ZZ': quarter_root, 'ZX': quarter_root,
                                 'YZ': _scale(half, -1), 'YX': _scale(half, -1)})
        self.assertEqual(second, {'IZ': _INV_SQRT2, 'ZY': _scale(quarter_root, -1),
                                  'ZX': _scale(quarter_root, -1),
                                  'IX': _scale(quarter_root, -1), 'IY': quarter_root})
        product_ab = _multiply(first, second)
        ordinary_trace = product_ab['II']
        signed_trace = (_ZERO, _ZERO)
        for pauli, coefficient in midpoint.items():
            real, imag = product_ab.get(pauli, (_ZERO, _ZERO))
            signed_trace = (_add(signed_trace[0], _mul(coefficient, real)),
                            _add(signed_trace[1], _mul(coefficient, imag)))
        self.assertEqual(ordinary_trace, (_q2(Fraction(-1, 2)), _ZERO))
        self.assertEqual(signed_trace, (_q2(0, Fraction(-1, 4)), _ZERO))
        # A flip with alpha*beta=+1 forces the signed trace to vanish;
        # alpha*beta=-1 forces the ordinary trace to vanish. Both fail.
        self.assertNotEqual(ordinary_trace[0], _ZERO)
        self.assertNotEqual(signed_trace[0], _ZERO)
        sector_plus = _add(ordinary_trace[0], signed_trace[0])
        sector_minus = _add(ordinary_trace[0], _scale(signed_trace[0], -1))
        self.assertNotEqual(sector_plus, sector_minus)
        self.assertNotEqual(sector_plus, _scale(sector_minus, -1))

    def test_literal_native_layers_and_dirty_clean_return(self):
        word_a, word_b, layers = _native_counterexample()
        self.assertEqual([len(layer) for _, layer in layers], [2, 2, 1, 2])
        for _, layer in layers:
            self.assertEqual(len({gate[1] for gate in layer}), len(layer))
            self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in layer))
        midpoint, _, first, second = _counterexample_expansions()
        expected_images = [_matrix(value) for value in (midpoint, first, second)]
        expected_axes = [(('YY', 1), ('ZX', 0)), (('ZY', 1), ('IY', 0)),
                         (('ZX', 0),), (('ZZ', 1), ('IZ', 0))]
        for (routing, _), images in zip(layers, expected_axes):
            route = _word_matrix(2, routing)
            for label, wire in images:
                np.testing.assert_allclose(route @ _pauli(2, {wire: Z}) @ route.conj().T,
                                           _matrix({label: _ONE}), atol=ATOL, rtol=0)
        base_complete = None
        for dirty in (0, 1):
            width = 2 + dirty
            encoder_a, encoder_b = (_word_matrix(width, word) for word in (word_a, word_b))
            clean_z, data_z = _pauli(width, {1: Z}), _pauli(width, {0: Z})
            actual = [encoder_a @ clean_z @ encoder_a.conj().T,
                      encoder_a @ data_z @ encoder_a.conj().T,
                      encoder_b @ data_z @ encoder_b.conj().T]
            for image, expected in zip(actual, expected_images):
                np.testing.assert_allclose(image, np.kron(np.eye(1 << dirty), expected),
                                           atol=ATOL, rtol=0)
            np.testing.assert_allclose(encoder_b @ clean_z @ encoder_b.conj().T,
                                       actual[0], atol=ATOL, rtol=0)
            complete = encoder_b.conj().T @ encoder_a
            if dirty == 0:
                base_complete = complete
            else:
                np.testing.assert_allclose(complete, np.kron(I2, base_complete),
                                           atol=ATOL, rtol=0)
            np.testing.assert_allclose(complete @ clean_z, clean_z @ complete,
                                       atol=ATOL, rtol=0)
            embedding = _clean_embedding(width, 1)
            projector = embedding @ embedding.conj().T
            np.testing.assert_allclose((np.eye(1 << width) - projector) @ complete @ embedding,
                                       0, atol=ATOL, rtol=0)
        # S is a traceless involution but conjugates XI to four Paulis.
        encoded = expected_images[0]
        np.testing.assert_allclose(encoded @ encoded, np.eye(4), atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.trace(encoded), 0)
        conjugated = encoded @ np.kron(X, I2) @ encoded
        four_terms = {'XI': _q2(Fraction(1, 2)), 'YX': _q2(Fraction(-1, 2)),
                      'YY': _q2(Fraction(1, 2)), 'YZ': _q2(Fraction(1, 2))}
        np.testing.assert_allclose(conjugated, _matrix(four_terms), atol=ATOL, rtol=0)
        self.assertFalse(_is_pauli(conjugated, _pauli_basis(2)))

    def test_counterexample_inner_clean_spread_is_one(self):
        first = _conjugate({'YI': _ONE}, [('ZX', 1), ('YY', 1)])
        second = _conjugate({'YI': _ONE}, [('ZX', 1)])
        expected = {'YI': _INV_SQRT2, 'XX': _scale(_INV_SQRT2, -1)}
        self.assertEqual(first, expected)
        self.assertEqual(second, expected)
        # The silent YY gate counts physically, but adds no clean spread.
        self.assertEqual(_pauli_product('YY', 'YI')[1] % 2, 0)
        self.assertEqual(_pauli_product('ZX', 'YI')[1] % 2, 1)

    def test_flat_inner_spread_and_outer_clifford_terms(self):
        width = 3
        routing = [('H', 1), ('CX', 1, 0), ('S', 2), ('CX', 0, 2)]
        basis = _pauli_basis(width)
        change = _word_matrix(width, routing)
        inner = _word_matrix(width, [('T', 0), ('TDG', 1), ('T', 2)])
        outer_route = [('H', 0), ('CX', 2, 0), ('H', 2), ('S', 1)]
        outer = _word_matrix(width, _adjoint(outer_route)
                             + [('TDG', 0), ('T', 1), ('T', 2)] + outer_route)
        for spread in range(width + 1):
            factors = {wire: X for wire in range(spread)}
            if spread < width:
                factors[spread] = Z
            source = _pauli(width, factors)
            core = change @ inner @ source @ inner.conj().T @ change.conj().T
            coefficients = np.array([np.trace(pauli @ core) / (1 << width) for pauli in basis])
            support = np.where(abs(coefficients) > ATOL)[0]
            self.assertEqual(len(support), 1 << spread)
            np.testing.assert_allclose(abs(coefficients[support]), 2 ** (-spread / 2),
                                       atol=ATOL, rtol=0)
            reconstructed = np.zeros_like(core)
            for index in support:
                image = outer @ basis[index] @ outer.conj().T
                _assert_clifford(self, image, width, basis)
                reconstructed += coefficients[index] * image
            np.testing.assert_allclose(reconstructed, outer @ core @ outer.conj().T,
                                       atol=ATOL, rtol=0)


    def test_bounded_spread_boundary_expansion_and_phase_duplication(self):
        word_a, word_b, layers = _native_counterexample()
        encoder_a, encoder_b = (_word_matrix(2, word) for word in (word_a, word_b))
        outer_a = _word_matrix(2, _adjoint(layers[1][0]) + layers[1][1] + layers[1][0])
        outer_b = _word_matrix(2, _adjoint(layers[3][0]) + layers[3][1] + layers[3][0])
        midpoint, _, image_a, image_b = _counterexample_expansions()
        encoded, first, second = map(_matrix, (midpoint, image_a, image_b))
        clean_z, data_z = _pauli(2, {1: Z}), _pauli(2, {0: Z})
        inner_b = outer_b.conj().T @ encoder_b
        core_clean_b = inner_b @ clean_z @ inner_b.conj().T
        core_data_b = inner_b @ data_z @ inner_b.conj().T
        projector = (np.eye(4) + encoded) / 2
        basis = _pauli_basis(2)
        spread_terms = [(1, outer_a @ _matrix({'YI': _ONE}) @ outer_a.conj().T),
                        (-1, outer_a @ _matrix({'XX': _ONE}) @ outer_a.conj().T)]
        np.testing.assert_allclose(sum(sign * term for sign, term in spread_terms) / np.sqrt(2),
                                   encoded, atol=ATOL, rtol=0)
        # sqrt(2) = zeta8 + conjugate(zeta8). Exact real parts turn
        # six coefficients 1/4 into twelve coefficients 2^(-5/2).
        common_coefficient = _q2(0, Fraction(1, 8))
        self.assertEqual(_mul(common_coefficient, _q2(0, 1)), _q2(Fraction(1, 4)))
        zeta = np.exp(1j * np.pi / 4)
        common_float = 2 ** (-5 / 2)
        for label in ('II', 'IY', 'ZI', 'ZY'):
            test_pauli = _matrix({label: _ONE})
            np.testing.assert_allclose(test_pauli @ outer_a, outer_a @ test_pauli,
                                       atol=ATOL, rtol=0)
            left = first @ test_pauli @ first
            sandwich = encoded @ left @ encoded
            _assert_clifford(self, left, 2, basis)
            _assert_clifford(self, sandwich, 2, basis)
            core = core_data_b @ (outer_b.conj().T @ test_pauli @ outer_b) @ core_data_b
            right_terms = [core, core_clean_b @ core, core @ core_clean_b,
                           core_clean_b @ core @ core_clean_b]
            bare_terms = [left, sandwich]
            mixed_terms = [sign * product for sign, term in spread_terms
                           for product in (term @ left, left @ term)]
            self.assertEqual(len(bare_terms) + len(mixed_terms) + len(right_terms),
                             2 + 2 ** (1 + 1) + 4)
            expansion = (sum(bare_terms) / 4 + sum(mixed_terms) / (4 * np.sqrt(2))
                         - outer_b @ sum(right_terms) @ outer_b.conj().T / 4)
            expected = projector @ (left - second @ test_pauli @ second) @ projector
            np.testing.assert_allclose(expansion, expected, atol=ATOL, rtol=0)
            # Left boundaries are bare Cliffords; the right boundaries
            # all retain the same outer phase operation E_B.
            duplicated_left = [phase * term for term in bare_terms for phase in (zeta, zeta.conjugate())]
            duplicated_left += mixed_terms
            duplicated_right = [-phase * term for term in right_terms
                                for phase in (zeta, zeta.conjugate())]
            self.assertEqual(len(duplicated_left) + len(duplicated_right), 16)
            for term in duplicated_left + duplicated_right:
                _assert_clifford(self, term, 2, basis)
            duplicated = common_float * (sum(duplicated_left)
                                          + outer_b @ sum(duplicated_right) @ outer_b.conj().T)
            np.testing.assert_allclose(duplicated, expected, atol=ATOL, rtol=0)
        # The general normalization uses phase duplication, then only
        # opposite-sign pairs to reach the power-of-two boundary cap.
        for spread in range(9):
            raw_count = 6 * 2 ** ((spread + 1) // 2) + 2 ** (spread + 1)
            count = 2 ** (spread + 3)
            self.assertLessEqual(raw_count, count)
            self.assertEqual((count - raw_count) % 2, 0)
            common = _ONE
            for _ in range(spread + 4):
                common = _mul(common, _INV_SQRT2)
            multiplier = (_q2(2 ** (spread // 2)) if spread % 2 == 0
                          else _q2(0, 2 ** (spread // 2)))
            self.assertEqual(_mul(common, multiplier), _q2(Fraction(1, 4)))
        # Finite-boundary normalization: b=t+4; log2(K)=t+3 contributes
        # four further t's to the recurrence constant.
        for spread, expected_constant in ((0, 73), (1, 134)):
            count, exponent = 2 ** (spread + 3), spread + 4
            constant = 7 * count + exponent + 4 * (spread + 3) + 1
            self.assertEqual(constant, expected_constant)
            self.assertEqual(constant, 7 * count + 5 * spread + 17)
            multiplier = 2 * (count ** 2 - 1)
            threshold = 0
            for _ in range(count):
                threshold = multiplier * threshold + constant
            self.assertEqual(threshold, constant * (multiplier ** count - 1) // (multiplier - 1))
            if spread == 0:
                self.assertEqual(threshold, 37100281773115459)

    def test_signed_counterexample_blocks_both_full_norm_completions(self):
        midpoint, _, image_a, image_b = _counterexample_expansions()
        encoded, first, original_second = map(_matrix, (midpoint, image_a, image_b))
        second = -original_second
        occupied = (np.eye(4) + encoded) / 2
        unused = (np.eye(4) - encoded) / 2
        occupied_norm = np.linalg.norm((first - second) @ occupied, ord=2)
        unused_norm = np.linalg.norm((first - second) @ unused, ord=2)
        self.assertAlmostEqual(occupied_norm, np.sqrt(1 - 1 / np.sqrt(2)), places=12)
        self.assertAlmostEqual(unused_norm, np.sqrt(1 + 1 / np.sqrt(2)), places=12)
        self.assertGreater(np.linalg.norm(first - second, ord=2), occupied_norm + 0.1)
        self.assertGreater(np.linalg.norm(first - encoded @ second, ord=2), occupied_norm + 0.1)
        # B' is itself a native two-layer image: prepend data X to UB.
        _, word_b, _ = _native_counterexample()
        encoder = _word_matrix(2, [('X', 0)] + word_b)
        np.testing.assert_allclose(encoder @ _pauli(2, {0: Z}) @ encoder.conj().T,
                                   second, atol=ATOL, rtol=0)
        np.testing.assert_allclose(encoder @ _pauli(2, {1: Z}) @ encoder.conj().T,
                                   encoded, atol=ATOL, rtol=0)

if __name__ == '__main__':
    unittest.main()
