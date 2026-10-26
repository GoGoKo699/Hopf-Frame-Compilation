"""Native fixtures for a charted, three-call one-clean rotation.

Small q=2,3,4 cases check algebra and every dirty input column; the uniform
precision theorem uses q>=5. These fixtures do not prove sublinear depth or
an improvement to the general frame frontier.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

import numpy as np

from tests.test_one_clean_compiler import (
    _coefficients, _literal_mask, _paired_weights, _sandwich_word,
)
from tests.test_operator_source_compiler import (
    I2, X, Y, Z, _adjoint, _pauli, _rotation, _sign_word, _word_matrix,
)
from tests.test_source_t_depth import _paired_loader
from tests.test_t_depth import _word


ATOL = 8e-11


def _positive_cliffords(word):
    """Keep literal scalars while using only H,S,CX and T/T-dagger."""
    result = []
    for name, *wires in word:
        if name == 'X':
            wire = wires[0]
            result += [('H', wire), ('S', wire), ('S', wire), ('H', wire)]
        elif name == 'Z':
            result += [('S', wires[0])] * 2
        elif name == 'SDG':
            result += [('S', wires[0])] * 3
        else:
            result.append((name, *wires))
    return result


def _chart_signs(q, alpha):
    """Numerical fixture for the proof's certified geometric digit rule."""
    means = (.75 * np.cos(alpha) - .5 * np.sin(alpha) - .5,
             .5 * np.cos(alpha) + np.sin(alpha))
    signs = np.zeros(2 * (q + 1), dtype=int)
    scale = 1 << (q - 1)
    for start, mean in zip((1, 2), means):
        integer = int(np.clip(np.rint((1 - mean) * scale / 2), 0, scale))
        signs[start:2 * q + 1:2] = _sign_word(integer, q - 1)
    return signs


def _seed_scalar(seed, mask, flag):
    inverse = _adjoint(seed)
    source = inverse + [('X', 0)] + seed
    control1 = inverse + [('CX', flag, 0)] + seed
    control0 = inverse + [('X', flag), ('CX', flag, 0), ('X', flag)] + seed
    return ([('H', flag)] + control1 + _adjoint(mask) + source
            + mask + control0 + [('H', flag)])


def _joint_word(q, signs):
    """Qbar=A0^dagger Bbar A0; tail gates occur only in the dressed mask."""
    target, flag = q + 1, q + 2
    schedule = _paired_loader(q)
    seed, tail = _word(schedule[:4]), _word(schedule[4:])
    dressed_mask = tail + _literal_mask(signs) + _adjoint(tail)
    fixed = _seed_scalar(seed, [('Z', 1)], flag)
    programmable = _seed_scalar(seed, dressed_mask, flag)
    cz = [('H', target), ('CX', flag, target), ('H', target)]
    left = cz + fixed + cz
    middle = [('CX', flag, target)] + programmable + [('CX', flag, target)]
    return _positive_cliffords(left + middle + _adjoint(left)), tail


def _three_call_word(word, flag):
    # The two R=-Z factors cancel their signs. The OAA word itself needs -I.
    minus_identity = [('X', flag), ('Z', flag), ('X', flag), ('Z', flag)]
    return _positive_cliffords(
        word + [('Z', flag)] + _adjoint(word) + [('Z', flag)] + word
        + minus_identity)


class JointPrecisionRotationTests(unittest.TestCase):
    def test_chart_coverage_and_available_tail_means(self):
        # Every angle has a nearest chart within pi/8, leaving a strict margin
        # to pi/6. That margin permits terminating certified chart selection.
        self.assertGreater(np.cos(np.pi / 8), np.sqrt(3) / 2)
        for index in range(8):
            center = index * np.pi / 4
            for offset in (-np.pi / 8, -.07, 0., .11, np.pi / 8):
                theta = center + offset
                candidates = [
                    np.arctan2(np.sin(theta - k * np.pi / 4),
                               np.cos(theta - k * np.pi / 4))
                    for k in range(8)
                ]
                alpha = min(candidates, key=abs)
                self.assertLess(abs(alpha), np.pi / 6)
        for alpha in (-np.pi / 6, -.17, 0., .22, np.pi / 6):
            u = .75 * np.cos(alpha) - .5 * np.sin(alpha) - .5
            v = .5 * np.cos(alpha) + np.sin(alpha)
            self.assertLess(abs(u), .4)
            self.assertLessEqual(abs(v), (2 + np.sqrt(3)) / 4 + 1e-15)
            plus, minus = .25 + u / 2, v / 4
            self.assertAlmostEqual(plus + minus, .5 * np.cos(alpha))
            self.assertAlmostEqual(-plus / 2 + 3 * minus / 2,
                                   .5 * np.sin(alpha))

    def test_exact_digit_error_rectangle_and_full_error_constant(self):
        # Exact arithmetic includes both signs of the allowed estimate error.
        for q in (5, 8, 16):
            e = Fraction(1, 1 << (q - 1))
            weights = [Fraction(1, 1 << (j + 1)) for j in range(q - 1)] + [e]
            for mean in (Fraction(-1), Fraction(-7, 19), Fraction(0),
                         Fraction(9, 17), Fraction(1)):
                for shift in (-e / 4, e / 4):
                    estimate = max(Fraction(-1), min(Fraction(1), mean + shift))
                    integer = round((1 - estimate) * (1 << (q - 2)))
                    bits = _sign_word(integer, q - 1)
                    actual = sum(w * (1 - 2 * bit) for w, bit in zip(weights, bits))
                    self.assertLessEqual(abs(actual - mean), 5 * e / 4)
        corners = []
        for a in (Fraction(-5, 4), Fraction(5, 4)):
            for b in (Fraction(-5, 8), Fraction(5, 8)):
                corners.append((a + b) ** 2 + (-a / 2 + 3 * b / 2) ** 2)
        self.assertEqual(max(corners), Fraction(925, 256))
        self.assertLess(13 * max(corners), 49)

    def test_native_dirty_gauge_and_three_call_full_isometry(self):
        for q, alpha in ((2, -.19), (3, .17), (4, -.48)):
            with self.subTest(q=q, alpha=alpha):
                core, target, flag = q + 1, q + 1, q + 2
                width, size = flag + 1, 1 << flag
                signs = _chart_signs(q, alpha)
                word, tail = _joint_word(q, signs)
                actual = _word_matrix(width, word)
                old = _word_matrix(width, _sandwich_word(
                    q, _literal_mask(signs), target, flag))
                gauge = _word_matrix(width, _positive_cliffords(tail))
                np.testing.assert_allclose(actual, gauge.conj().T @ old @ gauge,
                                           atol=ATOL, rtol=0)
                inverse = _word_matrix(width, _positive_cliffords(_adjoint(word)))
                np.testing.assert_allclose(inverse @ actual, np.eye(1 << width),
                                           atol=ATOL, rtol=0)

                weights, fixed = _paired_weights(q)
                s, t = _coefficients(weights, fixed, signs)
                block = np.kron(s * I2 + t * (X @ Z), np.eye(1 << core))
                np.testing.assert_allclose(actual[:size, :size], block,
                                           atol=ATOL, rtol=0)
                amplified = _word_matrix(width, _three_call_word(word, flag))
                radius2 = s * s + t * t
                np.testing.assert_allclose(amplified[:size, :size],
                                           (3 - 4 * radius2) * block,
                                           atol=ATOL, rtol=0)
                embedding = np.eye(1 << width)[:, :size]
                desired = np.kron(_rotation(alpha), np.eye(1 << core))
                error = amplified @ embedding - embedding @ desired
                h = s * np.cos(alpha) + t * np.sin(alpha)
                error2 = 2 - 2 * (3 - 4 * radius2) * h
                np.testing.assert_allclose(error.conj().T @ error,
                                           error2 * np.eye(size), atol=ATOL, rtol=0)
                delta2 = (s - .5 * np.cos(alpha)) ** 2 + (t - .5 * np.sin(alpha)) ** 2
                self.assertAlmostEqual(error2, 8 * (radius2 - .25) ** 2
                                       + 2 * (3 - 4 * radius2) * delta2)

    def test_literal_coarse_clifford_and_actual_inverse(self):
        for k in range(8):
            # Chronological Z,H gives the literal matrix HZ=Ry(pi/4).
            word = _positive_cliffords([('Z', 0), ('H', 0)] * k)
            actual = _word_matrix(1, word)
            np.testing.assert_allclose(actual, _rotation(k * np.pi / 4),
                                       atol=ATOL, rtol=0)
            inverse = _word_matrix(1, _positive_cliffords(_adjoint(word)))
            np.testing.assert_allclose(inverse @ actual, I2, atol=ATOL, rtol=0)

    def test_native_counts_and_paired_tail_depth(self):
        for q in (2, 3, 5, 11, 19):
            word, _ = _joint_word(q, _chart_signs(q, .17))
            amplified = _three_call_word(word, q + 2)
            for native in (amplified, _positive_cliffords(_adjoint(amplified))):
                self.assertTrue(all(gate[0] in ('H', 'S', 'CX', 'T', 'TDG')
                                    for gate in native))
                self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in native),
                                 24 * q + 84)
                self.assertLessEqual(max(wire for gate in native for wire in gate[1:]),
                                     q + 2)
            tail_bands = [gates for kind, gates in _paired_loader(q)[4:] if kind == 'T']
            self.assertEqual(len(tail_bands), q - 1)
            for band in tail_bands:
                targets = [gate[1] for gate in band]
                self.assertEqual(len(targets), len(set(targets)))
                self.assertEqual(len(targets), 2)
            # Twelve tail words and 27 four-T seed source/controlled-source calls.
            self.assertEqual(12 * len(tail_bands) + 27 * 4, 12 * q + 96)

    def test_programmed_precision_bound_in_the_uniform_range(self):
        for q in (5, 8, 16, 24):
            weights, fixed = _paired_weights(q)
            for alpha in (-np.pi / 6, -.31, 0., .17, np.pi / 6):
                s, t = _coefficients(weights, fixed, _chart_signs(q, alpha))
                radius2 = s * s + t * t
                delta2 = (s - .5 * np.cos(alpha)) ** 2 + (t - .5 * np.sin(alpha)) ** 2
                error2 = 8 * (radius2 - .25) ** 2 + 2 * (3 - 4 * radius2) * delta2
                self.assertLessEqual(error2, (7 * 2. ** (-q)) ** 2 + 1e-28)

    def test_third_level_calibration_gap_is_sharp_and_seed_escapes_it(self):
        # M=X_u(T_v CS_uv^dagger) is a Hermitian third-level involution.
        # v-2uv=-u+(u xor v), so the two T gates share one physical layer.
        cubic_word = [('CX', 0, 1), ('TDG', 0), ('T', 1),
                      ('CX', 0, 1), ('X', 0)]
        cubic = _word_matrix(2, _positive_cliffords(cubic_word))
        identity = np.eye(4)
        np.testing.assert_allclose(cubic, cubic.conj().T, atol=ATOL, rtol=0)
        np.testing.assert_allclose(cubic @ cubic, identity, atol=ATOL, rtol=0)

        # Independently check that every Pauli-generator image is Clifford.
        generators = [_pauli(2, {wire: factor})
                      for wire in (0, 1) for factor in (X, Z)]
        paulis = [np.kron(left, right)
                  for left in (I2, X, Y, Z) for right in (I2, X, Y, Z)]
        for generator in generators:
            clifford = cubic @ generator @ cubic.conj().T
            for probe in generators:
                image = clifford @ probe @ clifford.conj().T
                overlap = max(abs(np.trace(pauli.conj().T @ image) / 4)
                              for pauli in paulis)
                self.assertAlmostEqual(overlap, 1.)

        probe = _pauli(2, {1: X})
        product = cubic @ probe @ cubic @ probe
        real_part = (product + product.conj().T) / 2
        np.testing.assert_allclose(real_part, identity / np.sqrt(2),
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(np.linalg.matrix_power(product, 8), identity,
                                   atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(real_part - identity / 2, ord=2),
                               1 / np.sqrt(2) - .5)

        # The fixed native seed has calibration 1/2 instead. Its product has
        # order six, which the proof excludes at every finite hierarchy level.
        seed = _word(_paired_loader(2)[:4])
        seed_source = _word_matrix(2, _positive_cliffords(
            _adjoint(seed) + [('X', 0)] + seed))
        probe = _pauli(2, {1: Z})
        product = seed_source @ probe @ seed_source @ probe
        np.testing.assert_allclose((product + product.conj().T) / 2, identity / 2,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(np.linalg.matrix_power(product, 3), -identity,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(np.linalg.matrix_power(product, 6), identity,
                                   atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(product @ product - identity, ord=2), 1.)

    def test_universal_mask_scalar_closure_is_not_generic_clifford_behavior(self):
        seed = _word(_paired_loader(2)[:4])
        source = _word_matrix(2, _positive_cliffords(
            _adjoint(seed) + [('X', 0)] + seed))
        for left in (I2, X, Y, Z):
            for right in (I2, X, Y, Z):
                mask = np.kron(left, right)
                masked = mask @ source @ mask
                overlap = (source @ masked + masked @ source) / 2
                scalar = np.trace(overlap) / 4
                np.testing.assert_allclose(overlap, scalar * np.eye(4),
                                           atol=ATOL, rtol=0)

        # A Hermitian Clifford source can have commuting Pauli support instead.
        clifford = _word_matrix(2, [('H', 0), ('H', 1)])
        mask = np.kron(Z, Z)
        masked = mask @ clifford @ mask
        overlap = (clifford @ masked + masked @ clifford) / 2
        np.testing.assert_allclose(overlap, -np.kron(Y, Y), atol=ATOL, rtol=0)
        self.assertAlmostEqual(abs(np.trace(overlap)), 0.)
        self.assertAlmostEqual(np.linalg.norm(overlap, ord=2), 1.)


if __name__ == '__main__':
    unittest.main()
