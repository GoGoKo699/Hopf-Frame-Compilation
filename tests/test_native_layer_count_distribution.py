"""Finite layer-count bounds and correlated native Pauli-packet checks."""

from fractions import Fraction
from itertools import combinations, product
from random import Random
import unittest

import numpy as np

from tests.test_clifford_orbit_query_rank import _rank
from tests.test_operator_source_compiler import _adjoint, _word_matrix


def _pair(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width))
            ^ ((right & mask) & (left >> width))).bit_count() % 2


def _gram_rank(axes, width):
    return _rank(sum(_pair(a, b, width) << i for i, b in enumerate(axes))
                 for a in axes)


def _profiles(total, maximum):
    if total == 0:
        yield ()
    for k in range(min(total, maximum), 0, -1):
        for rest in _profiles(total - k, k):
            yield (k, *rest)


class NativeLayerCountDistributionTests(unittest.TestCase):
    def test_native_prefix_packet_packing_and_relative_overlap(self):
        # A 2-by-2 binary matrix, four physical wires. Native gate
        # matrices are compared with the independently derived Pauli expansion
        # and relative-matrix contraction, including every physical column.
        q = 2
        size = 1 << (2 * q)
        matrices = [cols for cols in product(range(1 << q), repeat=q)
                    if _rank(cols) == q]

        def apply(cols, word):
            value = 0
            for j, col in enumerate(cols):
                if word >> j & 1:
                    value ^= col
            return value

        def transpose(cols):
            return tuple(sum(((cols[j] >> i) & 1) << j for j in range(q))
                         for i in range(q))

        def inverse(cols):
            return tuple(next(x for x in range(1 << q)
                              if apply(cols, x) == 1 << j) for j in range(q))

        def section(x, z):
            matrix = np.zeros((size, size), dtype=complex)
            for column in range(size):
                matrix[column ^ x, column] = (-1) ** ((z & column).bit_count())
            return matrix

        probe = section((1 << q) - 1, 0)
        packets = {}
        for cols in matrices:
            frame = [('CX', i, q + j) for j, col in enumerate(cols)
                     for i in range(q) if col >> i & 1]
            frame += [('H', i) for i in range(q)]
            word = ([('T', i) for i in range(q)] + _adjoint(frame)
                    + [('T', q + j) for j in range(q)] + frame)
            native = _word_matrix(2 * q, word)
            packet = native @ probe @ native.conj().T
            expected = np.zeros_like(packet)
            marginal = np.zeros(1 << q)
            for u, v in product(range(1 << q), repeat=2):
                alpha = apply(transpose(cols), u)
                pauli = section(((1 << q) - 1) ^ apply(cols, v), u | (v << q))
                coefficient = (0 if v & ~alpha else
                               2 ** (-(q + alpha.bit_count()) / 2)
                               * 1j ** (u.bit_count() - v.bit_count()))
                expected += coefficient * pauli
                actual_coefficient = np.trace(pauli.conj().T @ packet) / size
                marginal[v] += abs(actual_coefficient) ** 2
            np.testing.assert_allclose(packet, expected, atol=2e-14, rtol=0)
            np.testing.assert_allclose(
                marginal,
                [0.25 ** v.bit_count() * 0.75 ** (q - v.bit_count())
                 for v in range(1 << q)], atol=2e-14, rtol=0)
            np.testing.assert_allclose(
                _word_matrix(2 * q, word + _adjoint(word)), np.eye(size),
                atol=2e-14, rtol=0)
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 2 * q)
            packets[cols] = packet

        relative_values = {}
        for left, right in product(matrices, repeat=2):
            relative = tuple(apply(transpose(right), col)
                             for col in transpose(inverse(left)))
            defect = tuple(col ^ (1 << j) for j, col in enumerate(relative))
            kernel_test = transpose(defect)
            expected_overlap = 0.0
            for alpha in range(1 << q):
                beta = apply(relative, alpha)
                common = sum(1 for v in range(1 << q)
                             if not v & ~(alpha & beta)
                             and apply(kernel_test, v) == 0)
                expected_overlap += common * 2 ** (
                    -q - (alpha.bit_count() + beta.bit_count()) / 2)
            overlap = np.trace(packets[left] @ packets[right]) / size
            self.assertAlmostEqual(overlap.real, expected_overlap, places=13)
            self.assertAlmostEqual(overlap.imag, 0, places=13)
            self.assertLessEqual(abs(overlap),
                                 0.75 ** _rank([a ^ b for a, b in zip(left, right)])
                                 + 2e-14)
            if defect in relative_values:
                self.assertAlmostEqual(overlap.real, relative_values[defect], places=13)
            relative_values[defect] = overlap.real
            if left != right:
                distance_squared = np.linalg.norm(
                    packets[left] - packets[right], 'fro') ** 2 / size
                self.assertGreaterEqual(distance_squared, 0.5 - 2e-14)

    def test_isotropic_layers_and_augmented_rank(self):
        rng = Random(773)
        for width in range(1, 8):
            for _ in range(50):
                layers = []
                for _ in range(rng.randint(1, 8)):
                    axes = [1 << i for i in range(rng.randint(1, width))]
                    for _ in range(4 * width):
                        pivot = rng.randrange(1, 1 << (2 * width))
                        axes = [a ^ pivot if _pair(a, pivot, width) else a for a in axes]
                    layers.append(axes)
                axes = [a for layer in layers for a in layer]
                t = len(axes)
                v = _rank(axes)
                twice_r = _gram_rank(axes, width)
                largest = sorted(map(len, layers), reverse=True)[:2]
                s = t - sum(largest)
                self.assertLessEqual(twice_r, 2 * (v - largest[0]))
                self.assertLessEqual(twice_r, s + max(0, 2 * v - t))
                boundary = [rng.randrange(1 << (2 * width)) for _ in range(6)]
                self.assertLessEqual(_gram_rank(axes + boundary, width),
                                     twice_r + 2 * len(boundary))

    def test_commuting_layer_encoding_repeated_signed_affine_axes(self):
        rng = Random(931)
        for width in range(1, 7):
            for _ in range(15):
                parities = [(rng.randrange(1 << width), rng.choice((-1, 1)), rng.randrange(2))
                            for _ in range(30)]
                coefs = [0] * (1 << width)
                for mask, sign, offset in parities:
                    coefs[0] += sign * offset
                    bits = [1 << j for j in range(width) if mask & (1 << j)]
                    for degree in range(1, 4):
                        for subset in combinations(bits, degree):
                            coefs[sum(subset)] += sign * (-1)**offset * (-2)**(degree - 1)
                for x in range(1 << width):
                    direct = sum(sign * (((x & mask).bit_count() % 2) ^ offset)
                                 for mask, sign, offset in parities)
                    model = sum(a for mask, a in enumerate(coefs) if x & mask == mask)
                    self.assertEqual(direct % 8, model % 8)

    def test_radial_profile_envelopes_at_all_breakpoints(self):
        for t in range(1, 22):
            for parts in _profiles(t, t):
                upper = min(Fraction(t, 4), Fraction(t, 2) - parts[0])
                if upper < 0:
                    continue
                s = t - sum(parts[:2])
                profile = Fraction(3 * t*t, 8) - sum(Fraction(k*k, 2) for k in parts)
                two = (Fraction(s*t, 2) + Fraction(s*s, 4) if 2*s <= t else
                       Fraction(3*t*t, 8) - Fraction((t-s)**2, 4))
                for r in {Fraction(0), upper, *(Fraction(k) for k in parts if k < upper)}:
                    penalty = sum((k*k - max(0, k-r)**2)/Fraction(2) for k in parts)
                    self.assertLessEqual(2*r*t - 2*r*r - penalty, min(profile, two))


if __name__ == '__main__':
    unittest.main()
