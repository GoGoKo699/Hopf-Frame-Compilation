"""Finite layer-count bounds and correlated native Pauli-packet checks."""

from fractions import Fraction
from itertools import combinations, product
from random import Random
import unittest

import numpy as np

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, native_columns, scale,
)
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

    def test_actual_band_insertion_response_and_return(self):
        # The odd-size companion of t^3+t+1 and its difference from I
        # are invertible. Emit their replacement in one common native
        # frame, independently of the two B words it must connect.
        q = 3
        width = 2 * q
        size = 1 << width
        old = (1, 2, 4)
        new = (2, 4, 3)
        self.assertEqual(_rank(new), q)
        self.assertEqual(_rank([a ^ b for a, b in zip(old, new)]), q)

        def apply(cols, word):
            value = 0
            for j, col in enumerate(cols):
                if word >> j & 1:
                    value ^= col
            return value

        def transpose(cols):
            n = len(cols)
            return tuple(sum(((cols[j] >> i) & 1) << j for j in range(n))
                         for i in range(n))

        def inverse(cols):
            return tuple(next(x for x in range(1 << len(cols))
                              if apply(cols, x) == 1 << j)
                         for j in range(len(cols)))

        def b_word(cols):
            frame = [('CX', i, q + j) for j, col in enumerate(cols)
                     for i in range(q) if col >> i & 1]
            frame += [('H', i) for i in range(q)]
            return (_adjoint(frame) + [('T', q + j) for j in range(q)]
                    + frame)

        def trace(matrix):
            return np.trace(matrix) / size

        # C=H_first U_{D^{-T}} maps the 2q distinct Z axes to the
        # combined new/old lists. Gaussian elimination only synthesizes
        # this Clifford frame; it does not use either response formula.
        axis_cols = tuple(col | (1 << (q + j))
                          for cols in (new, old) for j, col in enumerate(cols))
        self.assertEqual(_rank(axis_cols), width)
        rows = list(inverse(axis_cols))  # rows of D^{-T}
        elimination = []

        def add_row(control, target):
            rows[target] ^= rows[control]
            elimination.append(('CX', control, target))

        for pivot in range(width):
            source = next(i for i in range(pivot, width) if rows[i] >> pivot & 1)
            if source != pivot:
                add_row(source, pivot)
                add_row(pivot, source)
                add_row(source, pivot)
            for target in range(width):
                if target != pivot and rows[target] >> pivot & 1:
                    add_row(pivot, target)
        self.assertEqual(rows, [1 << j for j in range(width)])
        frame = list(reversed(elimination)) + [('H', i) for i in range(q)]
        layer = [('T', i) for i in range(q)] + [('TDG', q + i) for i in range(q)]
        bridge_word = _adjoint(frame) + layer + frame
        self.assertEqual(len({gate[1] for gate in layer}), width)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in bridge_word), 2 * q)

        a_word = [('T', i) for i in range(q)]
        a = _word_matrix(width, a_word)
        b_old, b_new = (_word_matrix(width, b_word(cols)) for cols in (old, new))
        bridge = _word_matrix(width, bridge_word)
        np.testing.assert_allclose(bridge, b_new @ b_old.conj().T, atol=4e-14, rtol=0)
        np.testing.assert_allclose(bridge @ b_old @ a, b_new @ a, atol=4e-14, rtol=0)
        return_word = (a_word + b_word(old) + bridge_word
                       + _adjoint(b_word(new)) + _adjoint(a_word))
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in return_word), 6 * q)
        # Exact Q(sqrt(2), i), every physical input column, including
        # arbitrary dirty superpositions: no prefix return is assumed.
        self.assertEqual(native_columns(width, return_word), eye(size))

        probe = _word_matrix(width, [('X', i) for i in range(q)])
        first_packet = a @ probe @ a.conj().T
        packets = [b @ first_packet @ b.conj().T for b in (b_old, b_new)]
        self.assertLessEqual(abs(trace(packets[1] @ packets[0])), 0.75 ** q + 4e-14)
        self.assertAlmostEqual(trace(packets[1] @ bridge @ packets[0]
                                     @ bridge.conj().T), 1, places=13)

        gamma = 2 ** -0.5
        responses = []
        for cols, packet in zip((old, new), packets):
            values = {}
            for w in range(1, 1 << q):
                pivot = (w & -w).bit_length() - 1
                parity = [('CX', i, pivot) for i in range(q)
                          if i != pivot and w >> i & 1]
                insertion = _word_matrix(
                    width, parity + [('T', pivot)] + _adjoint(parity))
                actual = trace(packet @ insertion @ packet @ insertion.conj().T)
                expected = ((1 + gamma) / 2 + (1 - gamma) / 2
                            * (-1) ** w.bit_count()
                            * 2 ** -apply(transpose(cols), w).bit_count())
                self.assertAlmostEqual(actual, expected, places=13)
                values[w] = actual
            responses.append(values)
        self.assertAlmostEqual(abs(responses[0][2] - responses[1][2]),
                               (1 - gamma) / 8, places=13)

        self.assertAlmostEqual(trace(first_packet @ packets[0]),
                               ((1 + gamma) / 2) ** q, places=13)
        np.testing.assert_allclose(b_old.conj().T @ packets[0] @ b_old,
                                   first_packet, atol=4e-14, rtol=0)

        # A nonzero parity polynomial with E=B_old keeps BOTH exterior
        # B words. Compare its scalar phase average with native gates.
        phase_average = 0.0
        omega = np.exp(1j * np.pi / 4)
        for z in range(1 << q):
            mean = sum(omega ** (2 * (apply(transpose(old), x) ^ z).bit_count()
                                - (apply(transpose(new), x) ^ z).bit_count())
                       for x in range(1 << q)) / (1 << q)
            phase_average += abs(mean) ** 2 / (1 << q)
        actual = trace(packets[1] @ b_old @ packets[0] @ b_old.conj().T)
        self.assertAlmostEqual(actual, phase_average, places=13)

    def test_terminal_quotient_initialized_profiles_and_clean_support(self):
        # Three physical wires: dirty 0, logical 1, clean 2. Exact native
        # columns are compared with a quotient-phase decoder on all eight
        # members of an affine Pauli family, including its phase carry.
        width, size = 3, 8
        i_power = (CO, CI, -CO, -CI)
        root_half = Q2(0, Fraction(1, 2))
        omega = CQ2(root_half, root_half)
        eighth_roots = [CO]
        for _ in range(7):
            eighth_roots.append(eighth_roots[-1] * omega)

        def pauli(w, x, z, phase=CO):
            result = [[CZ for _ in range(1 << w)] for _ in range(1 << w)]
            for column in range(1 << w):
                result[column ^ x][column] = phase * (-1) ** (z & column).bit_count()
            return result

        def conjugate(unitary, operator):
            return multiply(multiply(unitary, operator), adjoint(unitary))

        def trace(matrix):
            return sum((row[j] for j, row in enumerate(matrix)), CZ)

        observable = pauli(width, 3, 1, CI)  # Y_dirty X_logical
        quotient_phase = [[eighth_roots[(row >> 1) & 1] if row == col else CZ
                           for col in range(size)] for row in range(size)]
        swap = [('CX', 0, 2), ('CX', 2, 0), ('CX', 0, 2)]

        for full_support in (True, False):
            with self.subTest(full_support=full_support):
                if full_support:
                    q = 2
                    c_word = swap + [('H', 0), ('S', 1)]
                    # Prefix T_2 changes the actual dirty input. Its later
                    # inverse cancels; the fixed section y_2=0 saw neither.
                    middle = [('T', 2), ('T', 1), ('TDG', 2)]
                    logical = [('S', 1), ('T', 1), ('H', 1), ('T', 1),
                               ('H', 1), ('TDG', 1), ('SDG', 1)]
                    dirty_image = 2
                else:
                    q = 1
                    c_word = [('S', 1)]
                    # p(y)=y_1+2*y_0*y_2. The clean slice is y_2=0,
                    # so p'=y_1 agrees on the code and differs off it.
                    middle = [('T', 1), ('T', 0), ('T', 2),
                              ('CX', 2, 0), ('TDG', 0), ('CX', 2, 0)]
                    logical = [('T', 1)]
                    dirty_image = 0
                outer = [('T', j) for j in range(q)] + [('H', j) for j in range(q)]
                b_word = c_word + outer
                word = b_word + middle + _adjoint(b_word)
                native = native_columns(width, word)
                expected = native_columns(width, logical)
                code = [row[:4] for row in native]
                self.assertEqual(code, [row[:4] for row in expected])
                if full_support:
                    self.assertEqual(native, expected)
                    self.assertEqual(
                        native_columns(width, b_word + [('T', 2)] + _adjoint(b_word)),
                        native_columns(width, [('T', 0)]))
                    for t in range(4):
                        y = t  # fixed quotient section y_2=0
                        restricted_bands = [(y >> 2) & 1, (y >> 1) & 1,
                                            -((y >> 2) & 1)]
                        self.assertEqual(restricted_bands, [0, (t >> 1) & 1, 0])
                else:
                    self.assertNotEqual(native, expected)
                    self.assertNotEqual(quotient_phase,
                                        native_columns(width, middle))

                c_matrix = native_columns(width, c_word)
                outer_matrix = native_columns(width, outer)
                b_code = native_columns(width, b_word, columns=range(4))
                projector = multiply(b_code, adjoint(b_code))
                clean_image = pauli(width, 1, 0) if full_support else pauli(width, 0, 4)
                self.assertEqual(projector, scale(add(eye(size), conjugate(
                    outer_matrix, clean_image)), Fraction(1, 2)))
                if full_support:
                    self.assertEqual([projector[y][y] for y in range(size)],
                                     [CQ2(Q2(Fraction(1, 2),
                                              Fraction(-1 if y & 1 else 1, 4)))
                                      for y in range(size)])
                    self.assertNotEqual(projector[0][1], CZ)
                else:
                    self.assertEqual(projector,
                                     [[CO if row == col and row < 4 else CZ
                                       for col in range(size)] for row in range(size)])
                reconstructed_code = multiply(
                    adjoint(native_columns(width, b_word)),
                    multiply(quotient_phase, b_code))
                self.assertEqual(reconstructed_code, code)

                # Store the original phase bc and Clifford cocycle a
                # separately. Their sum has cubic high bit abc, so assuming
                # quadratic phase digits after conjugation would be wrong.
                recorded_observable = pauli(
                    width, (1 << dirty_image) | 2, (1 << dirty_image) | 2, -CO)
                self.assertEqual(conjugate(c_matrix, observable), recorded_observable)
                decoded_observable = conjugate(outer_matrix, recorded_observable)
                decoded_middle = multiply(multiply(adjoint(quotient_phase),
                                                    decoded_observable), quotient_phase)
                direct_observable = multiply(multiply(adjoint(code), observable), code)
                high_bits = []
                for a, b, c in product((0, 1), repeat=3):
                    original = pauli(width, b | (a << 1), c, i_power[b * c])
                    recorded = pauli(width, (b << dirty_image) | (a << 1),
                                     (c << dirty_image) | (a << 1),
                                     i_power[(b * c + a) % 4])
                    self.assertEqual(conjugate(c_matrix, original), recorded)
                    high_bits.append((b * c + a) // 2)
                    probe_code = [row[:4] for row in original[:4]]
                    direct = trace(multiply(probe_code, direct_observable)) / 4
                    transported = conjugate(outer_matrix, recorded)
                    decoded = trace(multiply(multiply(projector, transported),
                                             decoded_middle)) / 4
                    self.assertEqual(direct, decoded)
                    if a == b == c == 1:
                        expected_bias = (Q2(1) + root_half) / 2 if full_support else root_half
                        self.assertEqual(direct, CQ2(expected_bias))
                self.assertEqual(sum(high_bits) % 2, 1)

        # A single independent-axis native band with deficient crossing
        # invalidates the clean-projector support dichotomy outside the
        # stated transverse grammar: T(X_0)T(X_0 Z_1).
        frame = [('CX', 0, 1), ('H', 0)]
        band = _adjoint(frame) + [('T', 0), ('T', 1)] + frame
        columns = native_columns(2, band, columns=(0, 2))
        projector = multiply(columns, adjoint(columns))
        expected = [[CZ for _ in range(4)] for _ in range(4)]
        expected[0][0] = expected[1][1] = CQ2(Fraction(1, 2))
        expected[0][1], expected[1][0] = CI / 2, -CI / 2
        expected[2][2] = CO
        self.assertEqual(projector, expected)

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
