"""Exact native-coordinate and signed-defect checks for middle orbit removal.

These are algebra/native-defect fixtures.  They do not supply the endpoint
encoders, one-clean return, or globally Boolean occupied phase of a complete
query source.  In particular, a commuting Gram calculation alone is not the
literal cubic identity used by the growing-displacement argument.
"""
from itertools import product
import unittest

from tests.exact_matrices import CO, Gate, apply_word, native_columns
from tests.test_bounded_relation_native_matching import _rank
from tests.test_external_center_middle_correction import _linear_word
from tests.test_two_layer_galois_return import _inverse


WIDTH = 8
P = (1 << 6) - 1
R = (1 << 6) | (1 << 7)
V = 1
BASIS = tuple(1 << bit for bit in range(WIDTH))


def _dot(x, y):
    return (x & y).bit_count() & 1


def _beta(x, y):
    # Polar form of e2(x_0,...,x_5) + (p.x) x_0 + x_0.
    return (_dot(P, x) * _dot(P, y) ^ _dot(P & x, y)
            ^ _dot(P, x) * _dot(V, y) ^ _dot(V, x) * _dot(P, y))


def _displacement(probe, x):
    return ((P if probe & 1 and _dot(P, x) else 0)
            ^ (R if probe & 2 and _dot(R, x) else 0))


def _cubic(action, x, y, z):
    return ((action(x) & action(y) & action(z)).bit_count()
            ^ (x & y & z).bit_count()) & 1


def _wedge(flag, x, y, z):
    return (_beta(x, y) * _dot(flag, z)
            ^ _beta(x, z) * _dot(flag, y)
            ^ _beta(y, z) * _dot(flag, x))


def _native(word):
    return [Gate(name, tuple(qubits)) for name, *qubits in word]


def _roots():
    omega = native_columns(1, [('T', 0)])[1][1]
    result = [CO]
    for _ in range(7):
        result.append(result[-1] * omega)
    return result


class GrowingDisplacementNativeOrbitTests(unittest.TestCase):
    def test_simultaneous_cubic_identity_product_algebra_and_native_orbits(self):
        beta_columns = [sum(_beta(x, y) << j for j, y in enumerate(BASIS))
                        for x in BASIS]
        self.assertEqual(_rank(beta_columns), 4)
        self.assertEqual(_rank(beta_columns + [P, R]), 6)
        for direction, x in product((P, R), BASIS):
            self.assertEqual(_beta(direction, x), 0)
        # Both displacement generators are actual native-coordinate actions:
        # one six-axis transvection and one independent transposition.
        for probe in range(4):
            flag = P if probe & 1 else 0
            action = lambda x, u=probe: x ^ _displacement(u, x)
            self.assertEqual(_displacement(probe, V), flag)
            for x, y in product(BASIS, repeat=2):
                self.assertEqual(_dot(_displacement(probe, x), y),
                                 _dot(x, _displacement(probe, y)))
                self.assertEqual(_dot(action(x), action(y)), _dot(x, y))
            for other, x in product(range(4), BASIS):
                self.assertEqual(
                    _displacement(probe, _displacement(other, x)), 0)
            # Repeated as well as distinct indices: the full trilinear tensor.
            for x, y, z in product(BASIS, repeat=3):
                self.assertEqual(_cubic(action, x, y, z),
                                 _wedge(flag, x, y, z))

        kernel = {x for x in range(1 << WIDTH)
                  if not _dot(P, x) and not _dot(R, x)}
        self.assertEqual(len(kernel), 64)
        x, y = 1 | 2, 1 | 4
        self.assertIn(x, kernel)
        self.assertIn(y, kernel)
        self.assertEqual(_beta(x, y), 1)
        self.assertEqual(x & y, V)
        for x, y in product(kernel, repeat=2):
            self.assertEqual(_beta(x, y), _dot(P & x, y))
            self.assertIn((x & y) ^ (V if _beta(x, y) else 0), kernel)
        algebra = kernel | {x ^ V for x in kernel}
        self.assertEqual(len(algebra), 128)
        self.assertEqual(algebra,
                         {x for x in range(1 << WIDTH) if not _dot(R, x)})
        for x, y in product(algebra, repeat=2):
            self.assertIn(x & y, algebra)

        # Only the final two original native axes form a permutation orbit.
        permuted = [axis for axis in BASIS if all(
            (axis ^ _displacement(probe, axis)).bit_count() == 1
            for probe in range(4))]
        self.assertEqual(permuted, [1 << 6, 1 << 7])
        self.assertEqual(_rank([_displacement(u, x)
                                for u, x in product(range(4), BASIS)]), 2)
        self.assertEqual(_rank([_displacement(u, x) & P
                                for u, x in product(range(4), BASIS)]), 1)

    def test_signed_native_pair_defect_and_actual_reversed_word(self):
        roots = _roots()
        rows = [axis ^ _displacement(3, axis) for axis in BASIS]
        change = _linear_word(rows)
        # Eight independent commuting signed native axes.  X6 T6 X6 is
        # T_{-Z6}, so its scalar omega is retained in every phase below.
        band = ([('T', bit) for bit in range(6)]
                + [('X', 6), ('T', 6), ('X', 6), ('T', 7)])
        defect = change + _inverse(band) + _inverse(change) + band
        inverse = _inverse(defect)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in band), 8)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in defect), 16)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in inverse), 16)
        native_band, native_defect, native_inverse = map(
            _native, (band, defect, inverse))
        for state in range(1 << WIDTH):
            bits = [(state >> j) & 1 for j in range(WIDTH)]
            parity = _dot(P, state)
            quadratic = ((sum(bits[i] * bits[j]
                              for i in range(6) for j in range(i + 1, 6))) & 1)
            control = 1 ^ quadratic ^ parity * bits[0] ^ bits[0]
            exponent = (sum(bits[:6]) + 1 - bits[6] + bits[7]) % 8
            self.assertEqual(apply_word(native_band, {state: CO}),
                             {state: roots[exponent]})
            # The swap contributes the literal Clifford S6^dagger S7;
            # its signed native axes cannot be silently replaced by positives.
            phase = (4 * parity * control - 2 * bits[6] + 2 * bits[7]) % 8
            self.assertEqual(apply_word(native_defect, {state: CO}),
                             {state: roots[phase]})
            self.assertEqual(apply_word(native_inverse, {state: CO}),
                             {state: roots[-phase % 8]})

        # Every pair of axis orientations and T/T-dagger directions yields
        # the stated Clifford defect, including the band scalar constants.
        swap = [('CX', 0, 1), ('CX', 1, 0), ('CX', 0, 1)]
        for orientations in product(range(2), repeat=2):
            for directions in product((-1, 1), repeat=2):
                pair = []
                for bit, (orientation, direction) in enumerate(
                        zip(orientations, directions)):
                    if orientation:
                        pair.append(('X', bit))
                    pair.append(('T' if direction == 1 else 'TDG', bit))
                    if orientation:
                        pair.append(('X', bit))
                word = swap + _inverse(pair) + _inverse(swap) + pair
                coefficients = [direction * (1 - 2 * orientation)
                                for orientation, direction in
                                zip(orientations, directions)]
                for state in range(4):
                    x0, x1 = state & 1, state >> 1
                    band_phase = sum(direction * (((state >> bit) & 1) ^ orientation)
                                     for bit, (orientation, direction) in enumerate(
                                         zip(orientations, directions))) % 8
                    phase = ((coefficients[0] - coefficients[1]) * (x0 - x1)) % 8
                    self.assertEqual(phase % 2, 0)
                    self.assertEqual(apply_word(_native(pair), {state: CO}),
                                     {state: roots[band_phase]})
                    self.assertEqual(apply_word(_native(word), {state: CO}),
                                     {state: roots[phase]})
                    self.assertEqual(apply_word(_native(_inverse(word)), {state: CO}),
                                     {state: roots[-phase % 8]})

    def test_symmetric_square_zero_crosscoupling_fails_literal_cubic(self):
        def displacement(x):
            return (P if _dot(R, x) else 0) ^ (R if _dot(P, x) else 0)

        action = lambda x: x ^ displacement(x)
        self.assertEqual(_rank([displacement(x) for x in BASIS]), 2)
        self.assertEqual(displacement(V), R)
        for x, y in product(BASIS, repeat=2):
            self.assertEqual(_dot(displacement(x), y), _dot(x, displacement(y)))
            self.assertEqual(displacement(displacement(x)), 0)
            self.assertEqual(_dot(action(x), action(y)), _dot(x, y))
        # It has the same U, K, beta-support separation, and common candidate
        # v as above.  Those binary conditions do not imply native matching.
        failures = [(x, y, z) for x, y, z in product(BASIS, repeat=3)
                    if _cubic(action, x, y, z) != _wedge(R, x, y, z)]
        witness = (1, 2, 1 << 6)
        self.assertIn(witness, failures)
        self.assertEqual(_cubic(action, *witness), 1)
        self.assertEqual(_wedge(R, *witness), 0)


if __name__ == '__main__':
    unittest.main()
