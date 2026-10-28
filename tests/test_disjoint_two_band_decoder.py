"""Exact finite support fingerprints for disjoint-axis two-band decoders.

These checks exercise the proof's collision-free expansion and robust gap.
They are not a search for all-size decoders or a bound on overlapping bands.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import Q2
from tests.test_one_band_copy_feasibility import _pauli_product


_ZERO, _ONE = Q2(0), Q2(1)
_ISQRT2 = Q2(0, Fraction(1, 2))


def _pair(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((right & mask) & (left >> width)).bit_count()) & 1


def _span(axes):
    result = {0}
    for axis in axes:
        result |= {v ^ axis for v in result}
    return result


def _rotate(row, axis, orientation, width):
    result = {}
    for label, value in row.items():
        images = ((label, value),)
        if _pair(label, axis, width):
            other, phase = _pauli_product(axis, label, width)
            phase = (phase - 1) % 4
            assert phase in (0, 2)
            sign = orientation * (1 if phase == 0 else -1)
            images = ((label, value * _ISQRT2),
                      (other, sign * value * _ISQRT2))
        for other, coefficient in images:
            total = result.get(other, _ZERO) + coefficient
            if total != _ZERO:
                result[other] = total
            else:
                result.pop(other, None)
    return result


def _image(label, inner, outer, signs, width):
    row = {label: _ONE}
    for axis, sign in zip(inner + outer, signs):
        row = _rotate(row, axis, sign, width)
    return row


def _normalized_bands(inner, outer, signs, width):
    first = list(zip(inner, signs[:len(inner)]))
    second = list(zip(outer, signs[len(inner):]))
    moved = [(a, s) for a, s in first
             if not any(_pair(a, b, width) for b in outer)]
    kept = [(a, s) for a, s in first if (a, s) not in moved]
    return (tuple(a for a, _ in kept),
            tuple(a for a, _ in second + moved),
            tuple(s for _, s in kept + second + moved))


class DisjointTwoBandDecoderTests(unittest.TestCase):
    def test_partial_band_normalization_and_outer_axis_recovery(self):
        width = 3
        x = [1 << j for j in range(width)]
        z = [1 << (width + j) for j in range(width)]
        templates = (
            ((), ()),
            ((), (z[0], z[2])),
            ((z[0], z[1]), (x[0], x[1])),
            ((z[0], z[1]), (x[0],)),
            ((z[0], z[1]), (x[0] ^ x[1], z[2])),
            ((z[0], z[1], z[2]),
             (x[0] ^ x[1], x[1] ^ x[2], x[0])),
            ((z[0], z[1], z[2]), (x[0] ^ x[1],)),
        )
        cases = 0
        for inner, outer in templates:
            self.assertEqual(_span(inner) & _span(outer), {0})
            for axes in (inner, outer):
                self.assertEqual(len(_span(axes)), 1 << len(axes))
                self.assertTrue(all(not _pair(a, b, width)
                                    for a in axes for b in axes))
            for pattern in range(4):
                signs = tuple(-1 if (j + pattern) % 3 == 0 else 1
                              for j in range(len(inner) + len(outer)))
                first, second, reordered = _normalized_bands(
                    inner, outer, signs, width)
                recovered = set()
                for label in range(1 << (2 * width)):
                    image = _image(label, inner, outer, signs, width)
                    self.assertEqual(image, _image(label, first, second,
                                                   reordered, width))
                    active = [a for a in first if _pair(label, a, width)]
                    expected_squares = Counter()
                    for branch in _span(active):
                        weight = sum(_pair(label ^ branch, b, width)
                                     for b in second)
                        expected_squares[Fraction(1, 1 << (len(active) + weight))] += 1 << weight
                    actual_squares = Counter()
                    for coefficient in image.values():
                        square = coefficient * coefficient
                        self.assertEqual(square.b, 0)
                        actual_squares[square.a] += 1
                    self.assertEqual(actual_squares, expected_squares)
                    self.assertEqual(sum((c * c for c in image.values()), _ZERO), _ONE)
                    is_two = not active and sum(_pair(label, b, width)
                                                for b in second) == 1
                    self.assertEqual(len(image) == 2, is_two)
                    if is_two:
                        left, right = image
                        recovered.add(left ^ right)
                        self.assertEqual({abs(c) for c in image.values()}, {_ISQRT2})
                    else:
                        largest = sorted((abs(c) for c in image.values()),
                                         reverse=True)[:2]
                        overlap = sum(largest, _ZERO) * _ISQRT2
                        # Exact squared distance to the nearest signed flat two-term sum.
                        self.assertGreaterEqual(2 - 2 * overlap, Q2(1, Fraction(-1, 2)))
                    cases += 1
                self.assertEqual(recovered, set(second))
        self.assertEqual(cases, 1792)

    def test_gap_is_attained_and_two_term_supports_are_separated(self):
        # E=T(X) T(Z), P=X has magnitudes 1/sqrt(2),1/2,1/2.
        row = _image(1, (2,), (1,), (1, 1), 1)
        self.assertEqual(sorted(abs(c) for c in row.values()),
                         [Q2(Fraction(1, 2)), Q2(Fraction(1, 2)), _ISQRT2])
        top = sorted((abs(c) for c in row.values()), reverse=True)[:2]
        self.assertEqual(2 - 2 * sum(top, _ZERO) * _ISQRT2,
                         Q2(1, Fraction(-1, 2)))
        flats = []
        for left in range(4):
            for right in range(left + 1, 4):
                for signs in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                    flats.append({left: signs[0] * _ISQRT2,
                                  right: signs[1] * _ISQRT2})
        for i, left in enumerate(flats):
            for right in flats[:i]:
                distance2 = _ZERO
                for j in range(4):
                    difference = left.get(j, _ZERO) - right.get(j, _ZERO)
                    distance2 += difference * difference
                self.assertGreaterEqual(distance2, 1)

    def test_enlarged_field_and_moving_residual_keep_distinct_deviations(self):
        # Source: p=2 over GF(16), eight core wires. Target: four bank wires
        # and a genuine p=2 GF(4) shear on four moving residual wires.
        def multiply_field(left, right, degree, modulus):
            value = 0
            while right:
                if right & 1:
                    value ^= left
                right >>= 1
                left <<= 1
                if left & (1 << degree):
                    left ^= modulus
            return value

        def shear(label, address, degree, modulus):
            width = 2 * degree
            mask = (1 << degree) - 1
            a_u, a_v = label & mask, (label >> degree) & mask
            b_u, b_v = (label >> width) & mask, label >> (width + degree)
            a_v ^= multiply_field(address, a_u, degree, modulus)
            for j in range(degree):
                if (multiply_field(address, 1 << j, degree, modulus) & b_v).bit_count() & 1:
                    b_u ^= 1 << j
            return a_u | (a_v << degree) | (b_u << width) | (b_v << (width + degree))

        source = lambda v, x: shear(v, x, 4, 0b10011)

        def target(label, address):
            x, z = label & 255, label >> 8
            bank_x = sum(((x >> j) & 1) << (j ^ address) for j in range(4))
            bank_z = sum(((z >> j) & 1) << (j ^ address) for j in range(4))
            residual = shear((x >> 4) | ((z >> 4) << 4), address, 2, 0b111)
            return (bank_x | ((residual & 15) << 4)
                    | (bank_z << 8) | ((residual >> 4) << 12))

        for left in range(4):
            for right in range(4):
                for j in range(16):
                    first = source(1 << j, right) ^ (1 << j)
                    self.assertEqual(source(first, left) ^ first, 0)
        bank = 1  # X on bank wire zero, with no residual component.
        first = target(bank, 2) ^ bank
        second = target(first, 1) ^ first
        self.assertEqual(second, 15)  # X on every bank wire, never canceled.
        self.assertNotEqual(target(1 << 4, 1), 1 << 4)

    def test_shared_spans_can_have_actual_path_collisions(self):
        # Four mutually commuting axes with one relation; this lies outside the theorem.
        width = 3
        inner, outer = (8, 16), (12, 20)
        self.assertEqual(_span(inner) & _span(outer), {0, 24})
        label = 3
        row = _image(label, inner, outer, (1, 1, 1, 1), width)
        paths = sum(1 << sum(_pair(label ^ branch, b, width) for b in outer)
                    for branch in _span(inner))
        self.assertLess(len(row), paths)


if __name__ == '__main__':
    unittest.main()
