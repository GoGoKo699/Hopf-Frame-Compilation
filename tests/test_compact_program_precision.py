
"""Finite checks for the compact-program field reference and flag schedule."""
from itertools import product
from math import acos, pi
from random import Random
import unittest

import numpy as np


def _multiply(a, b, polynomial):
    degree = polynomial.bit_length() - 1
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a >> degree:
            a ^= polynomial
    return result


def _field_order(polynomial):
    size = 1 << (polynomial.bit_length() - 1)
    for generator in range(1, size):
        orbit = [1]
        for _ in range(size - 2):
            orbit.append(_multiply(orbit[-1], generator, polynomial))
        if len(set(orbit)) == size - 1:
            return orbit
    raise AssertionError("No primitive generator in the fixture field")


def _reference(orbit):
    size = len(orbit) + 1
    p = 1 / size
    phi = acos(1 - 1 / (2 * (1 - p)))
    marked_phase = np.exp(1j * phi)
    overlap = 1 - p + p * marked_phase
    reflection_phase = 1 - marked_phase / overlap
    uniform = np.ones(size, dtype=complex) / np.sqrt(size)
    index_state = uniform.copy()
    index_state[-1] *= marked_phase
    index_state += (reflection_phase - 1) * uniform * np.vdot(uniform, index_state)
    index_state *= np.exp(-1j * np.angle(index_state[0]))

    # Z rotations on the binary digits, with the literal common minus.
    denominator = size - 1
    digits = np.arange(size)
    for j in range(size.bit_length() - 1):
        sign = 1 - 2 * ((digits >> j) & 1)
        index_state *= np.exp(1j * pi * (1 << j) * sign / denominator)
    reference = np.zeros(size, dtype=complex)
    reference[orbit + [0]] = -index_state
    return reference


def _field_shears(a, b, beta, inverse, polynomial, active, reverse=False):
    if not active:
        return a, b
    if reverse:
        a, b = b, a
        a ^= _multiply(inverse, b, polynomial)
        b ^= _multiply(beta, a, polynomial)
        return a, b
    b ^= _multiply(beta, a, polynomial)
    a ^= _multiply(inverse, b, polynomial)
    return b, a


def _dense_frame(angles, height):
    size = 1 << height
    frame = np.eye(size)
    for depth, layer in enumerate(angles):
        remaining = height - 1 - depth
        target = 1 << remaining
        for prefix, theta in enumerate(layer):
            low = prefix << (remaining + 1)
            columns = [low, low | target]
            c, s = np.cos(theta), np.sin(theta)
            frame[columns, :] = np.array([[c, -s], [s, c]]) @ frame[columns, :]
    return frame


def _flagged_column(original, angles, height, enabled):
    """Retain stored initial suffix bits and written-prefix flags on each path."""
    guards = tuple(
        (original & ((1 << (height - 1 - depth)) - 1)) == 0
        for depth in range(height)
    )
    paths = {(original, guards, ()): 1.0}
    for depth, layer in enumerate(angles):
        remaining = height - 1 - depth
        target = 1 << remaining
        updated = {}
        for (address, saved, prefixes), amplitude in paths.items():
            prefix = address >> (remaining + 1)
            if prefixes and prefixes[-1] != prefix:
                raise AssertionError("A written prefix changed after its layer")
            following_guard = True if not remaining else saved[depth + 1]
            next_bit = 0 if not remaining else address & (target >> 1)
            erase = following_guard and not next_bit
            returned = saved[:depth] + (saved[depth] ^ erase,) + saved[depth + 1:]
            theta = layer[prefix]
            if enabled and saved[depth]:
                changed = -np.sin(theta) if address & target else np.sin(theta)
                terms = ((address, np.cos(theta)), (address ^ target, changed))
            else:
                terms = ((address, 1.0),)
            for output, coefficient in terms:
                child = prefixes + (output >> remaining,)
                key = output, returned, child
                updated[key] = updated.get(key, 0.0) + amplitude * coefficient
        paths = updated

    column = np.zeros(1 << height)
    for (output, saved, prefixes), amplitude in paths.items():
        expected = tuple(output >> (height - j) for j in range(1, height + 1))
        if any(saved) or expected != prefixes:
            raise AssertionError("A suffix or prefix flag fails to uncompute")
        column[output] += amplitude
    return column


class CompactProgramPrecisionTests(unittest.TestCase):
    def test_literal_field_reference_and_eigenphases(self):
        for polynomial in (0b111, 0b1011, 0b10011):
            orbit = _field_order(polynomial)
            reference = _reference(orbit)
            length = len(orbit)
            desired = np.zeros(length + 1, dtype=complex)
            desired[orbit] = np.exp(-2j * pi * np.arange(length) / length) / np.sqrt(length)
            np.testing.assert_allclose(reference, desired, atol=2e-14)
            for d, beta in enumerate(orbit):
                with self.subTest(polynomial=polynomial, exponent=d):
                    inverse = orbit[-d % length]
                    image = [_field_shears(a, 0, beta, inverse, polynomial, True)[0]
                             for a in range(length + 1)]
                    shifted = np.zeros_like(reference)
                    shifted[image] = reference
                    np.testing.assert_allclose(
                        shifted, np.exp(2j * pi * d / length) * reference, atol=2e-14
                    )

    def test_field_multiplier_returns_clean_temp_and_reverses_occupied_temp(self):
        polynomial = 0b1011
        orbit = _field_order(polynomial)
        for d, beta in enumerate(orbit):
            inverse = orbit[-d % len(orbit)]
            for a, b, enabled in product(range(8), range(8), (False, True)):
                out = _field_shears(a, b, beta, inverse, polynomial, enabled)
                if not enabled:
                    self.assertEqual(out, (a, b))
                elif b == 0:
                    self.assertEqual(out, (_multiply(beta, a, polynomial), 0))
                self.assertEqual(
                    _field_shears(*out, beta, inverse, polynomial, enabled, reverse=True),
                    (a, b),
                )

    def test_prefix_and_suffix_flags_return_on_every_logical_column(self):
        rng = np.random.default_rng(322)
        for height in range(1, 6):
            angles = tuple(rng.uniform(-pi, pi, size=1 << d) for d in range(height))
            frame = _dense_frame(angles, height)
            for enabled in (False, True):
                with self.subTest(height=height, enabled=enabled):
                    actual = np.column_stack([
                        _flagged_column(x, angles, height, enabled)
                        for x in range(1 << height)
                    ])
                    expected = frame if enabled else np.eye(1 << height)
                    np.testing.assert_allclose(actual, expected, atol=2e-14)


    def test_circulant_threshold_corners_cancel_dirty_counter_and_words(self):
        rng = Random(410)
        for count in range(1, 5):
            span = 1 << (2 * count + 1).bit_length()
            for threshold in range(count + 2):
                cutoff = threshold + count + 1
                rows = tuple(
                    sum(1 << j for j in range(span) if (j - i) % span >= cutoff)
                    for i in range(span)
                )

                def corner(left, right):
                    return sum((row & right).bit_count() for i, row in enumerate(rows)
                               if left >> i & 1) % 2

                words = (list(product(range(1 << span), repeat=2)) if count == 1
                         else [(rng.getrandbits(span), rng.getrandbits(span))
                               for _ in range(12)])
                for total, counter, character in product(
                        range(count + 1), range(span), (0, 1)):
                    moved = (counter + total + (count + 1) * character) % span
                    expected = character and total >= threshold
                    for left, right in words:
                        actual = corner(left, right)
                        next_left = left ^ (1 << counter)
                        actual ^= corner(next_left, right)
                        next_right = right ^ (1 << moved)
                        actual ^= corner(next_left, next_right)
                        next_left ^= 1 << counter
                        actual ^= corner(next_left, next_right)
                        next_right ^= 1 << moved
                        self.assertEqual((next_left, next_right), (left, right))
                        self.assertEqual(actual, expected)



if __name__ == "__main__":
    unittest.main()
