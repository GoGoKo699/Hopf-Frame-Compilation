"""Bounded signed checks for the native commutation-cover certificate.

The source theorem is analytic and arbitrary-width. This axis fixture has
four composite bad zero-character axes and a good middle axis outside the
first span; it is not an optimized native synthesis or a returned query.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

from tests.test_native_band_source_structure import (
    _band, _is_c3, _pair, _rank,
)
from tests.test_outer_axis_normalizer import _image


class NativeCommutationCover(unittest.TestCase):
    def test_composite_zero_axes_have_literal_complement_symmetric_phase(self):
        # A consists of the four Z axes; each B axis is the product of the
        # other three. The two lists are independent and have no duplicates.
        values = []
        for word in range(16):
            exponent = word.bit_count() + sum(
                (word & (15 ^ (1 << i))).bit_count() & 1 for i in range(4))
            self.assertEqual(exponent % 8, 0 if word in (0, 15) else 4)
            values.append(exponent % 8)
        # Complement symmetry makes the complete phase commute with XXXX.
        self.assertTrue(all(values[word] == values[word ^ 15]
                            for word in range(16)))

    def test_nonduplicate_source_signed_stripping_and_actual_inverse(self):
        width = 5
        a = [1 << (width + i) for i in range(width)]
        b = [((15 ^ (1 << i)) << width) for i in range(4)] + [16]
        c = [15, 16]
        for axes in (a, b, c):
            self.assertEqual(_rank(axes), len(axes))
            self.assertTrue(all(_pair(p, q, width) == 0
                                for p, q in itertools.combinations(axes, 2)))
        self.assertEqual(_rank(a) + _rank(b) - _rank(a + b), 4)
        b_zero, b_good = b[:4], b[4:]
        a_zero, a_strip = a[:4], a[4:]
        self.assertTrue(all(p not in a for p in b_zero))
        self.assertTrue(all(_pair(p, q, width) == 0
                            for p in b_zero for q in a))
        self.assertTrue(all(any(_pair(p, q, width) for q in c)
                            for p in b_zero))
        self.assertTrue(all(_pair(p, q, width) == 0
                            for p in b_good for q in c))
        self.assertEqual(_rank(a + b_good), len(a) + 1)
        self.assertEqual(_pair(a_strip[0], b_good[0], width), 1)

        # Both native orientations exercise the signed derivative. On the
        # fifth wire the source is S_X^{orientation} T_Z^{orientation}.
        for orientation in (1, -1):
            with self.subTest(orientation=orientation):
                a_band = _band(a, width, [orientation] * len(a))
                b_band = _band(b, width, [orientation] * len(b))
                c_band = _band(c, width, [orientation] * len(c))
                strip_band = _band(a_strip, width, [orientation])
                source = c_band @ b_band @ a_band
                stripped = source @ strip_band.conj().T
                diagonal_inner = b_band @ _band(
                    a_zero, width, [orientation] * len(a_zero))

                self.assertTrue(_is_c3(source, width))
                self.assertTrue(_is_c3(stripped, width))
                self.assertTrue(_is_c3(diagonal_inner, width))
                self.assertEqual(_image(source, a_strip[0], width),
                                 (a_strip[0] ^ b_good[0],
                                  int(orientation == 1)))
                np.testing.assert_allclose(
                    stripped, c_band @ diagonal_inner, atol=2e-14)

                # Reverse the actual phase factors and their orientations.
                # Matrix equality checks every column and literal scalar.
                inverse = (_band(list(reversed(a)), width,
                                 [-orientation] * len(a))
                           @ _band(list(reversed(b)), width,
                                   [-orientation] * len(b))
                           @ _band(list(reversed(c)), width,
                                   [-orientation] * len(c)))
                np.testing.assert_allclose(
                    inverse, source.conj().T, atol=2e-14)
                np.testing.assert_allclose(
                    inverse @ source, np.eye(1 << width), atol=2e-14)

    def test_good_middle_elimination_with_non_c3_intermediate(self):
        width = 3
        # An actual paired core on wires 0,1 and an additional cancellation
        # on wire 2. Here B_0 does not commute A, so the broader reduction
        # must not assume that E_B0 E_A is C3.
        a = [2, 10, 32]
        b_bad, b_good = [16, 24], [4]
        c = [2, 10, 4]
        a_zero, a_strip = a[:2], a[2:]
        for axes in (a, b_bad + b_good, c):
            self.assertEqual(_rank(axes), len(axes))
            self.assertTrue(all(_pair(p, q, width) == 0
                                for p, q in itertools.combinations(axes, 2)))
        self.assertTrue(all(any(_pair(p, q, width) for q in c)
                            for p in b_bad))
        self.assertTrue(any(_pair(p, q, width)
                            for p in b_bad for q in a))
        middle_bad = _band(b_bad, width)
        middle_good = _band(b_good, width)
        outer = _band(c, width)
        stripped_band = _band(a_strip, width)
        intermediate = middle_bad @ _band(a, width)
        self.assertFalse(_is_c3(intermediate, width))
        source = outer @ middle_good @ intermediate
        first_strip = source @ stripped_band.conj().T
        reduced = first_strip @ middle_good.conj().T
        self.assertTrue(_is_c3(source, width))
        self.assertTrue(_is_c3(first_strip, width))
        self.assertTrue(_is_c3(reduced, width))
        self.assertEqual(_image(source, a_strip[0], width), (36, 1))
        self.assertEqual(_image(first_strip, b_good[0], width), (4, 0))
        np.testing.assert_allclose(
            reduced, outer @ middle_bad @ _band(a_zero, width), atol=2e-14)
        np.testing.assert_allclose(
            source, middle_good @ reduced @ stripped_band, atol=2e-14)
        np.testing.assert_allclose(
            source.conj().T,
            stripped_band.conj().T @ reduced.conj().T @ middle_good.conj().T,
            atol=2e-14)


if __name__ == '__main__':
    unittest.main()
