"""Finite full-space identities used by the mask and sparse zero-test reductions."""

from itertools import product
import unittest

import numpy as np


def shift_mask(mask, offset, size):
    return sum(((mask >> (index ^ offset)) & 1) << index for index in range(size))


def binary_convolution(left, right, size):
    result = 0
    for index in range(size):
        if (left >> index) & 1:
            result ^= shift_mask(right, index, size)
    return result


def gf8_multiply(left, right):
    """Use X^3 + X + 1 to identify three-bit coordinates with GF(8)."""
    result = 0
    while right:
        if right & 1:
            result ^= left
        left <<= 1
        if left & 8:
            left ^= 0b1011
        right >>= 1
    return result


class ReturnedMaskAndZeroTest(unittest.TestCase):
    def test_hidden_mask_four_calls_on_dirty_banks(self):
        size = 4
        for mask, address, origin in product(range(16), range(size), range(size)):
            first = shift_mask(mask, origin, size)
            second = shift_mask(mask, origin ^ address, size)
            expected = (mask.bit_count() % 2) << address
            for left, right in ((0, 0), (0b0101, 0b1011), (0b1111, 0b1111)):
                corners = ((left, right), (left ^ first, right),
                           (left ^ first, right ^ second), (left, right ^ second))
                result = 0
                for a, b in corners:
                    result ^= binary_convolution(a, b, size)
                self.assertEqual(result, expected)

    def test_validator_accepts_unit_mask_and_rejects_wrong_masks(self):
        size = 4
        for address, vector in product(range(size), range(1 << size)):
            if vector.bit_count() % 2 == 0:
                continue  # rejected by the independent parity check
            rejections = 0
            for sample, query in product(range(1 << size), range(size)):
                selected = sum(((((i ^ address) & query).bit_count() % 2) << i)
                               for i in range(size))
                rejections += (sample & selected & vector).bit_count() % 2
            if vector == 1 << address:
                self.assertEqual(rejections, 0)
            else:
                self.assertGreaterEqual(4 * rejections, size * (1 << size))

    def test_full_operator_target_echo_with_address_motion(self):
        """Scramble all address/work outputs inside each final target cylinder."""
        rng = np.random.default_rng(743)
        dim = 16
        target_z = np.diag([(-1.0) ** (index % 2) for index in range(dim)])
        raw = np.eye(dim)
        angle = 0.07
        moving = np.array([[np.sin(angle), np.cos(angle)],
                           [np.cos(angle), -np.sin(angle)]])
        raw[:2, :2] = moving
        raw[2:4, 2:4] = moving
        scramble = np.zeros((dim, dim), dtype=complex)
        for bit in range(2):
            matrix = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
            scramble[bit::2, bit::2] = np.linalg.qr(matrix)[0]
        raw = scramble @ raw
        echo = raw.conj().T @ target_z @ raw @ target_z
        np.testing.assert_allclose(echo[:, 4:], np.eye(dim)[:, 4:], atol=1e-12)
        np.testing.assert_allclose(echo[4:, :4], 0, atol=1e-12)
        self.assertAlmostEqual(np.linalg.norm(echo[:4, :4] + np.eye(4), 2),
                               2 * np.sin(angle))

    def test_affine_trial_bound_for_arbitrary_origins(self):
        size = 8
        for origins in product(range(size), repeat=2):
            singles = 0
            for slope, intercept in product(range(size), repeat=2):
                shifted = [value ^ gf8_multiply(slope, index) ^ intercept
                           for index, value in enumerate(origins)]
                singles += shifted.count(0) == 1
            self.assertGreaterEqual(16 * singles, 3 * size**2)

    def test_sparse_dilations_and_orthogonal_active_subspaces(self):
        support = {0, 1}
        sumset = {x ^ y for x in support for y in support}
        for difference in range(1, 8):
            moved = [gf8_multiply(scale, difference) for scale in (1, 2)]
            self.assertTrue(any(delta not in sumset for delta in moved))
            inactive_delta = next(delta for delta in moved if delta not in sumset)
            self.assertTrue(support.isdisjoint({x ^ inactive_delta for x in support}))

        # On disjoint supports the only active child may rotate the origin
        # within that support, even when it is entangled with its carrier.
        rng = np.random.default_rng(418)
        size = 16  # four possible origins, two carrier bits
        ident = np.eye(size, dtype=complex)
        first, second = ident.copy(), ident.copy()
        for lower, child in ((0, first), (8, second)):
            sector = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
            child[lower:lower + 8, lower:lower + 8] = np.linalg.qr(sector)[0]
        phase = np.diag([(-1) ** ((i & 3) == 3) for i in range(size)])
        word = ident.copy()
        for gate in (phase, first, phase, second, phase, first.conj().T,
                     phase, second.conj().T):
            word = gate @ word
        np.testing.assert_allclose(word, ident, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
