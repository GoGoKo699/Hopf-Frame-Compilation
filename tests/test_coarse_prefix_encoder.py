"""Complete prefix layouts, exact ledgers, and leaked-output conjugation checks."""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.tree_boundary import coarse_encoder_ledger, prefix_encoder


ATOL = 2e-11


def _prefix_frame(depth, coins):
    frame = np.eye(1 << depth, dtype=complex)
    for stage in range(depth):
        bit = depth - stage - 1
        for address in range(1 << stage):
            first = address << (depth - stage)
            second = first | (1 << bit)
            frame[[first, second]] = coins[(1 << stage) + address - 1] @ frame[[first, second]]
    return frame


def _delimiter(height):
    destinations = [0]
    for node in range(1, 2 << height):
        depth = node.bit_length() - 1
        prefix = node - (1 << depth)
        destinations.append((2 * prefix + 1) << (height - depth))
    return destinations


def _simultaneous_layers(height, coins):
    joint = np.eye(2 << height, dtype=complex)
    for stage in range(height):
        bit = height - stage
        for address in range(1 << stage):
            for suffix in range(1 << bit):
                if suffix.bit_count() == 1:
                    first = (address << (bit + 1)) | suffix
                    second = first | (1 << bit)
                    joint[[first, second]] = coins[(1 << stage) + address - 1] @ joint[[first, second]]
    destinations = _delimiter(height)
    return joint[np.ix_(destinations, destinations)]


class CoarsePrefixEncoderTests(unittest.TestCase):
    def test_simultaneous_layout_including_generic_complex_coins(self):
        rng = np.random.default_rng(1701)
        for height in range(1, 6):
            for family in ('canonical', 'generic'):
                with self.subTest(height=height, family=family):
                    table = np.resize([0., .3, -1.6, .9, 0., 1.2], (1 << height) - 1)
                    coins = []
                    for value in table:
                        if family == 'generic':
                            coin, _ = np.linalg.qr(rng.normal(size=(2, 2))
                                                  + 1j * rng.normal(size=(2, 2)))
                        else:
                            g, h = 1 / (1 + 1j * value), value / (1 + 1j * value)
                            coin = np.array([[1j * h, -g.conjugate()],
                                             [g, -1j * h.conjugate()]])
                        coins.append(coin)
                    expected = np.eye(2 << height, dtype=complex)
                    for depth in range(height + 1):
                        start, stop = 1 << depth, 1 << (depth + 1)
                        expected[start:stop, start:stop] = _prefix_frame(depth, coins)
                    actual = _simultaneous_layers(height, coins)
                    np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(actual.conj().T @ actual, np.eye(2 << height),
                                               atol=ATOL, rtol=0)
                    if family == 'canonical':
                        np.testing.assert_allclose(actual, prefix_encoder(table), atol=ATOL, rtol=0)

    def test_delimiter_by_conditioned_bit_reversals(self):
        for height in range(1, 9):
            native = [0]
            for source in range(1, 2 << height):
                bits = list(format(source, f"0{height + 1}b")[::-1])
                last_one = len(bits) - 1 - bits[::-1].index('1')
                bits[:last_one] = bits[:last_one][::-1]
                native.append(int(''.join(bits), 2))
            self.assertEqual(native, _delimiter(height))
            self.assertEqual(sorted(native), list(range(2 << height)))

    def test_weight_one_minterms_on_occupied_target_and_signal(self):
        for length in range(1, 11):
            for suffix in range(1 << length):
                for target in (0, 1):
                    for signal in (0, 1):
                        actual = target
                        for position in range(length):
                            actual ^= int(signal == 1 and suffix == 1 << position)
                        self.assertEqual(actual, target ^ int(signal == 1 and suffix.bit_count() == 1))

    def test_geometric_full_operator_precision_and_workspace(self):
        for height in range(3, 65):
            size = 1 << height
            ledger = coarse_encoder_ledger(height)
            bits = (size + height - 1) // height
            self.assertEqual(ledger.coarse_bits, bits)
            self.assertEqual(ledger.total_table_rows, size - 1)
            self.assertLessEqual(ledger.dirty_qubits, size + height + 7)
            for stage, width in enumerate(ledger.source_widths):
                self.assertEqual(width, bits + height - stage + 7)
                self.assertEqual(width + 1 + stage + 1, ledger.dirty_qubits)
            self.assertEqual(ledger.error_ratio, Fraction(43, 64) * (1 - Fraction(1, size)))
            self.assertLess(ledger.error_ratio, Fraction(43, 64))
        for height in (-1, 0, 1, 2):
            with self.assertRaises(ValueError):
                coarse_encoder_ledger(height)

    def test_small_conjugation_includes_flag_leakage_and_middle_error(self):
        rng = np.random.default_rng(314159)
        insert = np.eye(4, dtype=complex)[:, :2]
        encoder, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
        extended = np.eye(4, dtype=complex)
        extended[:2, :2] = encoder
        hermitian = np.array([[.1, .2j, .7, .3j], [-.2j, -.2, .1j, .4],
                              [.7, -.1j, .3, .2], [-.3j, .4, .2, -.4]])
        eigenvalues, vectors = np.linalg.eigh(hermitian)
        perturbation = (vectors * np.exp(1j * .07 * eigenvalues)) @ vectors.conj().T
        actual = perturbation @ extended
        delta = np.linalg.norm(actual @ insert - insert @ encoder, 2)
        self.assertGreater(np.linalg.norm(actual[2:, :2]), 0)
        for angle in (.001, .03, .4):
            middle = np.diag([np.exp(1j * angle), np.exp(-1j * angle), 1, 1])
            logical = middle[:2, :2]
            kappa = np.linalg.norm(middle - np.eye(4), 2)
            lhs = actual.conj().T @ middle @ actual @ insert - insert @ encoder.conj().T @ logical @ encoder
            rhs = (actual.conj().T @ (middle - np.eye(4)) @ (actual @ insert - insert @ encoder)
                   + (actual.conj().T @ insert - insert @ encoder.conj().T)
                   @ (logical - np.eye(2)) @ encoder)
            np.testing.assert_allclose(lhs, rhs, atol=ATOL, rtol=0)
            self.assertLessEqual(np.linalg.norm(lhs, 2), 2 * kappa * delta + ATOL)
            approximate_middle = perturbation @ middle
            epsilon = np.linalg.norm(approximate_middle - middle, 2)
            output = (actual.conj().T @ approximate_middle @ actual @ insert
                      - insert @ encoder.conj().T @ logical @ encoder)
            self.assertLessEqual(np.linalg.norm(output, 2), epsilon + 2 * kappa * delta + ATOL)


if __name__ == '__main__':
    unittest.main()
