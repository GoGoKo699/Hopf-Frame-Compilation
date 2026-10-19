"""Exact native checks for a conditional non-Pauli clean-code decoder.

The decoder uses no extra wire. Its C3 description does not make its
conditional Clifford action constant rank or supply a free native band.
Every fixture retains all physical columns and its actual inverse word.
"""
import unittest

from tests.exact_matrices import (
    CZ, adjoint, eye, multiply, native_columns,
)
from tests.test_clean_code_clifford_rigidity import _pauli_label, _rank
from tests.test_projected_band_matching import _cq_root
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _primitive_inverse
from tests.test_two_layer_galois_return import _conjugate, _tau


def _decoder_word(k):
    """Wire 0 is the quantum core; all other wires are arbitrary controls."""
    e = k % 2
    fanout = [('CX', 0, bit) for bit in range(1, k)]
    first = [('TDG', bit) for bit in range(1, k)]
    second = ([('T', 0)] if not e else []) + [('T', bit) for bit in range(1, k)]
    return first + list(reversed(fanout)) + second + fanout


def _controlled_s_bank(m):
    """Literal product CS_(0,i), with two parallel native T bands."""
    fanout = [('CX', 0, bit) for bit in range(1, m + 1)]
    clifford = [('S', 0)] * (m // 2)
    first = ([('T', 0)] if m % 2 else []) + [('T', bit) for bit in range(1, m + 1)]
    return clifford + first + fanout + [('TDG', bit) for bit in range(1, m + 1)] + list(reversed(fanout))


class NonPauliCodeDecodingTests(unittest.TestCase):
    def test_conditional_decoder_matches_the_literal_canonical_band(self):
        for k in (2, 3):
            with self.subTest(k=k):
                e = k % 2
                m = k - 1
                word = _decoder_word(k)
                d = native_columns(k, word)
                self.assertEqual(native_columns(k, _primitive_inverse(word)), adjoint(d))
                self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 2 * k - 1 - e)
                self.assertEqual(sum(g[0] == 'CX' for g in word), 2 * m)
                self.assertEqual(len(_native_t_bands(word)), 2)
                # Actual active axes Z_t and Z_t Z_i are independent.
                # Fan-in computes their parities on the satellite wires.
                fanout = [('CX', 0, bit) for bit in range(1, k)]
                active_word = fanout + [('T', bit) for bit in range(k)] + list(reversed(fanout))
                active = native_columns(k, active_word)
                self.assertEqual(native_columns(k, _primitive_inverse(active_word)), adjoint(active))
                for column in range(1 << k):
                    target = column & 1
                    weight = (column >> 1).bit_count()
                    r = (k - e) // 2 - weight
                    for row in range(1 << k):
                        expected_d = _cq_root(2 * target * r) if row == column else CZ
                        expected_active = _cq_root(weight + 2 * target * r + e * target) if row == column else CZ
                        self.assertEqual(d[row][column], expected_d)
                        self.assertEqual(active[row][column], expected_active)
                # E_a=omega^w D T_t^e is literal on both core sectors.
                central_phase = native_columns(k, [('T', bit) for bit in range(1, k)])
                parity_t = native_columns(k, [('T', 0)] if e else [])
                self.assertEqual(active, multiply(central_phase, multiply(d, parity_t)))

    def test_decoder_is_a_paid_two_band_operation_on_arbitrary_controls(self):
        # No central control is initialized. Both actual inverse schedules
        # keep the same finite wire set and native T-depth.
        for k in (2, 3, 6, 7):
            word = _decoder_word(k)
            for actual in (word, _primitive_inverse(word)):
                self.assertEqual(sum(g[0] in ('T', 'TDG') for g in actual), 2 * k - 1 - k % 2)
                self.assertEqual(sum(g[0] == 'CX' for g in actual), 2 * (k - 1))
                self.assertEqual({wire for gate in actual for wire in gate[1:]}, set(range(k)))
                for band in _native_t_bands(actual):
                    self.assertEqual(len({gate[1] for gate in band}), len(band))

    def test_conditional_clifford_displacement_rank_grows_with_control_width(self):
        for m in (1, 2, 3, 4):
            with self.subTest(m=m):
                width = m + 1
                word = _controlled_s_bank(m)
                d = native_columns(width, word)
                self.assertEqual(native_columns(width, _primitive_inverse(word)), adjoint(d))
                for column in range(1 << width):
                    exponent = 2 * (column & 1) * (column >> 1).bit_count()
                    for row in range(1 << width):
                        self.assertEqual(d[row][column], _cq_root(exponent) if row == column else CZ)
                x = _tau(1, width)
                image = _conjugate(d, x)
                # This independently emitted Clifford fixes every literal
                # sign: J=X_t (product S_i)(product CZ_(t,i)).
                clifford_word = [('S', bit) for bit in range(1, width)]
                clifford_word += [gate for bit in range(1, width)
                                  for gate in [('H', bit), ('CX', 0, bit), ('H', bit)]]
                clifford_word += [('H', 0), ('S', 0), ('S', 0), ('H', 0)]
                self.assertEqual(image, native_columns(width, clifford_word))
                self.assertEqual(multiply(image, image), eye(1 << width))
                columns = []
                for bit in range(2 * width):
                    label = _pauli_label(_conjugate(image, _tau(1 << bit, width)), width)
                    self.assertIsNotNone(label)
                    columns.append(label)
                observed_rank = _rank([label ^ (1 << bit)
                                       for bit, label in enumerate(columns)])
                self.assertEqual(observed_rank, m + m % 2)


if __name__ == '__main__':
    unittest.main()
