"""Exact four-active-axis guard for the all-active pair normal form.

One five-wire source retains the literal central-controlled quarter-turn,
both signed-swap conventions, both clean sectors and actual inverse. It
is a structural fixture, not a growing Boolean family or degree witness.
"""
import unittest

from tests.exact_matrices import (
    CO, CI, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_original_middle_pair_insertions import _ghz, _pair_coordinates
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _conjugate, _tau


def _cz(first, second):
    return [('H', second), ('CX', first, second), ('H', second)]


def _quarter_turn(pauli, orientation):
    phase = CI if orientation == 1 else -CI
    return add(scale(eye(len(pauli)), (CO + phase) / 2),
               scale(pauli, (CO - phase) / 2))


class AllActiveSingletonClosureTests(unittest.TestCase):
    def test_four_active_axes_literal_central_control_and_inverse(self):
        width, size = 5, 32
        identity = eye(size)
        omega = native_columns(1, [('T', 0)])[1][1]
        active_word = [('T', wire) for wire in range(4)]
        active = native_columns(width, active_word)
        encoder_word = _ghz(4)
        encoder = native_columns(width, encoder_word)
        clean = _tau(1 << width, width)
        pa = _tau(15, width)
        b = _tau(15 | (16 << width), width)
        d = _tau(15 | (31 << width), width)
        gamma = _tau(15 << width, width)
        self.assertEqual(_conjugate(encoder, clean), pa)
        self.assertEqual(multiply(b, d), gamma)
        j = _conjugate(active, pa)
        self.assertEqual(_conjugate(j, b), d)
        self.assertEqual(_rank([15 | (16 << width), 15 | (31 << width)]), 2)

        # Q_R has sign q=1+r0+r1+r0*r1, rj=xj+x3. Its quadratic
        # polar matrix on the two residual central bits has rank two.
        parity = [('CX', 3, 0), ('CX', 3, 1)]
        q_word = (parity + [('S', 0), ('S', 0), ('S', 1), ('S', 1)]
                  + _cz(0, 1) + _inverse(parity))
        q = scale(native_columns(width, q_word), -1)
        # The mixed derivative on x3=0, x2=x0+x1 is nonzero;
        # its two-variable alternating polar matrix therefore has rank two.
        self.assertEqual(q[0][0] * q[3][3], -q[5][5] * q[6][6])
        self.assertEqual(multiply(q, q), identity)
        self.assertEqual(_conjugate(q, pa), pa)
        for column in range(size):
            r0 = ((column >> 0) ^ (column >> 3)) & 1
            r1 = ((column >> 1) ^ (column >> 3)) & 1
            expected = (-1) ** (1 ^ r0 ^ r1 ^ (r0 & r1))
            self.assertEqual(q[column][column], CO * expected)

        plus = scale(add(identity, gamma), CO / 2)
        minus = scale(add(identity, scale(gamma, -1)), CO / 2)
        h = _conjugate(adjoint(active), b)
        qh = multiply(q, b)  # beta=0, so the actual P_h is B.
        self.assertEqual(multiply(plus, h), multiply(plus, qh))
        self.assertNotEqual(h, qh)  # Do not extend H across Gamma=-1.

        base_coordinates = _pair_coordinates(4) + _cz(0, 4)
        x1 = [('H', 1), ('S', 1), ('S', 1), ('H', 1)]
        for swap_sign, epsilon in ((1, 1), (-1, -1)):
            with self.subTest(s=swap_sign, epsilon=epsilon):
                coordinates = ([] if swap_sign == 1 else x1) + base_coordinates
                frame = native_columns(width, coordinates)
                signed_d = scale(d, swap_sign)
                self.assertEqual(_conjugate(frame, _tau(1 << width, width)), b)
                self.assertEqual(_conjugate(frame, _tau(2 << width, width)), signed_d)
                self.assertEqual(_conjugate(j, b), scale(signed_d, swap_sign))
                pair_band = [('T' if power == 1 else 'TDG', wire)
                             for wire, power in enumerate((epsilon, swap_sign * epsilon))]
                pair_word = _inverse(coordinates) + pair_band + coordinates
                pair = native_columns(width, pair_word)
                eta = (1 - swap_sign) // 2
                prefactor = CO if eta == 0 else (omega.conj() if epsilon == 1 else omega)
                branch_phase = omega if epsilon == 1 else omega.conj()
                expected_pair = scale(add(
                    multiply(plus, _quarter_turn(b, epsilon)),
                    scale(minus, branch_phase)), prefactor)
                self.assertEqual(pair, expected_pair)

                source_word = active_word + pair_word + _inverse(active_word)
                source = native_columns(width, source_word)
                expected_source = scale(add(
                    multiply(plus, _quarter_turn(qh, epsilon)),
                    scale(minus, branch_phase)), prefactor)
                self.assertEqual(source, expected_source)
                self.assertEqual(multiply(minus, source),
                                 scale(minus, prefactor * branch_phase))
                self.assertEqual(native_columns(width, _inverse(source_word)),
                                 adjoint(source))
                self.assertEqual(multiply(adjoint(source), source), identity)

                # All 32 physical columns, including the unused clean block.
                complete_word = encoder_word + source_word + _inverse(encoder_word)
                complete = native_columns(width, complete_word)
                self.assertEqual(complete, _conjugate(adjoint(encoder), expected_source))
                self.assertEqual(native_columns(width, _inverse(complete_word)),
                                 adjoint(complete))
                self.assertEqual(_conjugate(complete, clean), clean)
                self.assertEqual([len(band) for band in _native_t_bands(complete_word)],
                                 [4, 2, 4])
                # The occupied h direction is the arbitrary fifth wire Z,
                # so this fixture retains a noncentral controlled quarter-turn.
                logical_h = _conjugate(adjoint(encoder), b)
                for column in range(0, size, 2):
                    self.assertEqual(logical_h[column][column],
                                     CO * (-1 if column & 16 else 1))


if __name__ == '__main__':
    unittest.main()
