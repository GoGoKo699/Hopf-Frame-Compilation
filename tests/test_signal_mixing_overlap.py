"""Bounded exact receipts for a genuinely overlapping two-band source.

The k=1 and k=3 fixtures retain every logical column, the literal phase,
actual inverses and complete signal return. They do not assert a faster
indicator, minimal T-depth, or an unrestricted native-depth obstruction.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts


def _overlapping_source(k):
    """Chronological positive-S native words; signal a=k is the high bit."""
    signal = k
    parity = [('CX', j, signal) for j in range(k)]
    phase = parity + [('S', signal)] + parity
    band = ([('H', j) for j in range(k + 1)]
            + [('T', j) for j in range(k + 1)]
            + [('H', j) for j in range(k + 1)])
    encoder_a = parity + band
    encoder_b = [('H', signal)] + _inverse(phase) + band + phase
    return encoder_a, encoder_b, encoder_a + _inverse(encoder_b)


def _quadratic_x_involution(k):
    size = 1 << k
    diagonal = [[CZ for _ in range(size)] for _ in range(size)]
    linear, constant = {
        1: (1, 0), 3: (0, 0), 5: (1, 1), 7: (0, 1),
    }[k % 8]
    for z in range(size):
        weight = z.bit_count()
        value = (weight * (weight - 1) // 2
                 + linear * weight + constant) % 2
        diagonal[z][z] = -CO if value else CO
    hadamards = native_columns(k, [('H', j) for j in range(k)])
    return multiply(multiply(hadamards, diagonal), hadamards)


class SignalMixingOverlapTests(unittest.TestCase):
    def test_complete_transition_reflections_inverse_and_return(self):
        root = CQ2(Q2(0, 1))
        half = CQ2(Q2(Fraction(1, 2)))
        half_root = Q2(0, Fraction(1, 2))
        omega = CQ2(half_root, half_root)
        for k in (1, 3):
            with self.subTest(k=k):
                width, size = k + 1, 1 << k
                encoder_a, encoder_b, source = _overlapping_source(k)
                matrix = native_columns(width, source)
                self.assertEqual(_counts(source), (2 * width, 9 * k + 9))
                self.assertEqual(_counts(_inverse(source)),
                                 (2 * width, 9 * k + 9))
                self.assertEqual(native_columns(width, _inverse(source)),
                                 adjoint(matrix))

                # Every clean-code logical column, without quotienting a phase.
                block = [row[:size] for row in matrix[:size]]
                logical = scale(block, root)
                reflection = _quadratic_x_involution(k)
                expected = scale(add(eye(size), scale(reflection, -CI)),
                                 (CO - CI) * half)
                self.assertEqual(logical, expected)
                self.assertEqual(multiply(logical, logical),
                                 scale(reflection, -CO))
                self.assertEqual(multiply(adjoint(logical), logical), eye(size))

                clean_z = native_columns(width, [('Z', k)])
                image_a = native_columns(width, encoder_a)
                image_b = native_columns(width, encoder_b)
                reflection_a = multiply(multiply(image_a, clean_z),
                                        adjoint(image_a))
                reflection_b = multiply(multiply(image_b, clean_z),
                                        adjoint(image_b))
                self.assertEqual(add(multiply(reflection_a, reflection_b),
                                     multiply(reflection_b, reflection_a)),
                                 scale(eye(2 * size), CZ))

                recovered = (source + [('S', k)] + _inverse(source)
                             + [('S', k)] + source)
                self.assertEqual(_counts(recovered),
                                 (6 * width, 27 * k + 29))
                self.assertEqual(_counts(_inverse(recovered)),
                                 (6 * width, 27 * k + 33))
                recovered_matrix = native_columns(width, recovered)
                for row in range(2 * size):
                    for column in range(size):
                        expected = (omega * logical[row][column]
                                    if row < size else CZ)
                        self.assertEqual(recovered_matrix[row][column], expected)
                self.assertEqual(native_columns(width, _inverse(recovered)),
                                 adjoint(recovered_matrix))


if __name__ == '__main__':
    unittest.main()
