"""Exact bounded receipts for the repaired source and quarter-turn screen.

Every logical column and actual inverse is retained. Small Pauli checks
support the literal coefficient identity, not an unrestricted depth bound.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts
from tests.test_signal_mixing_overlap import _overlapping_source


def _paulis(width):
    size = 1 << width
    for x in range(size):
        for z in range(size):
            phase = (CO, CI, -CO, -CI)[(x & z).bit_count() % 4]
            matrix = [[CZ for _ in range(size)] for _ in range(size)]
            for column in range(size):
                matrix[column ^ x][column] = (
                    -phase if (z & column).bit_count() % 2 else phase)
            yield matrix


def _coefficient(pauli, matrix):
    size = len(matrix)
    return sum((pauli[row][column].conj() * matrix[row][column]
                for row in range(size) for column in range(size)), CZ) / size


class DressedSignalOverlapTests(unittest.TestCase):
    def test_involutory_repair_keeps_native_inverse_phase_and_return(self):
        k, width, size = 3, 4, 8
        encoder_a, encoder_b, source = _overlapping_source(k)
        # Chronological S,H repeats realize the literal scalar (HS)^3.
        omega_word = [('S', 0), ('H', 0)] * 3
        parity_word = [('S', j) for j in range(k) for _ in range(2)]
        dressing = parity_word + omega_word
        repaired = source + dressing
        repaired_b = _inverse(dressing) + encoder_b
        self.assertEqual(repaired, encoder_a + _inverse(repaired_b))
        self.assertEqual(_counts(repaired), (2 * (k + 1), 11 * k + 15))
        self.assertEqual(_counts(_inverse(repaired)),
                         (2 * (k + 1), 11 * k + 21))

        matrix = native_columns(width, repaired)
        self.assertEqual(native_columns(width, _inverse(repaired)),
                         adjoint(matrix))
        logical = scale([row[:size] for row in matrix[:size]], CQ2(Q2(0, 1)))
        self.assertEqual(logical, adjoint(logical))
        self.assertEqual(multiply(logical, logical), eye(size))
        original = native_columns(width, source)
        original_logical = scale([row[:size] for row in original[:size]],
                                 CQ2(Q2(0, 1)))
        half_root = Q2(0, Fraction(1, 2))
        omega = CQ2(half_root, half_root)
        parity = native_columns(k, parity_word)
        self.assertEqual(logical, scale(multiply(parity, original_logical), omega))

        recovered = (repaired + [('S', k)] + _inverse(repaired)
                     + [('S', k)] + repaired)
        self.assertEqual(_counts(recovered), (6 * (k + 1), 33 * k + 53))
        self.assertEqual(_counts(_inverse(recovered)),
                         (6 * (k + 1), 33 * k + 63))
        recovered_matrix = native_columns(width, recovered)
        for row in range(2 * size):
            for column in range(size):
                expected = omega * logical[row][column] if row < size else CZ
                self.assertEqual(recovered_matrix[row][column], expected)
        self.assertEqual(native_columns(width, _inverse(recovered)),
                         adjoint(recovered_matrix))

    def test_quarter_turn_coefficients_include_nonzero_trace_clifford(self):
        omega_word = [('S', 0), ('H', 0)] * 3
        controlled_s = [('T', 0), ('T', 1), ('CX', 0, 1),
                        ('TDG', 1), ('CX', 0, 1)]
        cz_word = [('H', 1), ('CX', 0, 1), ('H', 1)]
        # E Z E†=H. The literal phase in E cancels in E S E†.
        encoder = [('S', 0)] * 3 + [('H', 0), ('T', 0), ('H', 0), ('S', 0)]
        about_h = _inverse(encoder) + [('S', 0)] + encoder + _inverse(omega_word)
        about_cz = controlled_s + _inverse(omega_word)
        half = CQ2(Q2(Fraction(1, 2)))
        inverse_root = CQ2(Q2(0, Fraction(1, 2)))
        fixtures = ((1, [('H', 0)], about_h, Q2(Fraction(1, 2))),
                    (2, cz_word, about_cz, Q2(Fraction(1, 4))))
        for width, reflection_word, word, squared_ceiling in fixtures:
            with self.subTest(width=width):
                size = 1 << width
                reflection = native_columns(width, reflection_word)
                matrix = native_columns(width, word)
                self.assertEqual(matrix,
                                 scale(add(eye(size), scale(reflection, -CI)),
                                       inverse_root))
                self.assertEqual(native_columns(width, _inverse(word)), adjoint(matrix))
                if width == 2:
                    self.assertEqual(sum((reflection[j][j] for j in range(size)), CZ),
                                     CO * 2)
                paulis = list(_paulis(width))
                for pauli in paulis:
                    transported = multiply(multiply(matrix, pauli), adjoint(matrix))
                    reflected = multiply(multiply(reflection, pauli), reflection)
                    cross = add(scale(multiply(pauli, reflection), CI),
                                scale(multiply(reflection, pauli), -CI))
                    self.assertEqual(transported,
                                     scale(add(add(pauli, reflected), cross), half))
                    coefficients = [_coefficient(probe, transported) for probe in paulis]
                    nonzero = [value for value in coefficients if value != CZ]
                    if len(nonzero) == 1:
                        self.assertIn(nonzero[0], (CO, -CO))
                    else:
                        self.assertTrue(all(value.abs2() <= squared_ceiling
                                            for value in nonzero))


if __name__ == '__main__':
    unittest.main()
