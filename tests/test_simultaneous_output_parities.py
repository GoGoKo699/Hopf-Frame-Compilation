"""Actual one-band saturation of the scalar-commuting parity screen.

All physical and logical columns retain literal phases and actual inverses.
The family is not an indicator compiler or evidence for an asymptotic bound.
"""
import unittest

from tests.exact_matrices import (
    CI, CO, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_majorana_walsh_query import _inverse


def _conjugate(unitary, matrix):
    return multiply(multiply(unitary, matrix), adjoint(unitary))


class SimultaneousOutputParityTests(unittest.TestCase):
    def test_both_scalar_parity_spaces_vanish_for_a_complete_identity(self):
        outputs = 2
        width = outputs + 1
        # C_A Z_j C_A^dagger=Z_a X_j, while C_B gives X_a X_j.
        ca = ([('CX', 0, j) for j in range(1, width)]
              + [('H', j) for j in range(1, width)])
        cb = ([('H', j) for j in range(1, width)]
              + [('CX', j, 0) for j in range(1, width)])
        band = [('T', j) for j in range(1, width)]
        ua = _inverse(ca) + band + ca
        ub = [('H', 0)] + _inverse(cb) + band + cb
        source = ua + _inverse(ub)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in source), 2 * outputs)
        self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in source),
                         8 * outputs + 1)
        uam, ubm = native_columns(width, ua), native_columns(width, ub)
        unitary = native_columns(width, source)
        self.assertEqual(native_columns(width, _inverse(source)), adjoint(unitary))
        accepted = [[unitary[row][column] for column in range(0, 1 << width, 2)]
                    for row in range(0, 1 << width, 2)]
        half_root = native_columns(1, [('H', 0)])[0][0]
        self.assertEqual(accepted, scale(eye(1 << outputs), half_root))
        za = native_columns(width, [('Z', 0)])
        xa = native_columns(width, [('H', 0), ('Z', 0), ('H', 0)])
        self.assertEqual(_conjugate(uam, za), za)
        self.assertEqual(_conjugate(ubm, za), xa)
        for encoder, opposite_clean_image in ((ubm, za), (uam, xa)):
            scalar_parities = []
            for mask in range(1 << outputs):
                parity = native_columns(width, [('Z', j + 1) for j in range(outputs)
                                                if mask & (1 << j)])
                probe = _conjugate(encoder, parity)
                dressed = multiply(multiply(opposite_clean_image, probe),
                                   opposite_clean_image)
                if dressed in (probe, scale(probe, -CO)):
                    scalar_parities.append(mask)
            self.assertEqual(scalar_parities, [0])

    def test_literal_complete_transition_support_and_every_output_parity(self):
        for outputs in (2, 3):
            with self.subTest(outputs=outputs):
                width = outputs + 1
                # Signal a=0. C maps each original Z_j to X_a X_j.
                c = ([('H', j) for j in range(1, width)]
                     + [('CX', j, 0) for j in range(1, width)])
                band = [('T', j) for j in range(1, width)]
                e = _inverse(c) + band + c
                ub = [('H', 0)] + e
                source = _inverse(ub)  # U_A=I.
                self.assertEqual(sum(g[0] in ('T', 'TDG') for g in source),
                                 outputs)
                self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in source),
                                 4 * outputs + 1)
                ubm = native_columns(width, ub)
                unitary = native_columns(width, source)
                self.assertEqual(unitary, adjoint(ubm))
                self.assertEqual(native_columns(width, _inverse(source)), ubm)

                logical_word = [gate for j in range(outputs)
                                for gate in [('H', j), ('TDG', j), ('H', j)]]
                logical = native_columns(outputs, logical_word)
                accepted = [[unitary[row][column]
                             for column in range(0, 1 << width, 2)]
                            for row in range(0, 1 << width, 2)]
                half_root = native_columns(1, [('H', 0)])[0][0]
                self.assertEqual(accepted, scale(logical, half_root))
                self.assertEqual(multiply(adjoint(logical), logical),
                                 eye(1 << outputs))

                za = native_columns(width, [('Z', 0)])
                xa = native_columns(width, [('H', 0), ('Z', 0), ('H', 0)])
                self.assertEqual(_conjugate(ubm, za), xa)
                actual_support_word = _conjugate(unitary, za)
                expected_support_word = xa
                for j in range(1, width):
                    za_xj = native_columns(
                        width, [('Z', 0), ('H', j), ('Z', j), ('H', j)])
                    factor = scale(add(eye(1 << width), scale(za_xj, -CI)),
                                   half_root)
                    expected_support_word = multiply(expected_support_word, factor)
                self.assertEqual(actual_support_word, expected_support_word)

                # Exact Pauli support is a full output-mask cube; no nonzero
                # parity has scalar commutation with the actual clean image A.
                scalar_parities = []
                for mask in range(1 << outputs):
                    parity = native_columns(
                        width, [('Z', j + 1) for j in range(outputs)
                                if mask & (1 << j)])
                    probe = _conjugate(ubm, parity)
                    dressed = multiply(multiply(za, probe), za)
                    if dressed in (probe, scale(probe, -CO)):
                        scalar_parities.append(mask)
                self.assertEqual(scalar_parities, [0])


if __name__ == '__main__':
    unittest.main()
