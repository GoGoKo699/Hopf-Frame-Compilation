"""Exact complete-column receipts for overlapping two-band source probes.

These bounded native words certify the contraction, actual inverse and
non-Clifford residual guard. They supply no asymptotic indicator compiler.
"""
from itertools import product
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_majorana_walsh_query import _inverse


def _conjugate(unitary, observable):
    return multiply(multiply(unitary, observable), adjoint(unitary))


def _clean_block(matrix):
    # The sole signal is wire zero, hence even indices are its zero sector.
    return [[matrix[row][column] for column in range(0, len(matrix), 2)]
            for row in range(0, len(matrix), 2)]


class OverlappingSourceProbesTests(unittest.TestCase):
    def test_overlapping_parity_contraction_keeps_all_columns_and_inverse(self):
        # Operator V_A=TH and V_B=HTH on the same logical wire. Their
        # noncommuting T axes overlap; no disjoint-bank merge is used.
        va = [('H', 1), ('T', 1)]
        vb = [('H', 1), ('T', 1), ('H', 1)]
        ua = [('CX', 1, 0)] + va
        ub = [('CX', 1, 0), ('H', 0)] + vb
        source = ua + _inverse(ub)
        unitary = native_columns(2, source)
        self.assertEqual(native_columns(2, _inverse(source)), adjoint(unitary))

        w = native_columns(1, [('H', 0), ('T', 0)]
                           + _inverse([('H', 0), ('T', 0), ('H', 0)]))
        phi_w = [[entry * (-CO if row & column else CO)
                  for column, entry in enumerate(values)]
                 for row, values in enumerate(w)]
        half_root = native_columns(1, [('H', 0)])[0][0]
        self.assertEqual(_clean_block(unitary), scale(phi_w, half_root))
        # The contraction is complete even when this particular source is
        # not flat; success on just one input column would miss this defect.
        self.assertNotEqual(multiply(adjoint(phi_w), phi_w), eye(2))
        self.assertEqual([[entry * (-CO if row & column else CO)
                           for column, entry in enumerate(values)]
                          for row, values in enumerate(phi_w)], w)

    def test_flat_source_probe_is_unitary_but_need_not_be_clifford(self):
        # Here s=0. The source is flat, with complete normalized logical
        # transition L=V_B^dagger V_A and no hidden initialized work.
        ua = [('CX', 1, 0), ('H', 1), ('T', 1)]
        ub = [('H', 0), ('H', 1), ('T', 1), ('H', 1)]
        source = ua + _inverse(ub)
        uam, ubm = native_columns(2, ua), native_columns(2, ub)
        unitary = native_columns(2, source)
        self.assertEqual(native_columns(2, _inverse(source)), adjoint(unitary))
        za = native_columns(2, [('Z', 0)])
        zj = native_columns(2, [('Z', 1)])
        a, b = _conjugate(uam, za), _conjugate(ubm, za)
        c = _conjugate(ubm, zj)
        self.assertEqual(add(multiply(a, b), multiply(b, a)), zeros(4))
        aca = multiply(multiply(a, c), a)
        self.assertNotEqual(aca, c)
        self.assertNotEqual(aca, scale(c, -CO))
        plus = scale(add(c, aca), CO / 2)
        minus = scale(add(c, scale(aca, -CO)), CO / 2)
        probe = add(plus, multiply(minus, b))
        self.assertEqual(adjoint(probe), probe)
        self.assertEqual(multiply(probe, probe), eye(4))
        self.assertEqual(multiply(probe, a), multiply(a, probe))

        half_root = native_columns(1, [('H', 0)])[0][0]
        logical = scale(_clean_block(unitary), CO / half_root)
        self.assertEqual(multiply(adjoint(logical), logical), eye(2))
        expected = multiply(multiply(adjoint(logical),
                                     native_columns(1, [('Z', 0)])), logical)
        transported = multiply(multiply(adjoint(uam), probe), uam)
        self.assertEqual(_clean_block(transported), expected)
        # Exact signal return for the mathematical transported extension.
        self.assertTrue(all(transported[row][column] == CZ
                            for column in (0, 2) for row in (1, 3)))

        # T_j cannot be treated as a Clifford: one Pauli has a non-Pauli
        # image. Exhaust the 16 Hermitian Paulis with both literal signs.
        pauli_words = (
            [], [('H', 0), ('Z', 0), ('H', 0)],
            [('S', 0)] * 3 + [('H', 0), ('Z', 0), ('H', 0), ('S', 0)],
            [('Z', 0)],
        )
        image = _conjugate(probe, za)
        for first, second in product(pauli_words, repeat=2):
            word = first + [(gate, 1) for gate, _ in second]
            pauli = native_columns(2, word)
            self.assertNotEqual(image, pauli)
            self.assertNotEqual(image, scale(pauli, -CO))

        recovered = source + [('S', 0)] + _inverse(source) + [('S', 0)] + source
        recovered_matrix = native_columns(2, recovered)
        self.assertEqual(native_columns(2, _inverse(recovered)),
                         adjoint(recovered_matrix))
        omega = native_columns(1, [('T', 0)])[1][1]
        self.assertEqual(_clean_block(recovered_matrix), scale(logical, omega))
        self.assertTrue(all(recovered_matrix[row][column] == CZ
                            for column in (0, 2) for row in (1, 3)))


if __name__ == '__main__':
    unittest.main()
