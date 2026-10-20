"""A complete one-clean query can be promoted despite nonmonomial off-code action.

All native columns, literal phases and actual inverses are checked.  The
bounded sector-cover and one-signal recovery receipts do not supply a
faster source or an asymptotic depth result.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import CI, CO, CQ2, CZ, Q2, adjoint, native_columns
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _ccz, _counts, _cz, _x


class DirtyPromotionOfCleanQueriesTests(unittest.TestCase):
    def test_two_clean_compression_recovery_pays_flag_and_literal_correction(self):
        # x=0, arbitrary old output T=1, clean A=(2,3), new flag c=4.
        # Only A0 mixes in U.  The recovery phase must nevertheless test
        # BOTH A bits; on A1=1 its raw action remains nonmonomial.
        source = [('H', 2), ('CX', 0, 1)]
        zero_test = (_x(2) + _x(3) + [('H', 4)] + _ccz((2, 3, 4))
                     + [('H', 4)] + _x(3) + _x(2))
        phase = zero_test + [('S', 4)] + _inverse(zero_test)
        raw = source + phase + _inverse(source) + phase + source
        half_root = Q2(0, Fraction(1, 2))
        omega = CQ2(half_root, half_root)
        omega_cubed = CQ2(-half_root, half_root)
        inverse_root = CQ2(half_root)

        for orientation in (zero_test, _inverse(zero_test)):
            self.assertEqual(_counts(orientation), (7, 28))
        self.assertEqual(_counts(raw), (28, 120))
        self.assertEqual(_counts(_inverse(raw)), (28, 124))
        raw_matrix = native_columns(5, raw)
        self.assertEqual(native_columns(5, _inverse(raw)), adjoint(raw_matrix))

        # Exact full-column formula, including every unpromised A,c input.
        for column in range(32):
            query_row = column ^ (2 if column & 1 else 0)
            a0, a1, flag = bool(column & 4), bool(column & 8), bool(column & 16)
            for row in range(32):
                expected = CZ
                if not a1:
                    scalar = ((omega * (-CI if a0 else CO)) if flag else
                              (omega_cubed * (CI if a0 else CO)))
                    if row == query_row:
                        expected = scalar
                elif row == query_row or row == (query_row ^ 4):
                    # On A1=1, the full word is Z_c H_A0 tensor Q.
                    sign = -1 if flag else 1
                    if a0 and row & 4:
                        sign = -sign
                    expected = inverse_root * sign
                self.assertEqual(raw_matrix[row][column], expected)

        # F=(S†H)^9 acts as omega^{-3}I on an existing wire.  Its
        # actual inverse combines adjacent S powers, with no new phase.
        correction = ([('H', 4)] + [('S', 4)] * 3) * 9
        self.assertEqual(_counts(correction), (0, 36))
        self.assertEqual(_counts(_inverse(correction)), (0, 18))
        corrected = raw + correction
        self.assertEqual(_counts(corrected), (28, 156))
        self.assertEqual(_counts(_inverse(corrected)), (28, 142))
        corrected_matrix = native_columns(5, corrected)
        self.assertEqual(native_columns(5, _inverse(corrected)), adjoint(corrected_matrix))
        for column in range(4):
            query_row = column ^ (2 if column & 1 else 0)
            for row in range(32):
                self.assertEqual(corrected_matrix[row][column],
                                 CO if row == query_row else CZ)

    def test_phase_s_recovery_keeps_literal_source_phase_and_actual_inverse(self):
        # Logical x=0 and arbitrary old output T=1, signal a=2.
        # U=H_a tensor CX(x,T) has the complete block Q/sqrt(2).
        source = [('H', 2), ('CX', 0, 1)]
        source_matrix = native_columns(3, source)
        half_root = Q2(0, Fraction(1, 2))
        inverse_root = CQ2(half_root)
        omega = CQ2(half_root, half_root)
        for column in range(4):
            expected = column ^ (2 if column & 1 else 0)
            for row in range(4):
                self.assertEqual(source_matrix[row][column],
                                 inverse_root if row == expected else CZ)

        # Chronologically U,S,U^dagger,S,U.  No scalar correction is
        # emitted: the entire word is omega S_a^dagger tensor Q, and
        # the clean-code phase omega cancels in the promotion conjugations.
        recovered = source + [('S', 2)] + _inverse(source) + [('S', 2)] + source
        recovered_matrix = native_columns(3, recovered)
        inverse = _inverse(recovered)
        self.assertEqual(_counts(recovered), (0, 8))
        self.assertEqual(_counts(inverse), (0, 12))
        self.assertEqual(native_columns(3, inverse), adjoint(recovered_matrix))
        for column in range(8):
            expected = column ^ (2 if column & 1 else 0)
            phase = omega * (-CI if column & 4 else CO)
            for row in range(8):
                self.assertEqual(recovered_matrix[row][column],
                                 phase if row == expected else CZ)

    def test_complete_sector_cover_removes_nonmonomial_off_code_action(self):
        # x=(0,1), old output T=2, source-clean A=3, new output Z=4,
        # arbitrary toggle d=5.  The entire promoted word is zero clean.
        # B = S H T H S^dagger rotates Z to H, with its literal phase.
        rotation = [('S', 0)] * 3 + [('H', 0), ('T', 0),
                                     ('H', 0), ('S', 0)]
        controlled_h = _inverse(rotation) + _cz(3, 0) + rotation
        source = [('H', 2)] + _ccz((0, 1, 2)) + [('H', 2)] + controlled_h
        source_matrix = native_columns(4, source)
        self.assertEqual(native_columns(4, _inverse(source)), adjoint(source_matrix))
        for orientation in (source, _inverse(source)):
            self.assertEqual(_counts(orientation), (9, 27))
            self.assertTrue(all(gate[0] in ('H', 'S', 'CX', 'T', 'TDG')
                                for gate in orientation))
        # All clean-code columns implement the full dirty-output query.
        for column in range(8):
            expected = column ^ (4 if column & 3 == 3 else 0)
            for row in range(16):
                self.assertEqual(source_matrix[row][column],
                                 CO if row == expected else CZ)
        # Off-code action is nonmonomial and cannot be assumed identity.
        self.assertEqual(sum(source_matrix[row][8] != CZ for row in range(16)), 2)
        self.assertNotEqual(source_matrix, adjoint(source_matrix))

        zero_test = _x(5) + [('CX', 3, 5)]
        controlled_probe = _ccz((5, 2, 4))
        conjugated_probe = source + controlled_probe + _inverse(source)

        def echo(word):
            # Chronological C,R,C^dagger,R^dagger, hence operator R†C†RC.
            return zero_test + word + _inverse(zero_test) + _inverse(word)

        selective = echo(controlled_probe) + echo(conjugated_probe)
        self.assertEqual(_counts(selective), (64, 168))
        selective_matrix = native_columns(6, selective)
        for column in range(64):
            phase = -CO if (not column & 8 and column & 3 == 3 and column & 16) else CO
            for row in range(64):
                self.assertEqual(selective_matrix[row][column], phase if row == column else CZ)

        # Cyclic one-bit Gray cover: S0, X_A, S0, X_A.  H_Z converts
        # the complete phase into a query on arbitrary Z, returning T,A,d.
        promoted = ([('H', 4)] + selective + _x(3) + selective
                    + _x(3) + [('H', 4)])
        for orientation in (promoted, _inverse(promoted)):
            self.assertEqual(_counts(orientation), (128, 346))
            matrix = native_columns(6, orientation)
            for column in range(64):
                expected = column ^ (16 if column & 3 == 3 else 0)
                for row in range(64):
                    self.assertEqual(matrix[row][column], CO if row == expected else CZ)


if __name__ == '__main__':
    unittest.main()
