"""Four-row native lookup with literal phases and arbitrary core inputs.

Tiny Pauli-mask kernels cover all coherent address and core columns.  At
q=5, eight complete 256-dimensional native sectors cover the full eleven-
wire rotation operator.  A flattened native word separately propagates
coherent addresses and external-reference columns.  Fine precision uses
exact coefficient certificates and gate counts, not large simulation.
"""
from dataclasses import replace
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.native_residual_lookup import (
    _four_row_mask,
    emit_four_row_residual_table,
    emit_four_row_rotation_table,
)
from compiler_robust_hopf.native_residual_rotation import _mask
from compiler_robust_hopf.residual_table_preprocessing import CoefficientInterval
from compiler_robust_hopf.rotation_programming import program_rotation
from tests.test_native_residual_rotation import _algebra_rotation, _majorana_triplet, _native_matrix
from tests.test_native_residual_table import ALPHABET, _inverse, _restrict_sector


ATOL = 1.2e-10
COEFFICIENTS = ((0, 0), (1, 0), (0, 1), (Fraction(-3, 4), Fraction(1, 2)))


def _point(value):
    value = Fraction(value)
    return CoefficientInterval(value, value)


def _signs_for_pauli(xs, zs, core=3):
    """Invert the Majorana sign-to-Pauli relation independently, bit by bit."""
    prefix = 0
    signs = []
    for wire in range(core):
        first = int(wire in zs) ^ prefix
        second = first ^ int(wire in xs)
        signs.extend((first, second))
        prefix ^= int(wire in xs)
    return tuple(signs)


def _quadratic_supports(programs):
    # The degree-two ANF coefficient is the parity of the four truth rows.
    # Recover literal row masks rather than using the emitter's ANF helper.
    masks = tuple(_mask(program.signs) for program in programs)
    result = []
    for axis in ("X", "Z"):
        parity = set()
        for mask in masks:
            parity.symmetric_difference_update(wires[0] for name, wires in mask if name == axis)
        result.append(tuple(sorted(parity)))
    return tuple(result)


def _apply_complete_sectors(sectors, columns):
    result = np.empty_like(columns)
    for enable in (0, 1):
        for row in range(4):
            block = (4 * enable + row) * 256
            result[block:block + 256] = sectors[row, enable] @ columns[block:block + 256]
    return result


class NativeResidualLookupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = emit_four_row_residual_table(COEFFICIENTS, 5, enable_value=0)
        table = cls.table
        cls.native = {}
        cls.algebra_errors = []
        cache, algebra = {}, {}
        for factor, rotation in enumerate(table.rotations[:2]):
            axis = "z" if factor == 0 else "y"
            for row in range(4):
                program = rotation.programs[row]
                if (axis, program) not in algebra:
                    algebra[axis, program] = _algebra_rotation(program, _majorana_triplet(program), axis)
                for enable in (0, 1):
                    fixed = {table.addresses[0]: row & 1, table.addresses[1]: row >> 1,
                             table.enable: enable}
                    reduced, phase_power, width = _restrict_sector(rotation.gates, table.enable + 1, fixed)
                    key = reduced, phase_power
                    if key not in cache:
                        cache[key] = ((1 + 1j) / np.sqrt(2)) ** phase_power * _native_matrix(reduced, width)
                    actual = cache[key]
                    expected = algebra[axis, program] if enable == 0 else np.eye(256)
                    cls.algebra_errors.append(np.max(abs(actual - expected)))
                    cls.native[factor, row, enable] = actual
        cls.combined = {}
        for row in range(4):
            for enable in (0, 1):
                outer, middle = cls.native[0, row, enable], cls.native[1, row, enable]
                cls.combined[row, enable] = outer @ middle @ outer

    def test_pauli_anf_kernels_on_all_coherent_address_and_core_columns(self):
        # Cases cover neither, each single, and both quadratic axes, plus
        # a multiwire support. The overlapping final X/Z factor has odd
        # intersection in the last two cases, detecting order/phase errors.
        base = (((), ()), ((0,), (1,)), ((1,), (0,)))
        cases = (
            base + (((0, 1), (0, 1)),),
            base + (((0, 1, 2), (0, 1)),),
            base + (((0, 1), (0, 1, 2)),),
            base + (((0, 1, 2), (0, 1, 2)),),
            (((), ()), ((), ()), ((), ()), ((0, 2), (0,))),
        )
        for row_paulis in cases:
            with self.subTest(rows=row_paulis):
                signs = tuple(_signs_for_pauli(xs, zs) for xs, zs in row_paulis)
                mask = _four_row_mask(signs, (3, 4))
                expected = np.zeros((32, 32), dtype=complex)
                for row, row_signs in enumerate(signs):
                    expected[row * 8:(row + 1) * 8, row * 8:(row + 1) * 8] = _native_matrix(_mask(row_signs), 3)
                np.testing.assert_allclose(_native_matrix(mask.gates, 5), expected, atol=2e-14, rtol=0)
                np.testing.assert_allclose(_native_matrix(_inverse(mask.gates), 5), expected.conj().T,
                                           atol=2e-14, rtol=0)
                supports = []
                for axis in (0, 1):
                    parity = set()
                    for paulis in row_paulis:
                        parity.symmetric_difference_update(paulis[axis])
                    supports.append(tuple(sorted(parity)))
                self.assertEqual(mask.quadratic_supports, tuple(supports))
                k = sum(bool(support) for support in supports)
                self.assertEqual(mask.quadratic_axis_count, k)
                self.assertEqual(sum(name in ("T", "TDG") for name, _ in mask.gates), 7 * k)
                self.assertEqual(mask.t_count, 7 * k)
                # The native pivots lie in the arbitrary precision core;
                # there is no clean/dirty helper or target-flag allocation.
                self.assertTrue(all(0 <= wire < 5 for _, wires in mask.gates for wire in wires))

    def test_eight_complete_native_sectors_and_negative_enable(self):
        table = self.table
        self.assertEqual(table.addresses, (8, 9))
        self.assertEqual((table.signal, table.target, table.enable), (6, 7, 10))
        self.assertLess(max(self.algebra_errors), ATOL)
        self.assertIs(table.rotations[0], table.rotations[2])
        self.assertIs(table.programs[0], table.programs[2])
        self.assertEqual(table.gates, table.rotations[0].gates + table.rotations[1].gates + table.rotations[0].gates)
        for factor, rotation in enumerate(table.rotations[:2]):
            self.assertEqual(rotation.quadratic_axis_count, 2)
            self.assertEqual(rotation.quadratic_supports, _quadratic_supports(rotation.programs))
            positive = emit_four_row_rotation_table(rotation.programs, "z" if factor == 0 else "y")
            wrapper = (("X", (table.enable,)),)
            self.assertEqual(rotation.gates, wrapper + positive.gates + wrapper)
        identity = np.eye(256)
        for row in range(4):
            for enable in (0, 1):
                with self.subTest(row=row, enable=enable):
                    actual = self.combined[row, enable]
                    np.testing.assert_allclose(actual.conj().T @ actual, identity, atol=4 * ATOL, rtol=0)
                    if enable == 1:
                        np.testing.assert_allclose(actual, identity, atol=3 * ATOL, rtol=0)
                    else:
                        # A borrowed synthesis signal is an arbitrary
                        # coherent input, not a hidden clean initialization.
                        flip = np.arange(256) ^ (1 << table.signal)
                        np.testing.assert_allclose(actual[flip], actual[:, flip], atol=3 * ATOL, rtol=0)
        self.assertGreater(np.linalg.norm((self.combined[0, 0] - identity)[:, 0]), 1)

    def test_flat_native_lookup_on_coherent_address_reference_inputs_and_inverse(self):
        table = self.table
        inputs = np.zeros((2048, 2), dtype=complex)
        for block in range(8):
            inputs[256 * block + (17 * block + 1) % 256, 0] = 1j ** block / np.sqrt(8)
            inputs[256 * block + (53 * block + 128) % 256, 1] = (-1) ** block / np.sqrt(8)
        actual = _native_matrix(table.gates, table.enable + 1, inputs)
        expected = _apply_complete_sectors(self.combined, inputs)
        np.testing.assert_allclose(actual, expected, atol=4 * ATOL, rtol=0)
        np.testing.assert_allclose(_native_matrix(_inverse(table.gates), table.enable + 1, actual),
                                   inputs, atol=5 * ATOL, rtol=0)

    def test_fine_precision_bounds_layout_and_mask_dependent_t_ledger(self):
        for q in (5, 16, 80, 129):
            with self.subTest(q=q):
                table = emit_four_row_residual_table(COEFFICIENTS, q)
                self.assertEqual(table.core, tuple(range(q + 1)))
                self.assertEqual((table.signal, table.target, table.addresses, table.enable),
                                 (q + 1, q + 2, (q + 3, q + 4), q + 5))
                self.assertEqual(table.nqubits, q + 6)
                self.assertEqual(table.dirty_qubits, q + 2)
                self.assertLess(table.operator_error_bound, Fraction(130, 1 << q))
                expected_error = max(row.completion_error_bound for row in table.rotation_data)
                expected_error += sum(max(program.operator_error_bound for program in programs)
                                      for programs in table.programs)
                self.assertEqual(table.operator_error_bound, expected_error)
                for rotation in table.rotations:
                    supports = _quadratic_supports(rotation.programs)
                    k = sum(bool(support) for support in supports)
                    self.assertEqual(rotation.quadratic_supports, supports)
                    self.assertEqual(rotation.quadratic_axis_count, k)
                    self.assertEqual(rotation.operator_error_bound,
                                     max(program.operator_error_bound for program in rotation.programs))
                    self.assertEqual(sum(name in ("T", "TDG") for name, _ in rotation.mask_gates), 7 * k)
                    self.assertEqual(rotation.t_count, 180 * q + 210 + 70 * k)
                    self.assertEqual(sum(name in ("T", "TDG") for name, _ in rotation.gates), rotation.t_count)
                    self.assertFalse(any(wire in (table.signal, table.target, table.enable)
                                         for _, wires in rotation.mask_gates for wire in wires))
                kz, ky = (rotation.quadratic_axis_count for rotation in table.rotations[:2])
                self.assertEqual(table.t_count, 540 * q + 630 + 70 * (2 * kz + ky))
                self.assertLessEqual(table.t_count, 540 * q + 1050)
                self.assertEqual(sum(name in ("T", "TDG") for name, _ in table.gates), table.t_count)
                if q == 5:
                    self.assertEqual(table.t_count, 540 * q + 1050)
                for name, wires in table.gates:
                    self.assertIn(name, ALPHABET)
                    self.assertEqual(len(wires), 2 if name == "CX" else 1)
                    self.assertEqual(len(wires), len(set(wires)))
                    self.assertTrue(all(0 <= wire < table.nqubits for wire in wires))

    def test_rejects_invalid_rows_predicates_axes_and_forged_programs(self):
        program = program_rotation(_point(1), _point(0), 5)
        other = program_rotation(_point(0), _point(1), 6)
        for rows in ((program,), (program,) * 3, (program,) * 5, (program,) * 3 + (other,)):
            with self.subTest(length=len(rows)), self.assertRaises(ValueError):
                emit_four_row_rotation_table(rows)
        for predicate in (2, -1, True, False, 0., 1., Fraction(1)):
            with self.subTest(enable=predicate), self.assertRaises(ValueError):
                emit_four_row_rotation_table((program,) * 4, enable_value=predicate)
            with self.subTest(enable=predicate), self.assertRaises(ValueError):
                emit_four_row_residual_table(COEFFICIENTS, 5, enable_value=predicate)
        with self.assertRaises(ValueError):
            emit_four_row_rotation_table((program,) * 4, "x")
        with self.assertRaises(ValueError):
            emit_four_row_rotation_table((program,) * 3 + (replace(program, s=program.s + 1),))
        for rows in (COEFFICIENTS[:3], COEFFICIENTS + ((0, 0),),
                     COEFFICIENTS[:3] + ((1, 1),), COEFFICIENTS[:3] + ((Fraction(1, 3), 0),)):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                emit_four_row_residual_table(rows, 5)
        for q in (True, False, 5., Fraction(5), 4):
            with self.subTest(q=q), self.assertRaises(ValueError):
                emit_four_row_residual_table(COEFFICIENTS, q)
        for value in (True, False, .5):
            with self.subTest(coefficient=value), self.assertRaises(TypeError):
                emit_four_row_residual_table(COEFFICIENTS[:3] + ((value, 0),), 5)


if __name__ == "__main__":
    unittest.main()
