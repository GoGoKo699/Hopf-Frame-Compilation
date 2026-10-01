"""Literal two-row, enabled residual tables on arbitrary borrowed inputs.

At q=5, four complete 256-dimensional native sector matrices cover every
column of the ten-wire operator. Controls are restricted only gate by
gate: enable flips are tracked and its diagonal phases retained. No
1024-dimensional matrix or fine-precision statevector is constructed.
"""
from dataclasses import replace
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.native_residual_rotation import _mask, _toffoli
from compiler_robust_hopf.native_residual_table import (
    _table_mask,
    emit_residual_table,
    emit_rotation_table,
)
from compiler_robust_hopf.residual_table_preprocessing import CoefficientInterval
from compiler_robust_hopf.rotation_programming import program_rotation
from tests.test_native_residual_rotation import (
    _algebra_rotation,
    _majorana_triplet,
    _native_matrix,
)


ATOL = 9e-11
ALPHABET = {"X", "Z", "H", "S", "SDG", "T", "TDG", "CX"}
PHASE_POWERS = {"Z": 4, "S": 2, "SDG": -2, "T": 1, "TDG": -1}


def _point(value):
    value = Fraction(value)
    return CoefficientInterval(value, value)


def _inverse(gates):
    inverse_names = {"T": "TDG", "TDG": "T", "S": "SDG", "SDG": "S"}
    return tuple((inverse_names.get(name, name), wires)
                 for name, wires in reversed(gates))


def _restrict_sector(gates, width, fixed):
    """Restrict only wires that are demonstrably fixed throughout each gate.

    X on a fixed wire changes its tracked value. Its later diagonal gates
    contribute a literal scalar; no phase alignment is performed. H or a
    CNOT targeting such a wire is rejected, rather than silently ignored.
    """
    current = dict(fixed)
    free = [wire for wire in range(width) if wire not in fixed]
    relabel = {wire: index for index, wire in enumerate(free)}
    reduced = []
    phase_power = 0
    for name, wires in gates:
        if name not in ALPHABET:
            raise AssertionError(f"Unexpected elementary gate {name}.")
        if name == "CX":
            control, target = wires
            if target in fixed:
                raise AssertionError("A purported sector control is targeted by CNOT.")
            if control in fixed:
                if current[control]:
                    reduced.append(("X", (relabel[target],)))
            else:
                reduced.append((name, (relabel[control], relabel[target])))
        elif wires[0] in fixed:
            wire = wires[0]
            if name == "X":
                current[wire] ^= 1
            elif name in PHASE_POWERS:
                phase_power += current[wire] * PHASE_POWERS[name]
            else:
                raise AssertionError("A purported sector control is put in superposition.")
        else:
            reduced.append((name, (relabel[wires[0]],)))
    if current != fixed:
        raise AssertionError("The declared sector controls were not restored.")
    return tuple(reduced), phase_power % 8, len(free)


def _sector_matrix(gates, width, fixed, columns=None):
    reduced, phase_power, free_width = _restrict_sector(gates, width, fixed)
    phase = ((1 + 1j) / np.sqrt(2)) ** phase_power
    return phase * _native_matrix(reduced, free_width, columns)


class NativeResidualTableTests(unittest.TestCase):
    def test_phase_exact_toffoli_and_two_row_pauli_mask_kernels(self):
        # The enabling CCX must be literal, not a relative-phase Toffoli.
        for left, right, target in ((0, 1, 2), (2, 1, 0)):
            actual = _native_matrix(_toffoli(left, right, target), 3)
            expected = np.zeros((8, 8), dtype=complex)
            for index in range(8):
                flip = ((index >> left) & 1) * ((index >> right) & 1)
                expected[index ^ (flip << target), index] = 1
            np.testing.assert_allclose(actual, expected, atol=2e-14, rtol=0)
            self.assertEqual(sum(name in ("T", "TDG")
                                 for name, _ in _toffoli(left, right, target)), 7)

        programs = (
            program_rotation(_point(Fraction(3, 5)), _point(Fraction(4, 5)), 5),
            program_rotation(_point(Fraction(-5, 13)), _point(Fraction(12, 13)), 5),
        )
        address = 6
        gates = _table_mask(programs, address)
        self.assertFalse(any(name != "CX" and address in wires for name, wires in gates))
        self.assertFalse(any(name == "CX" and wires[1] == address for name, wires in gates))
        self.assertFalse(any(name in ("T", "TDG") for name, _ in gates))
        for value, program in enumerate(programs):
            actual = _sector_matrix(gates, 7, {address: value})
            literal = _native_matrix(_mask(program.signs), 6)
            np.testing.assert_allclose(actual, literal, atol=2e-14, rtol=0)
            np.testing.assert_allclose(
                _sector_matrix(_inverse(gates), 7, {address: value}),
                literal.conj().T, atol=2e-14, rtol=0)

    def test_complete_native_enabled_residual_sectors_and_shared_outer_table(self):
        coefficients = ((Fraction(-3, 4), Fraction(1, 2)),
                        (Fraction(0), Fraction(0)))
        table = emit_residual_table(coefficients, 5, enable_value=0)
        outer, middle, repeated = table.rotations
        self.assertIs(outer, repeated)
        self.assertIs(table.programs[0], table.programs[2])
        self.assertEqual(table.gates, outer.gates + middle.gates + outer.gates)
        self.assertEqual((table.address, table.enable, table.target), (8, 9, 7))
        identity = np.eye(256)
        matrices = {}
        for index, rotation in enumerate((outer, middle)):
            self.assertEqual(rotation.programs, table.programs[index])
            # Negative enable consists of complete-word X wrappers, not
            # a false assumption that enable is fixed to its input value.
            positive = emit_rotation_table(rotation.programs, "z" if index == 0 else "y")
            wrapper = (("X", (table.enable,)),)
            self.assertEqual(rotation.gates, wrapper + positive.gates + wrapper)
            for address in (0, 1):
                expected_active = _algebra_rotation(
                    rotation.programs[address],
                    _majorana_triplet(rotation.programs[address]),
                    "z" if index == 0 else "y")
                for enable in (0, 1):
                    with self.subTest(stage=index, address=address, enable=enable):
                        actual = _sector_matrix(
                            rotation.gates, table.enable + 1,
                            {table.address: address, table.enable: enable})
                        expected = expected_active if enable == 0 else identity
                        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                        matrices[index, address, enable] = actual

        # The literal block equalities above cover every input column and
        # relative phase across all four sectors, hence coherent controls
        # and arbitrary external references; no per-block phase is removed.
        for address in (0, 1):
            for enable in (0, 1):
                zword, yword = matrices[0, address, enable], matrices[1, address, enable]
                combined = zword @ yword @ zword
                np.testing.assert_allclose(combined.conj().T @ combined,
                                           identity, atol=3 * ATOL, rtol=0)
                if enable == 1:
                    np.testing.assert_allclose(combined, identity, atol=3 * ATOL, rtol=0)
                else:
                    # Signal-X symmetry is the all-input, borrowed-signal
                    # guarantee's hypothesis, not a zero-signal assumption.
                    flip = np.arange(256) ^ (1 << table.signal)
                    np.testing.assert_allclose(combined[flip], combined[:, flip],
                                               atol=3 * ATOL, rtol=0)
                    if address == 1:
                        # An enabled U(0) row is a pi/2 rotation, unlike
                        # the exact inactive identity on the same address.
                        self.assertGreater(np.linalg.norm((combined - identity)[:, 0]), 1)

        # Check a literal inverse on arbitrary target/signal/core columns.
        selected = [0, 1, 17, 63, 64, 127, 128, 255]
        active = matrices[0, 1, 0]
        inverse_columns = _sector_matrix(
            _inverse(outer.gates), table.enable + 1,
            {table.address: 1, table.enable: 0}, active[:, selected])
        np.testing.assert_allclose(inverse_columns, identity[:, selected],
                                   atol=2 * ATOL, rtol=0)

    def test_sector_restriction_retains_fixed_wire_phases_and_rejects_mixing(self):
        # Dropping a fixed enable's T phase or failing to track X wrappers
        # would turn this nontrivial literal phase into the identity.
        word = (("X", (1,)), ("T", (1,)), ("X", (1,)))
        np.testing.assert_allclose(_sector_matrix(word, 2, {1: 0}),
                                   ((1 + 1j) / np.sqrt(2)) * np.eye(2),
                                   atol=2e-15, rtol=0)
        np.testing.assert_allclose(_sector_matrix(word, 2, {1: 1}),
                                   np.eye(2), atol=2e-15, rtol=0)
        for bad in ((("H", (1,)),), (("CX", (0, 1)),), (("X", (1,)),)):
            with self.assertRaises(AssertionError):
                _restrict_sector(bad, 2, {1: 0})

    def test_high_precision_certificates_layout_and_linear_native_ledger(self):
        for q in (5, 16, 80, 129):
            coefficients = ((Fraction(-3, 4), Fraction(1, 2)),
                            (Fraction(0), Fraction(0)))
            table = emit_residual_table(coefficients, q)
            self.assertEqual(table.core, tuple(range(q + 1)))
            self.assertEqual((table.signal, table.target, table.address, table.enable),
                             (q + 1, q + 2, q + 3, q + 4))
            self.assertIs(table.rotations[0], table.rotations[2])
            self.assertIs(table.programs[0], table.programs[2])
            self.assertLess(table.operator_error_bound, Fraction(130, 1 << q))
            expected_error = max(data.completion_error_bound for data in table.rotation_data)
            expected_error += sum(max(program.operator_error_bound for program in programs)
                                  for programs in table.programs)
            self.assertEqual(table.operator_error_bound, expected_error)
            self.assertGreater(table.operator_error_bound, 0)
            for rotation in table.rotations:
                self.assertEqual(rotation.operator_error_bound,
                                 max(program.operator_error_bound for program in rotation.programs))
                self.assertEqual(sum(name in ("T", "TDG") for name, _ in rotation.gates),
                                 180 * q + 210)
                self.assertLessEqual(len(rotation.gates), 1720 * q)
            self.assertLessEqual(len(table.gates), 5160 * q)
            # Native T counts include every enabling Toffoli, rather than
            # treating conditioned sources as uncharged calls.
            native_t = sum(name in ("T", "TDG") for name, _ in table.gates)
            # Fifteen scalar words per rotation, each with two native
            # seven-T enabling Toffolis, add exactly 210 T gates.
            self.assertEqual(native_t, 540 * q + 630)
            for name, wires in table.gates:
                self.assertIn(name, ALPHABET)
                self.assertEqual(len(wires), 2 if name == "CX" else 1)
                self.assertEqual(len(wires), len(set(wires)))
                self.assertTrue(all(0 <= wire <= table.enable for wire in wires))
                self.assertFalse(name != "CX" and table.address in wires)
                if name == "CX":
                    self.assertNotIn(wires[1], (table.address, table.enable))

    def test_rejects_invalid_rows_predicates_axes_and_forged_programs(self):
        program = program_rotation(_point(1), _point(0), 5)
        other_q = program_rotation(_point(0), _point(1), 6)
        for rows in ((program,), (program, program, program), (program, other_q)):
            with self.assertRaises(ValueError):
                emit_rotation_table(rows)
        with self.assertRaises(ValueError):
            emit_rotation_table((program, program), "x")
        for value in (2, -1, True, False, 0.0, 1.0, Fraction(1)):
            with self.subTest(enable=value), self.assertRaises(ValueError):
                emit_rotation_table((program, program), enable_value=value)
            with self.subTest(enable=value), self.assertRaises(ValueError):
                emit_residual_table(((0, 0), (0, 0)), 5, enable_value=value)
        with self.assertRaises(ValueError):
            emit_rotation_table((program, replace(program, s=program.s + 1)))
        with self.assertRaises(ValueError):
            emit_residual_table(((0, 0),), 5)
        with self.assertRaises(ValueError):
            emit_residual_table(((0, 0), (1, 1)), 5)
        for q in (True, False, 5.0, Fraction(5), 4):
            with self.subTest(q=q), self.assertRaises(ValueError):
                emit_residual_table(((0, 0), (0, 0)), q)
        for value in (0.5, True, False):
            with self.subTest(coefficient=value), self.assertRaises(TypeError):
                emit_residual_table(((value, 0), (0, 0)), 5)
        with self.assertRaises(ValueError):
            emit_residual_table(((Fraction(1, 3), 0), (0, 0)), 5)


if __name__ == "__main__":
    unittest.main()
